"""BFQ input selection and scheduling around gcf-tools domain validation."""

from __future__ import annotations

import csv
import hashlib
import json
import logging
import os
import re
import shutil
import uuid

from dataclasses import dataclass, replace
from pathlib import Path

from bcl2fastq_pipeline import initial_holds, project_scheduling
from bcl2fastq_pipeline.config import parse_custom_options
from bcl2fastq_pipeline.state import StateError, utcnow

log = logging.getLogger(__name__)


@dataclass(frozen=True)
class InputSelection:
    sample_sheet: Path
    submission_form: Path
    pairs: tuple = ()
    sample_sheets: tuple = ()
    submission_forms: tuple = ()
    multi_pair: bool = False
    ignored_duplicates: tuple = ()
    run_info: Path | None = None

    def describe(self) -> str:
        if self.multi_pair:
            lines = ["Project-paired flowcell inputs:"]
            for pair in self.pairs:
                lines.append(
                    f"{pair.project_id}: SampleSheet: {pair.sample_sheet}\n"
                    f"{pair.project_id}: Submission form: {pair.submission_form}"
                )
            if not self.pairs:
                lines.extend(str(path) for path in self.sample_sheets + self.submission_forms)
            return "\n".join(lines)
        return f"SampleSheet: {self.sample_sheet}\nSubmission form: {self.submission_form}"


_PROJECT = re.compile(r"(?<![A-Za-z0-9])GCF-\d{4}-\d{3}(?![A-Za-z0-9])", re.IGNORECASE)


def _project_key(path):
    keys = _PROJECT.findall(path.name)
    return keys[0].upper() if len(keys) == 1 else None


def _input_files(directory):
    return (
        tuple(sorted(directory.glob("SampleSheet*.csv"))),
        tuple(sorted(directory.glob("*Sample-Submission-Form*.xlsx"))),
    )


def _multi_mode(sheets, forms):
    """Keep historical single-pair selection; activate pairing for multiple sheets."""
    hashes = set()
    keys = set()
    for path in sheets:
        key = _project_key(path)
        if key:
            keys.add(key)
        try:
            hashes.add(hashlib.sha256(path.read_bytes()).hexdigest())
        except OSError:
            hashes.add(str(path))
    if not keys and len(forms) <= 1:
        # Historical instruments may include an unrelated generic sheet next to
        # the BFQ opted-in sheet. Keep that selection contract for one form.
        return False
    return (
        len(keys) > 1
        or len(hashes) > 1
        or (not sheets and len({_project_key(path) for path in forms} - {None}) > 1)
    )


def _project_selection(sheets, forms):
    from configmaker.project_pairs import resolve_project_pairs  # noqa: PLC0415

    pairing = resolve_project_pairs(sheets, forms)
    # Keep candidates even on failure: aggregate diagnostics must include invalid
    # and unmatched files rather than validating only successfully resolved pairs.
    return InputSelection(
        sheets[0] if sheets else Path("SampleSheet.csv"),
        forms[0] if forms else Path("Sample-Submission-Form.xlsx"),
        tuple(
            replace(
                pair,
                sample_sheet=Path(pair.sample_sheet),
                submission_form=Path(pair.submission_form),
            )
            for pair in pairing.pairs
        ),
        tuple(sheets),
        tuple(forms),
        True,
        tuple(pairing.duplicates),
    )


class PreflightValidationError(StateError):
    """A structured input failure, also consumable by BFQ's error reporter."""

    def __init__(self, validation_result, report_path=None):
        self.validation_result = validation_result
        self.report_path = Path(report_path) if report_path else None
        detail = validation_result.render_text()
        if self.report_path:
            detail += f"\nComplete structured report: {self.report_path}"
        super().__init__(f"Input preflight failed.\n{detail}")


def has_custom_options_marker(path: Path) -> bool:
    """Recognize explicit BFQ opt-in even when unrelated bytes break CSV decoding."""
    try:
        content = path.read_bytes()
    except OSError:
        return False
    for line in content.splitlines():
        first = line.split(b",", 1)[0].removeprefix(b"\xef\xbb\xbf").strip().strip(b'"')
        if first.lower() == b"[customoptions]":
            return True
    return False


def _candidate(directory: Path, canonical: str, pattern: str, *, sheet=False, curated=False):
    canonical_path = directory / canonical
    # A curated canonical file always wins, including when unreadable or malformed.
    # A broken symlink must likewise be reported, not silently replaced from source.
    if curated and (canonical_path.exists() or canonical_path.is_symlink()):
        return canonical_path
    candidates = sorted(directory.glob(pattern))
    if sheet:
        malformed_opt_in = None
        for path in candidates:
            try:
                options, _ = parse_custom_options(path)
            except (OSError, UnicodeError, csv.Error):
                if malformed_opt_in is None and has_custom_options_marker(path):
                    malformed_opt_in = path
                continue
            if options:
                return path
        if malformed_opt_in is not None:
            return malformed_opt_in
    if canonical_path.exists() or canonical_path.is_symlink():
        return canonical_path
    return candidates[0] if candidates else None


def select_run_inputs(source_path, output_path, *, refresh=False, projects=None) -> InputSelection:
    """Resolve each effective input independently without modifying any file.

    Output-side inputs are curated. Fill only missing inputs from the instrument;
    explicit refresh chooses both instrument inputs. Missing paths are returned so
    the shared validator can explain the problem in its normal diagnostic format.
    """
    source = Path(source_path)
    output = Path(output_path)
    source_files = _input_files(source)
    output_files = _input_files(output)
    if (
        _multi_mode(*output_files)
        or (refresh and _multi_mode(*source_files))
        or (not (output / "SampleSheet.csv").exists() and _multi_mode(*source_files))
    ):
        chosen = []
        for original, curated in zip(source_files, output_files, strict=True):
            if refresh and projects is None:
                chosen.append(original)
                continue
            selected = set(projects or ()) if refresh else set()
            keys = {_project_key(path) for path in original + curated} - {None}
            effective = []
            for key in sorted(keys):
                source_candidates = [path for path in original if _project_key(path) == key]
                curated_candidates = [path for path in curated if _project_key(path) == key]
                effective.extend(
                    source_candidates
                    if key in selected
                    else curated_candidates or source_candidates
                )
            # Generic files take part in duplicate/conflict diagnostics; a curated
            # generic sheet never silently replaces a project-qualified input.
            effective.extend(path for path in curated or original if _project_key(path) is None)
            chosen.append(tuple(effective))
        selection = _project_selection(*chosen)
        geometry = next(
            (
                directory / "RunInfo.xml"
                for directory in (output, source)
                if (directory / "RunInfo.xml").is_file()
            ),
            None,
        )
        return replace(selection, run_info=geometry)
    chosen = []
    for canonical, pattern, sheet in (
        ("SampleSheet.csv", "SampleSheet*.csv", True),
        ("Sample-Submission-Form.xlsx", "*Sample-Submission-Form*.xlsx", False),
    ):
        path = (
            None if refresh else _candidate(output, canonical, pattern, sheet=sheet, curated=True)
        )
        if path is None:
            path = _candidate(source, canonical, pattern, sheet=sheet)
        chosen.append(path if path is not None else source / canonical)
    return InputSelection(*chosen)


def copy_run_inputs(selection: InputSelection, output_path) -> InputSelection:
    """Copy selected inputs to canonical downstream paths without changing content."""
    output = Path(output_path)
    output.mkdir(parents=True, exist_ok=True)
    if selection.multi_pair:
        result = validate_selection(selection)
        if not result.ok:
            raise PreflightValidationError(result)
        effective_names = {
            path.name
            for pair in selection.pairs
            for path in (pair.sample_sheet, pair.submission_form)
        }
        replacing = any(
            path.parent != output
            for pair in selection.pairs
            for path in (pair.sample_sheet, pair.submission_form)
        )
        obsolete = [
            path
            for paths in _input_files(output)
            for path in paths
            if replacing and path.name not in effective_names
        ]
        for pair in selection.pairs:
            for source in (pair.sample_sheet, pair.submission_form):
                destination = output / source.name
                if source != destination:
                    temporary = output / f".{destination.name}.{os.getpid()}.tmp"
                    try:
                        shutil.copy2(source, temporary)
                        os.replace(temporary, destination)
                    finally:
                        temporary.unlink(missing_ok=True)
        if selection.run_info is not None:
            destination = output / "RunInfo.xml"
            if selection.run_info != destination:
                shutil.copy2(selection.run_info, destination)
        if obsolete:
            # Explicit replacement must not leave obsolete canonical inputs as
            # competing candidates. Preserve their bytes outside discovery.
            retained = output / "provenance" / "input-replacements" / uuid.uuid4().hex
            retained.mkdir(parents=True)
            for path in obsolete:
                os.replace(path, retained / path.name)
        sheets = tuple(output / pair.sample_sheet.name for pair in selection.pairs)
        forms = tuple(output / pair.submission_form.name for pair in selection.pairs)
        copied = replace(
            _project_selection(sheets, forms),
            ignored_duplicates=selection.ignored_duplicates,
            run_info=(output / "RunInfo.xml") if selection.run_info else None,
        )
        provenance = result.to_dict()
        manifest = output / "project-inputs.json"
        if not manifest.exists() or any(
            path.parent != output
            for pair in selection.pairs
            for path in (pair.sample_sheet, pair.submission_form)
        ):
            manifest.write_text(
                json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8"
            )
        return copied
    destinations = []
    for source, name in (
        (selection.sample_sheet, "SampleSheet.csv"),
        (selection.submission_form, "Sample-Submission-Form.xlsx"),
    ):
        destination = output / name
        if source != destination and source.exists():
            temporary = output / f".{name}.{os.getpid()}.tmp"
            try:
                # Replacing the canonical file avoids partial copies and never
                # overwrites a source reached through an old output-side symlink.
                shutil.copy2(source, temporary)
                os.replace(temporary, destination)
            finally:
                temporary.unlink(missing_ok=True)
        destinations.append(destination)
    return InputSelection(*destinations)


def validate_selection(
    selection: InputSelection, *, run_info=None, mismatch_policy=None, libprep_config=None
):
    """Always validate current bytes; successful reports are never used as a cache."""
    # flowcell-manager's entry point removes its executable directory from sys.path
    # before domain imports, avoiding a neighbouring configmaker.py script shadowing
    # the package. Metadata-independent commands (including help) stay lightweight.
    from configmaker.validation import validate_inputs  # noqa: PLC0415

    if selection.multi_pair:
        from configmaker.project_pairs import validate_project_pairs  # noqa: PLC0415

        from bcl2fastq_pipeline.demultiplexing import (  # noqa: PLC0415
            DemultiplexingPlanError,
            plan,
        )

        run_info = run_info or selection.run_info
        result = validate_project_pairs(
            selection.sample_sheets,
            selection.submission_forms,
            run_info=run_info,
            mismatch_policy=mismatch_policy,
            libprep_config=libprep_config,
        )
        for duplicate in selection.ignored_duplicates:
            if duplicate not in result.duplicates:
                result.duplicates.append(duplicate)
        if result.ok:
            try:
                plan(
                    selection,
                    run_info=run_info,
                    force_bcl2fastq=bool(os.environ.get("FORCE_BCL2FASTQ")),
                )
            except DemultiplexingPlanError as error:
                result.add("error", "demux.unsupported_plan", str(error))
        return initial_holds.add_preflight_plan(
            selection, project_scheduling.add_preflight_plan(selection, result)
        )
    result = validate_inputs([selection.sample_sheet], [selection.submission_form])
    return initial_holds.add_preflight_plan(
        selection, project_scheduling.add_preflight_plan(selection, result)
    )


def require_valid_inputs(  # noqa: PLR0913
    source_path, output_path, *, refresh=False, projects=None, show_schedule=True, state=None
):
    """Read-only preparation guard before an operator command removes results."""
    selection = select_run_inputs(source_path, output_path, refresh=refresh, projects=projects)
    run_info = Path(source_path) / "RunInfo.xml"
    libprep_config = None
    if selection.multi_pair:
        from bcl2fastq_pipeline.workflow_config import AUTHORITATIVE_CONFIG  # noqa: PLC0415

        if AUTHORITATIVE_CONFIG.is_file():
            libprep_config = AUTHORITATIVE_CONFIG
    result = (
        validate_selection(
            selection,
            run_info=run_info if run_info.exists() else None,
            libprep_config=libprep_config,
        )
        if selection.multi_pair
        else validate_selection(selection)
    )
    print(selection.describe())
    print(result.render_text())
    if not result.ok:
        raise PreflightValidationError(result)
    if show_schedule:
        project_scheduling.print_plan(result.summary.get("project_schedule"), projects=projects)
    initial_holds.print_plan(
        result.summary.get("project_initial_holds"), state=state, projects=projects
    )
    return selection, result


def freeze_execution_inputs(cfg, store):
    """Bind analysis to attempt-owned copies of validated bytes, including legacy pairs."""
    selection = getattr(cfg.run, "input_selection", None)
    if not selection or getattr(cfg.run, "selected_projects", None):
        return
    state = store.read(cfg.run.run_id)
    directory = (
        cfg.static.paths.manager_dir / "execution-inputs" / cfg.run.run_id / str(state["attempt"])
    )
    if not selection.multi_pair:
        # Legacy demultiplexers read this path inside their tool container; keep
        # it under the already mounted output tree, rather than the manager mount.
        directory = cfg.output_path / "provenance" / "execution-inputs" / str(state["attempt"])
        if directory.parent.is_symlink() or directory.parent.parent.is_symlink():
            raise StateError("Execution-input provenance directories must not be symlinks")
    copied = copy_run_inputs(selection, directory)
    geometry = []
    if cfg.run.flowcell_path is not None:
        source_info = cfg.run.flowcell_path / "RunInfo.xml"
        if source_info.is_file():
            content = source_info.read_bytes()
            (directory / "RunInfo.xml").write_bytes(content)
            if source_info.read_bytes() != content:
                raise StateError("Physical run geometry changed while freezing execution inputs")
            geometry.append(
                {
                    "path": str(directory / "RunInfo.xml"),
                    "sha256": hashlib.sha256(content).hexdigest(),
                }
            )
    before = validate_selection(selection).to_dict()["inputs"]
    after = validate_selection(copied).to_dict()["inputs"]
    identities = lambda items: sorted((Path(item["path"]).name, item["sha256"]) for item in items)
    # Duplicate candidates are recorded separately in the original report; compare
    # exactly the effective pair files used during execution.
    effective_names = {
        path.name
        for pair in (selection.pairs or (selection,))
        for path in (pair.sample_sheet, pair.submission_form)
    }
    if identities(
        [item for item in before if Path(item["path"]).name in effective_names]
    ) != identities([item for item in after if Path(item["path"]).name in effective_names]):
        raise StateError("Project inputs changed while freezing execution inputs; retry validation")
    cfg.run.input_selection = copied
    cfg.run.sample_sheet = copied.sample_sheet
    cfg.run.sample_submission_form = copied.submission_form

    def record(current):
        current["attempts"][-1]["execution_inputs"] = {
            "directory": str(directory),
            "inputs": after,
            "geometry": geometry,
        }
        return current

    store.mutate(cfg.run.run_id, record)


def run_preflight(cfg, store, stage):
    """Validate under the caller's execution lease and persist the complete report."""
    selection = getattr(cfg.run, "input_selection", None) or InputSelection(
        cfg.run.sample_sheet or cfg.output_path / "SampleSheet.csv",
        cfg.run.sample_submission_form or cfg.output_path / "Sample-Submission-Form.xlsx",
    )
    log.info("Input preflight before %s for %s\n%s", stage, cfg.run.run_id, selection.describe())
    run_info = cfg.run.flowcell_path / "RunInfo.xml" if cfg.run.flowcell_path else None
    frozen_info = selection.sample_sheet.parent / "RunInfo.xml"
    if frozen_info.is_file():
        run_info = frozen_info
    if run_info is not None and not run_info.exists():
        run_info = None
    # Preserve legacy callers and validator test doubles, while multi-pair mode
    # additionally validates physical run geometry across all projects.
    result = (
        validate_selection(selection, run_info=run_info)
        if selection.multi_pair
        else validate_selection(selection)
    )
    report_path = cfg.static.paths.report_dir / f"{cfg.run.run_id}.input-preflight.json"
    report = result.to_dict()
    report.update(checked_at=utcnow(), restart_boundary=stage, report_path=str(report_path))
    # State is durable independently of the output directory, including failed input
    # reports, and keeps the validator version and hashes of exactly the parsed bytes.
    store.record_input_preflight(cfg.run.run_id, stage, report)
    temporary = report_path.with_name(f".{report_path.name}.{os.getpid()}.tmp")
    saved_report_path = None
    try:
        report_path.parent.mkdir(parents=True, exist_ok=True)
        temporary.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        os.replace(temporary, report_path)
        saved_report_path = report_path
    except OSError as error:
        # The durable state copy is sufficient; never replace an actionable input
        # diagnosis with a secondary sidecar failure or link to an older report.
        log.warning(
            "Cannot write input-preflight sidecar %s: %s; full report retained in state",
            report_path,
            error,
        )
        report.update(report_path=None, report_write_error=str(error))
        store.record_input_preflight(cfg.run.run_id, stage, report)
    finally:
        try:
            temporary.unlink(missing_ok=True)
        except OSError:
            log.warning("Cannot remove temporary input-preflight report %s", temporary)
    log.info("%s", result.render_text())
    if not result.ok:
        raise PreflightValidationError(result, saved_report_path)
    return result
