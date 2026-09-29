"""BFQ input selection and scheduling around gcf-tools domain validation."""

from __future__ import annotations

import csv
import json
import logging
import os
import shutil

from dataclasses import dataclass
from pathlib import Path

from bcl2fastq_pipeline.config import parse_custom_options
from bcl2fastq_pipeline.state import StateError, utcnow

log = logging.getLogger(__name__)


@dataclass(frozen=True)
class InputSelection:
    sample_sheet: Path
    submission_form: Path

    def describe(self) -> str:
        return f"SampleSheet: {self.sample_sheet}\nSubmission form: {self.submission_form}"


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


def select_run_inputs(source_path, output_path, *, refresh=False) -> InputSelection:
    """Resolve each effective input independently without modifying any file.

    Output-side inputs are curated. Fill only missing inputs from the instrument;
    explicit refresh chooses both instrument inputs. Missing paths are returned so
    the shared validator can explain the problem in its normal diagnostic format.
    """
    source = Path(source_path)
    output = Path(output_path)
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


def validate_selection(selection: InputSelection):
    """Always validate current bytes; successful reports are never used as a cache."""
    # flowcell-manager's entry point removes its executable directory from sys.path
    # before domain imports, avoiding a neighbouring configmaker.py script shadowing
    # the package. Metadata-independent commands (including help) stay lightweight.
    from configmaker.validation import validate_inputs  # noqa: PLC0415

    return validate_inputs([selection.sample_sheet], [selection.submission_form])


def require_valid_inputs(source_path, output_path, *, refresh=False):
    """Read-only preparation guard before an operator command removes results."""
    selection = select_run_inputs(source_path, output_path, refresh=refresh)
    result = validate_selection(selection)
    print(selection.describe())
    print(result.render_text())
    if not result.ok:
        raise PreflightValidationError(result)
    return selection, result


def run_preflight(cfg, store, stage):
    """Validate under the caller's execution lease and persist the complete report."""
    selection = InputSelection(
        cfg.run.sample_sheet or cfg.output_path / "SampleSheet.csv",
        cfg.run.sample_submission_form or cfg.output_path / "Sample-Submission-Form.xlsx",
    )
    log.info("Input preflight before %s for %s\n%s", stage, cfg.run.run_id, selection.describe())
    result = validate_selection(selection)
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
