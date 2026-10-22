"""Structured, reproducible conversion inputs for project-specific input pairs.

BCL Convert >=4.1 supports the seven per-sample settings below on v2 sheets.
Other disagreements fail before conversion; they are never merged last-file-wins.
"""

import csv
import hashlib
import io
import json
import re
import xml.etree.ElementTree as ET

from pathlib import Path
from types import SimpleNamespace

from bcl2fastq_pipeline import tenx_demultiplexing
from bcl2fastq_pipeline.config import parse_custom_options

PER_SAMPLE_SETTINGS = {
    "OverrideCycles",
    "BarcodeMismatchesIndex1",
    "BarcodeMismatchesIndex2",
    "AdapterRead1",
    "AdapterRead2",
    "AdapterBehavior",
    "AdapterStringency",
}
GLOBAL_SETTINGS = PER_SAMPLE_SETTINGS | {
    "MinimumTrimmedReadLength",
    "MaskShortReads",
    "MinimumAdapterOverlap",
    "TrimUMI",
    "CreateFastqForIndexReads",
    "NoLaneSplitting",
    "FastqCompressionFormat",
    "FindAdaptersWithIndels",
    "IndependentIndexCollisionCheck",
    "SoftwareVersion",
}
DATA_COLUMNS = {"Lane", "Sample_ID", "Sample_Name", "Sample_Project", "index", "index2"}
ANNOTATION_COLUMNS = {
    "Sample_Plate",
    "Sample_Well",
    "I7_Index_ID",
    "I5_Index_ID",
    "Index_Plate_Well",
    "Description",
    "Sample_Group",
    "Sample_Ref",
    "Sample_Control",
}
CANONICAL = {key.casefold(): key for key in GLOBAL_SETTINGS | DATA_COLUMNS}
DEFAULTS = {
    "BarcodeMismatchesIndex1": "1",
    "BarcodeMismatchesIndex2": "1",
    "AdapterRead1": "na",
    "AdapterRead2": "na",
    "AdapterBehavior": "trim",
    "AdapterStringency": "0.9",
}


class DemultiplexingPlanError(ValueError):
    """Effective inputs cannot be represented safely by the supported backend."""


def _sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _sheet(path):
    sections = {}
    section = None
    with Path(path).open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle)
        for cells in reader:
            if not any(cell.strip() for cell in cells):
                continue
            tag = cells[0].strip()
            if (
                tag.startswith("[")
                and tag.endswith("]")
                and all(not cell.strip() for cell in cells[1:])
            ):
                section = tag[1:-1].casefold()
                if section in sections:
                    raise DemultiplexingPlanError(f"{path}: repeated [{section}] section")
                sections[section] = []
            elif section:
                sections[section].append((reader.line_num, cells))
    data = sections.get("bclconvert_data", sections.get("data", []))
    if not data:
        raise DemultiplexingPlanError(f"{path}: no demultiplexing data rows")
    header = [CANONICAL.get(value.strip().casefold(), value.strip()) for value in data[0][1]]
    # Excel trailing commas are not named columns.
    named = [(i, name) for i, name in enumerate(header) if name]
    if len({name for _i, name in named}) != len(named):
        raise DemultiplexingPlanError(f"{path}: duplicate data column names")
    rows = [
        (number, {name: cells[i] if i < len(cells) else "" for i, name in named})
        for number, cells in data[1:]
    ]
    for number, row in rows:
        unsupported = {
            key
            for key, value in row.items()
            if key not in DATA_COLUMNS | PER_SAMPLE_SETTINGS | ANNOTATION_COLUMNS and value.strip()
        }
        if unsupported:
            raise DemultiplexingPlanError(
                f"{path}:{number}: unsupported populated demultiplexer data columns {sorted(unsupported)}. "
                "Use supported per-sample settings; global conversion settings belong in [Settings] "
                "or [BCLConvert_Settings]. Keep other metadata in the submission form."
            )
    settings = {}
    for name in ("settings", "bclconvert_settings"):
        for number, cells in sections.get(name, []):
            key = CANONICAL.get(cells[0].strip().casefold(), cells[0].strip())
            value = cells[1].strip() if len(cells) > 1 else ""
            if key in settings and settings[key] != value:
                raise DemultiplexingPlanError(f"{path}:{number}: conflicting setting {key}")
            settings[key] = value
    # These are documented BCL Convert equivalents of legacy adapter settings.
    for old in ("Adapter", "TrimAdapter"):
        if old in settings:
            adapter = settings.pop(old)
            for key in ("AdapterRead1", "AdapterRead2"):
                if key in settings and settings[key] != adapter:
                    raise DemultiplexingPlanError(f"{path}: conflicting {old}/{key} settings")
                settings[key] = adapter
    return settings, rows


def _physical_mask(run_info):
    if not run_info or not Path(run_info).is_file():
        return None
    reads = ET.parse(run_info).findall(".//Reads/Read")
    return ";".join(
        ("I" if read.attrib.get("IsIndexedRead") == "Y" else "Y") + read.attrib["NumCycles"]
        for read in reads
    )


def _combined_plan(selection, *, run_info=None, force_bcl2fastq=False):
    """Return one combined pass, or explain why the effective settings conflict."""
    pairs = sorted(selection.pairs, key=lambda pair: pair.project_id)
    if not pairs:
        return None
    if force_bcl2fastq:
        raise DemultiplexingPlanError(
            "FORCE_BCL2FASTQ is unsupported for multiple project input pairs: its command-line "
            "mismatch and cycle policies are not represented by this planner. Unset FORCE_BCL2FASTQ "
            "and use the supported BCL Convert backend. Legacy single-pair conversion is unchanged."
        )
    projects = {}
    sources = []
    all_rows = []
    for pair in pairs:
        custom, _ = parse_custom_options(Path(pair.sample_sheet))
        settings, rows = _sheet(pair.sample_sheet)
        if getattr(pair, "legacy_sheet", False):
            rows = [
                (number, row)
                for number, row in rows
                if row.get("Sample_Project") == pair.project_id
            ]
        try:
            profile, settings, rows = tenx_demultiplexing.configure(
                custom, settings, rows, run_info
            )
        except tenx_demultiplexing.ProfileError as error:
            raise DemultiplexingPlanError(str(error)) from error
        if profile["backend"] != "bcl-convert":
            if len({str(member.sample_sheet) for member in pairs}) > 1:
                raise DemultiplexingPlanError(
                    "Legacy backends require an explicit separate DemultiplexingPass per project"
                )
            legacy = _legacy_plan(pair, custom, profile)
            legacy["projects"] = sorted(member.project_id for member in pairs)
            legacy["sources"] = [
                {**legacy["sources"][0], "project": member.project_id} for member in pairs
            ]
            return legacy
        unsupported = set(settings) - GLOBAL_SETTINGS
        if unsupported and not force_bcl2fastq:
            raise DemultiplexingPlanError(
                f"{pair.project_id}: unsupported BCL Convert settings {sorted(unsupported)} "
                f"in {pair.sample_sheet}; replace obsolete settings explicitly before retrying."
            )
        projects[pair.project_id] = settings
        sources.append(
            {
                "project": pair.project_id,
                "sample_sheet": str(pair.sample_sheet),
                "sample_sheet_sha256": _sha256(pair.sample_sheet),
                "submission_form": str(pair.submission_form),
                "submission_form_sha256": _sha256(pair.submission_form),
                "custom": custom,
                "settings": settings,
                "profile": profile,
            }
        )
        for number, row in rows:
            if row.get("Sample_Project") not in (None, "", pair.project_id):
                raise DemultiplexingPlanError(
                    f"{pair.sample_sheet}:{number}: Sample_Project does not match {pair.project_id}"
                )
            row["Sample_Project"] = pair.project_id
            all_rows.append((pair.project_id, str(pair.sample_sheet), number, row))
    keys = set().union(*(settings for settings in projects.values()))
    per_sample = set().union(*(set(row) & PER_SAMPLE_SETTINGS for *_source, row in all_rows))
    if force_bcl2fastq and per_sample:
        raise DemultiplexingPlanError(
            "FORCE_BCL2FASTQ cannot represent per-sample settings; use BCL Convert >=4.1."
        )
    globals_ = {}
    for key in keys:
        values = {settings.get(key, DEFAULTS.get(key, "")) for settings in projects.values()}
        if len(values) == 1 and key not in per_sample:
            globals_[key] = values.pop()
        elif key in PER_SAMPLE_SETTINGS and not force_bcl2fastq:
            per_sample.add(key)
        else:
            details = ", ".join(
                f"{p}={settings.get(key, '<default>')!r}" for p, settings in projects.items()
            )
            raise DemultiplexingPlanError(
                f"Cannot combine projects with different {key} settings: {details}. "
                "Use matching settings or an explicitly supported separate-pass backend plan."
            )
    if force_bcl2fastq and per_sample:
        raise DemultiplexingPlanError(
            "FORCE_BCL2FASTQ cannot represent per-sample conversion settings; use BCL Convert >=4.1."
        )
    version_two = bool(per_sample) and not force_bcl2fastq
    columns = [
        key
        for key in ("Lane", "Sample_ID", "Sample_Name", "Sample_Project", "index", "index2")
        if any(key in row for *_source, row in all_rows)
        # v2 rejects Sample_Name unless name-based output is explicitly enabled.
        # BFQ names reads by Sample_ID; keep the optional name in source provenance.
        and not (version_two and key == "Sample_Name")
    ]
    columns += sorted(per_sample)
    buffer = io.StringIO(newline="")
    writer = csv.writer(buffer, lineterminator="\n")
    writer.writerows(
        [
            ["[Header]"],
            [
                "FileFormatVersion" if version_two else "IEMFileVersion",
                "2" if version_two else "4",
            ],
            [],
            ["[BCLConvert_Settings]" if version_two else "[Settings]"],
        ]
    )
    writer.writerows(sorted(globals_.items()))
    writer.writerows([[], ["[BCLConvert_Data]" if version_two else "[Data]"], columns])
    mappings = []
    for project, source, source_row, row in all_rows:
        effective = dict(row)
        for key in per_sample:
            value = row.get(key) or projects[project].get(key) or DEFAULTS.get(key)
            if key == "OverrideCycles" and not value:
                value = _physical_mask(run_info)
            if value is None:
                raise DemultiplexingPlanError(
                    f"{project}: {key} needs an explicit value or a valid RunInfo.xml before combining."
                )
            if key.startswith("BarcodeMismatchesIndex"):
                index = "index" if key.endswith("1") else "index2"
                if not row.get(index) or row.get(index, "").casefold() == "na":
                    value = "na"
            effective[key] = value
        writer.writerow([effective.get(key, "") for key in columns])
        mappings.append(
            {
                "derived_row": buffer.getvalue().count("\n"),
                "project": project,
                "source": source,
                "source_row": source_row,
                "sample_id": row.get("Sample_ID"),
                "effective_parameters": {
                    key: value for key, value in effective.items() if not key.startswith("_")
                },
                "index_provenance": row.get("_index_provenance"),
            }
        )
    content = buffer.getvalue()
    return {
        "schema_version": 1,
        "tool": "bcl2fastq" if force_bcl2fastq else "bcl-convert",
        "passes": 1,
        "sample_sheet_version": 2 if version_two else 1,
        "projects": sorted(projects),
        "sources": sources,
        "global_settings": globals_,
        "per_sample_settings": sorted(per_sample),
        "rows": mappings,
        "sample_sheet_content": content,
        "sample_sheet_sha256": hashlib.sha256(content.encode()).hexdigest(),
    }


def conversion_pairs(selection):
    if getattr(selection, "pairs", ()):
        return selection.pairs
    if not getattr(selection, "sample_sheet", None):
        return ()
    _settings, rows = _sheet(selection.sample_sheet)
    projects = sorted({row.get("Sample_Project") for _number, row in rows} - {None, ""})
    return tuple(
        SimpleNamespace(
            project_id=project,
            sample_sheet=selection.sample_sheet,
            submission_form=selection.submission_form,
            legacy_sheet=True,
        )
        for project in projects
    )


def needs_plan(selection):
    if getattr(selection, "pairs", ()):
        return True
    custom, _ = parse_custom_options(Path(selection.sample_sheet))
    return "10X Genomics" in custom.get("Libprep", "") or bool(custom.get("Demultiplexer"))


def _legacy_plan(pair, custom, profile):
    content = Path(pair.sample_sheet).read_text(encoding="utf-8-sig")
    return {
        "schema_version": 1,
        "tool": profile.get("legacy_tool", profile["backend"]),
        "backend": profile["backend"],
        "passes": 1,
        "sample_sheet_version": 1,
        "projects": [pair.project_id],
        "rows": [],
        "global_settings": {},
        "per_sample_settings": [],
        "sources": [
            {
                "project": pair.project_id,
                "sample_sheet": str(pair.sample_sheet),
                "sample_sheet_sha256": _sha256(pair.sample_sheet),
                "submission_form": str(pair.submission_form),
                "submission_form_sha256": _sha256(pair.submission_form),
                "custom": custom,
                "profile": profile,
                "settings": _sheet(pair.sample_sheet)[0],
            }
        ],
        "sample_sheet_content": content,
        "sample_sheet_sha256": hashlib.sha256(content.encode()).hexdigest(),
    }


def plan(selection, *, run_info=None, force_bcl2fastq=False):
    """Combine compatible projects; separate passes require explicit source labels."""
    pairs = conversion_pairs(selection)
    if not pairs:
        return None
    grouped = {}
    for pair in pairs:
        custom, _ = parse_custom_options(Path(pair.sample_sheet))
        label = custom.get("DemultiplexingPass", "combined")
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", label):
            raise DemultiplexingPlanError(
                "DemultiplexingPass must be an alphanumeric label with optional dash/underscore"
            )
        grouped.setdefault(label, []).append(pair)
    plans = []
    for label, members in sorted(grouped.items()):
        item = _combined_plan(
            SimpleNamespace(pairs=members), run_info=run_info, force_bcl2fastq=force_bcl2fastq
        )
        item["pass_id"] = label
        if run_info:
            item["run_info"] = {"path": str(run_info), "sha256": _sha256(run_info)}
        plans.append(item)
    if len(plans) == 1:
        return plans[0]
    return {
        "schema_version": 1,
        "tool": "multiple",
        "passes": len(plans),
        "pass_plans": plans,
        "projects": sorted(p.project_id for p in pairs),
        "sources": [source for item in plans for source in item["sources"]],
    }


def pass_plans(result):
    return result.get("pass_plans", [result])


def prepare(cfg):
    """Freeze the derived sheet and its provenance inside this attempt's output."""
    import os  # noqa: PLC0415

    from bcl2fastq_pipeline.workflow_config import physical_run_info  # noqa: PLC0415

    selection = cfg.run.input_selection
    result = plan(
        selection,
        run_info=physical_run_info(cfg),
        force_bcl2fastq=bool(os.environ.get("FORCE_BCL2FASTQ")),
    )
    destination = cfg.output_path / "provenance" / "demultiplexing"
    destination.mkdir(parents=True, exist_ok=True)
    for item in pass_plans(result):
        sheet = destination / f"SampleSheet.{item['pass_id']}.csv"
        sheet.write_text(item["sample_sheet_content"], encoding="utf-8")
        item["derived_sample_sheet"] = str(sheet)
    for pair in conversion_pairs(selection):
        source = destination / "sources" / pair.project_id
        source.mkdir(parents=True, exist_ok=True)
        for path in (Path(pair.sample_sheet), Path(pair.submission_form)):
            (source / path.name).write_bytes(path.read_bytes())
    for item in pass_plans(result):
        for source in item["sources"]:
            profile = source.get("profile", {})
            if profile.get("source"):
                content = Path(profile["source"]).read_bytes()
                if hashlib.sha256(content).hexdigest() != profile["sha256"]:
                    raise DemultiplexingPlanError(
                        "Demultiplexing profile changed while freezing conversion inputs"
                    )
                (destination / "demultiplexing.config").write_bytes(content)
    sheet = Path(pass_plans(result)[0]["derived_sample_sheet"])
    cfg.run.demultiplexing_plan = result
    cfg.run.sample_sheet = sheet
    record(cfg)
    return sheet


def record(cfg, *, version=None, command=None, status="prepared", item=None, image=None):  # noqa: PLR0913
    result = cfg.run.demultiplexing_plan
    target = item if item is not None else result
    target.update({"status": status})
    if version is not None:
        target["tool_version"] = version
        target["tool_image"] = image or version
    if command is not None:
        target["command"] = command
    destination = cfg.output_path / "provenance" / "demultiplexing" / "plan.json"
    destination.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def delivery_provenance(cfg, project, *, refresh=False):
    """Keep retained FASTQs tied to their own conversion on downstream retries."""
    destination = cfg.output_path / "provenance" / "project-demultiplexing" / project
    canonical = cfg.static.paths.output_dir / cfg.run.run_id
    if not refresh:
        for retained in (
            destination,
            canonical / "provenance" / "project-demultiplexing" / project,
        ):
            manifest = retained / "plan.json"
            if not manifest.exists():
                continue
            try:
                data = json.loads(manifest.read_text())
                if (
                    data["project"] != project
                    or _sha256(retained / "SampleSheet.project.csv")
                    != data["project_sample_sheet_sha256"]
                ):
                    raise ValueError("project identity or derived sheet hash differs")
            except (OSError, ValueError, KeyError, TypeError) as error:
                raise DemultiplexingPlanError(
                    f"Cannot recover retained demultiplexing provenance from {manifest}: {error}"
                ) from error
            return retained
    source = cfg.output_path / "provenance" / "demultiplexing" / "plan.json"
    if not source.is_file():
        source = (
            cfg.static.paths.output_dir
            / cfg.run.run_id
            / "provenance"
            / "demultiplexing"
            / "plan.json"
        )
    if not source.is_file():
        return None
    result = json.loads(source.read_text())
    result = next((item for item in pass_plans(result) if project in item["projects"]), result)
    if project not in result["projects"]:
        raise DemultiplexingPlanError(
            f"Project {project} is absent from conversion provenance {source}"
        )
    destination.mkdir(parents=True, exist_ok=True)
    raw_rows = list(csv.reader(io.StringIO(result["sample_sheet_content"])))
    marker = "[BCLConvert_Data]" if result.get("sample_sheet_version", 2) == 2 else "[Data]"
    index = raw_rows.index([marker]) + 1
    header = raw_rows[index]
    project_column = header.index("Sample_Project")
    own_rows = [row for row in raw_rows[index + 1 :] if row and row[project_column] == project]
    sheet = destination / "SampleSheet.project.csv"
    with sheet.open("w", encoding="utf-8", newline="") as handle:
        csv.writer(handle, lineterminator="\n").writerows(raw_rows[: index + 1] + own_rows)
    retained = {
        key: value
        for key, value in result.items()
        if key not in {"sample_sheet_content", "sources", "rows", "projects", "command"}
    }
    retained.update(
        {
            "project": project,
            "sources": [item for item in result["sources"] if item["project"] == project],
            "rows": [item for item in result["rows"] if item["project"] == project],
            "project_sample_sheet_sha256": _sha256(sheet),
        }
    )
    (destination / "plan.json").write_text(json.dumps(retained, indent=2, sort_keys=True) + "\n")
    return destination
