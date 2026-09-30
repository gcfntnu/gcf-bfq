"""Sequencer-only QC, independent of Excel metadata and project analysis.

All generated inputs, summaries and reports belong to one demultiplexing execution
and live in ``Stats/sequencing_qc``. State and email delivery are caller concerns.
Counts below are clusters (read pairs for paired-end runs), not R1+R2 records.
"""

from __future__ import annotations

import csv
import json
import logging
import math
import shlex
import shutil
import subprocess
import xml.etree.ElementTree as ET

from collections import defaultdict
from html import escape
from pathlib import Path

import yaml

from bcl2fastq_pipeline.interop import prepare_index_metrics, run_interop_csv

log = logging.getLogger(__name__)
REPORT_DIRECTORY = Path("Stats") / "sequencing_qc"
REPORT_FILENAME = "sequencing_qc.html"


def _number(value):
    """InterOp fields may contain an uncertainty, or an unavailable value."""
    try:
        value = float(str(value).split("+/-")[0].strip())
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def _count(value):
    number = _number(value)
    if number is None or number < 0 or not number.is_integer():
        raise ValueError(f"Invalid read count: {value!r}")
    return int(number)


def _percent(numerator, denominator):
    return 100 * numerator / denominator if numerator is not None and denominator else None


def _display(value, suffix=""):
    if value is None:
        return "Unavailable"
    return (f"{value:,}" if isinstance(value, int) else f"{value:,.2f}") + suffix


def _sheet_samples(path, warnings):
    """Read Illumina [Data]/[BCLConvert_Data] and the supported plain 10x CSV."""
    samples = {}
    if path is None or not Path(path).is_file():
        warnings.append("SampleSheet unavailable; planned samples cannot be enumerated.")
        return samples
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        header = None
        for values in csv.reader(handle):
            if not values or not any(value.strip() for value in values):
                continue
            first = values[0].strip()
            if first.startswith("["):
                header = None
                continue
            names = [value.strip().lower().replace("_", "") for value in values]
            if "sampleid" in names or "sample" in names:
                header = names
                continue
            if header is None:
                continue
            row = dict(zip(header, (value.strip() for value in values), strict=False))
            sid = row.get("sampleid") or row.get("sample")
            if not sid:
                continue
            project = row.get("sampleproject") or row.get("project") or "Unspecified"
            key = (project, sid)
            if key not in samples:
                samples[key] = {
                    "project": project,
                    "sample_id": sid,
                    "sample_name": row.get("samplename") or sid,
                    "expected": True,
                    "assigned_reads": None,
                    "sheet_rows": [],
                }
            samples[key]["sheet_rows"].append(
                {
                    "lane": row.get("lane", ""),
                    "index": "+".join(
                        value for value in (row.get("index"), row.get("index2")) if value
                    ),
                }
            )
    if not samples:
        warnings.append(
            "SampleSheet contains no readable sample rows; planned samples unavailable."
        )
    return samples


def _assign(samples, sid, project, reads, *, name=None, lane=None, index=None):  # noqa: PLR0913
    candidates = [
        key
        for key, sample in samples.items()
        if sid in (sample["sample_id"], sample["sample_name"])
        and (not project or project in ("default", sample["project"]))
    ]
    if len(candidates) > 1 and lane is not None:
        candidates = [
            key
            for key in candidates
            if any(
                (not row["lane"] or row["lane"] == str(lane))
                and (
                    not index
                    or not row["index"]
                    or row["index"].upper() == index.replace("-", "+").upper()
                )
                for row in samples[key].get("sheet_rows", [])
            )
        ]
    if len(candidates) == 1:
        key = candidates[0]
    else:
        # An ambiguous identifier cannot safely be attributed to a planned sample.
        key = (project or "Unspecified", sid)
        samples.setdefault(
            key,
            {
                "project": key[0],
                "sample_id": sid,
                "sample_name": name or sid,
                "expected": False,
                "assigned_reads": None,
            },
        )
    current = samples[key]["assigned_reads"]
    samples[key]["assigned_reads"] = reads + (current if current is not None else 0)
    samples[key].setdefault("observed_lanes", set()).add(str(lane))


def _candidate_dirs(output):
    roots = [
        output / "Stats",
        output / "Reports",
        output / "Reports" / "legacy" / "Stats",
        output / "outs" / "fastq_path" / "Stats",
        output / "fastq_path" / "Stats",
        output / "outs" / "Stats",
    ]
    # mkfastq may add a flowcell directory underneath fastq_path (or the
    # configured --output-dir). Bound discovery to these documented layouts;
    # do not recursively scan FASTQ trees or prior report snapshots.
    for pattern in (
        "outs/fastq_path/*/Stats",
        "fastq_path/*/Stats",
        "*/outs/fastq_path/Stats",
        "*/outs/fastq_path/*/Stats",
        "*/Stats",
    ):
        roots.extend(sorted(output.glob(pattern)))
    return list(dict.fromkeys(roots))


def _find_source(output, tool):
    roots = _candidate_dirs(output)
    csv_path = next(
        (
            root / "Demultiplex_Stats.csv"
            for root in roots
            if (root / "Demultiplex_Stats.csv").is_file()
        ),
        None,
    )
    json_path = next(
        (root / "Stats.json" for root in roots if (root / "Stats.json").is_file()), None
    )
    xml_path = next(
        (
            root / "DemultiplexingStats.xml"
            for root in roots
            if (root / "DemultiplexingStats.xml").is_file()
        ),
        None,
    )
    name = str(tool[0] if isinstance(tool, (list, tuple)) and tool else tool or "").lower()
    module = "bclconvert" if "bcl-convert" in name or "bclconvert" in name else "bcl2fastq"
    if not name:
        module = "bclconvert" if csv_path else "bcl2fastq"
    source = csv_path if module == "bclconvert" else json_path or xml_path
    if source is None:
        # BFQ_TEST and legacy runs may retain a different actual demultiplexer output.
        source = csv_path or json_path or xml_path
        module = "bclconvert" if source == csv_path and source else "bcl2fastq"
    if source is None:
        raise RuntimeError(
            "No demultiplexer statistics found in Reports/Stats or supported mkfastq outputs"
        )
    return module, source


def _csv_demux(path, samples):
    lanes = {}
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not {"Lane", "SampleID", "# Reads"}.issubset(reader.fieldnames or []):
            raise ValueError(f"Missing Lane, SampleID or # Reads columns in {path}")
        for row in reader:
            sid, count = row["SampleID"], _count(row["# Reads"])
            lane = lanes.setdefault(row["Lane"], {"total_reads": 0, "undetermined_reads": None})
            lane["total_reads"] += count
            if sid.lower() == "undetermined":
                lane["undetermined_reads"] = (lane["undetermined_reads"] or 0) + count
            else:
                _assign(
                    samples,
                    sid,
                    row.get("Sample_Project"),
                    count,
                    lane=row["Lane"],
                    index=row.get("Index"),
                )
    unknown = None
    unknown_path = path.parent / "Top_Unknown_Barcodes.csv"
    if unknown_path.is_file():
        unknown = []
        with unknown_path.open(encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle):
                unknown.append(
                    {
                        "lane": row["Lane"],
                        "index": row["index"],
                        "index2": row.get("index2", ""),
                        "reads": _count(row["# Reads"]),
                    }
                )
    return lanes, unknown


def _json_demux(path, samples):
    data = json.loads(path.read_text())
    lanes = {}
    for result in data.get("ConversionResults", []):
        total, undetermined = (
            result.get("TotalClustersPF"),
            result.get("Undetermined", {}).get("NumberReads"),
        )
        determined = 0
        for sample in result.get("DemuxResults", []):
            count = _count(sample["NumberReads"])
            determined += count
            _assign(
                samples,
                sample["SampleId"],
                sample.get("SampleProject"),
                count,
                name=sample.get("SampleName"),
                lane=result["LaneNumber"],
                index=(sample.get("IndexMetrics") or [{}])[0].get("IndexSequence"),
            )
        if total is None and undetermined is not None:
            total = determined + _count(undetermined)
        lanes[str(result["LaneNumber"])] = {
            "total_reads": _count(total) if total is not None else None,
            "undetermined_reads": _count(undetermined) if undetermined is not None else None,
        }
    unknown = None
    if "UnknownBarcodes" in data:
        unknown = [
            {"lane": str(lane["Lane"]), "index": barcode, "index2": "", "reads": _count(count)}
            for lane in data["UnknownBarcodes"]
            for barcode, count in lane["Barcodes"].items()
        ]
    return lanes, unknown


def _xml_demux(path, samples):
    root = ET.parse(path).getroot()
    lanes = {}
    for project in root.findall(".//Project"):
        pid = project.get("name")
        for sample in project.findall("Sample"):
            sid = sample.get("name")
            # The all barcode includes the individual barcodes; never count both.
            barcodes = sample.findall("Barcode[@name='all']") or sample.findall("Barcode")
            for barcode in barcodes:
                for lane in barcode.findall("Lane"):
                    count = _count(lane.findtext("BarcodeCount"))
                    lane_id = lane.get("number")
                    current = lanes.setdefault(
                        lane_id, {"total_reads": None, "undetermined_reads": None}
                    )
                    if pid == "all" and sid == "all":
                        current["total_reads"] = count
                    elif pid == "default" and sid == "all":
                        current["undetermined_reads"] = count
                    elif sid != "all" and pid != "all" and sid.lower() != "undetermined":
                        _assign(samples, sid, pid, count, lane=lane_id)
    return lanes, None


def _run_metadata(cfg, inputs, warnings):
    output = Path(cfg.output_path)
    source_root = Path(cfg.run.flowcell_path) if cfg.run.flowcell_path else output
    runinfo = output / "RunInfo.xml"
    source = source_root / "RunInfo.xml"
    if not runinfo.is_file() and source.is_file():
        shutil.copy2(source, runinfo)
    if not runinfo.is_file():
        raise RuntimeError(
            "RunInfo.xml is unavailable in both output and instrument run directories"
        )
    shutil.copy2(runinfo, inputs / runinfo.name)
    run = ET.parse(runinfo).getroot().find("Run")
    if run is None:
        raise ValueError(f"No Run element in {runinfo}")
    geometry = []
    for read in run.findall("Reads/Read"):
        geometry.append(
            ("I" if read.get("IsIndexedRead") == "Y" else "R") + read.get("NumCycles", "?")
        )
    parameters = next(
        (
            root / filename
            for root in (output, source_root)
            for filename in ("RunParameters.xml", "runParameters.xml")
            if (root / filename).is_file()
        ),
        None,
    )
    if parameters:
        shutil.copy2(parameters, inputs / "RunParameters.xml")
    else:
        warnings.append("RunParameters.xml unavailable; instrument settings are not included.")
    return " / ".join(geometry) or "Unavailable"


def _interop(cfg, inputs, warnings):
    output = Path(cfg.output_path)
    prepare_index_metrics(output)
    for command, filename, required in (
        ("interop_summary", "interop_summary.csv", output / "InterOp"),
        (
            "interop_index-summary",
            "interop_index-summary.csv",
            output / "InterOp" / "IndexMetricsOut.bin",
        ),
    ):
        destination = inputs / filename
        destination.unlink(missing_ok=True)
        if not required.exists():
            warnings.append(f"{command} unavailable: {required.name} is missing.")
            continue
        try:
            run_interop_csv(command, output, destination, inputs)
        except (OSError, RuntimeError, subprocess.SubprocessError) as error:
            warnings.append(f"{command} unavailable: {error}")


def _interop_metrics(path):
    """Read non-index, whole-lane rows; preserve the InterOp units (M reads/Gbp)."""
    lanes, summary, header, read_name, read_index = {}, {}, None, None, 0
    if not path.is_file():
        return lanes, summary
    for row in csv.reader(path.read_text().splitlines()):
        if not row:
            continue
        first = row[0].strip()
        if first == "Level":
            header, read_name = row, "summary"
        elif first.startswith("Read ") and len(row) == 1:
            header = None
            if "(I)" in first:
                read_name = None
            else:
                read_index += 1
                read_name = f"R{read_index}"
        elif first == "Lane":
            header = row
        elif first.startswith("Extracted"):
            break
        elif header and read_name == "summary":
            values = dict(zip(header, row, strict=False))
            if first == "Non-indexed":
                summary = {
                    "yield_gb": _number(values.get("Yield")),
                    "q30_percent": _number(values.get("%>=Q30")),
                }
            elif first == "Total":
                # Total can include indexes, so label this separately.
                summary["total_yield_gb"] = _number(values.get("Yield"))
                read_name, header = None, None
        elif header and read_name and first.isdigit():
            values = dict(zip(header, row, strict=False))
            if values.get("Surface", "-").strip() != "-":
                continue
            lane = lanes.setdefault(first, {})
            lane[f"{read_name}_q30_percent"] = _number(values.get("%>=Q30"))
            if read_index == 1:
                for source, target in (
                    ("Density", "density_k_mm2"),
                    ("Reads", "reads_m"),
                    ("Cluster PF", "pf_percent"),
                    ("Aligned", "phix_percent"),
                ):
                    lane[target] = _number(values.get(source))
    return lanes, summary


def _table(headers, rows):
    return (
        "<table><thead><tr>"
        + "".join(f"<th>{escape(header)}</th>" for header in headers)
        + "</tr></thead><tbody>"
        + "".join(
            "<tr>" + "".join(f"<td>{escape(str(value))}</td>" for value in row) + "</tr>"
            for row in rows
        )
        + "</tbody></table>"
    )


def _summaries(result):
    samples = result["samples"]
    lines = [
        f"Read geometry: {result['read_geometry']}",
        f"Demultiplexer statistics: {result['demultiplexer']}",
    ]
    planned = defaultdict(int)
    for sample in samples:
        if sample["expected"]:
            planned[sample["project"]] += 1
    lines.extend(
        f"{project}: {count} unique planned samples in SampleSheet."
        for project, count in sorted(planned.items())
    )
    lines.extend(
        [
            f"Total PF clusters (read pairs for paired-end runs): {_display(result['total_reads'])}",
            f"Undetermined: {_display(result['undetermined_reads'])} ({_display(result['undetermined_percent'], '%')})",
            f"Non-index yield: {_display(result.get('yield_gb'), ' Gbp')}",
            f"Expected samples with explicit zero assigned reads: {result['zero_read_sample_count']}",
            f"Expected samples with unavailable assignment counts: {result['unavailable_sample_count']}",
        ]
    )
    columns = [
        ("Lane", "lane"),
        ("Density (K/mm²)", "density_k_mm2"),
        ("Reads (M)", "reads_m"),
        ("PF (%)", "pf_percent"),
        ("PhiX (%)", "phix_percent"),
        ("R1 Q30 (%)", "R1_q30_percent"),
        ("R2 Q30 (%)", "R2_q30_percent"),
        ("Undetermined (%)", "undetermined_percent"),
    ]
    rows = [
        [lane["lane"], *[_display(lane.get(key)) for _, key in columns[1:]]]
        for lane in result["lane_metrics"]
    ]
    text = "\n".join(lines)
    html = "<p>" + escape(text).replace("\n", "<br>\n") + "</p><h3>Flowcell metrics</h3>"
    if rows:
        html += _table([title for title, _ in columns], rows)
        text += "\nFlowcell metrics:\n" + "\n".join(
            "; ".join(f"{title}: {value}" for (title, _), value in zip(columns, row, strict=True))
            for row in rows
        )
    else:
        html += "<p>Lane metrics unavailable.</p>"
        text += "\nLane metrics unavailable."
    if result["unknown_indexes"] is None:
        unknown_text = "Unknown-index information unavailable."
    elif not result["unknown_indexes"]:
        unknown_text = "Unknown-index statistics contain no barcodes."
    else:
        unknown_text = "Top unknown indexes (up to 5):\n" + "\n".join(
            f"Lane {item['lane']}: {item['index']}{'+' + item['index2'] if item['index2'] else ''}: {item['reads']:,} clusters"
            for item in sorted(
                result["unknown_indexes"], key=lambda item: item["reads"], reverse=True
            )[:5]
        )
    text += "\n" + unknown_text
    html += "<p>" + escape(unknown_text).replace("\n", "<br>\n") + "</p>"
    if result["warnings"]:
        warning_text = "\n".join(result["warnings"])
        text += "\nUnavailable inputs / report notes:\n" + warning_text
        html += (
            "<h3>Unavailable inputs / report notes</h3><p>"
            + escape(warning_text).replace("\n", "<br>\n")
            + "</p>"
        )
    return text, html


def _configuration(result):
    # An owned custom section keeps expected zero/unavailable samples visible even
    # when a MultiQC demultiplexer module suppresses them. Never inherit workflow
    # custom data, sample filters or external upload credentials.
    data = {
        f"{sample['project']} / {sample['sample_id']}": {
            "project": sample["project"],
            "sample_id": sample["sample_id"],
            "expected": "Yes" if sample["expected"] else "No (statistics only)",
            "assigned_reads": sample["assigned_reads"]
            if sample["assigned_reads"] is not None
            else "Unavailable",
        }
        for sample in result["samples"]
    }
    return {
        "title": f"Sequencing QC — {result['run_id']}",
        "intro_text": "Genomics Core Facility, NTNU. Sequencing and demultiplexing QC before analysis. Counts are clusters (read pairs for paired-end runs). No automatic QC gate is applied.",
        "megaqc_upload": False,
        "megaqc_url": None,
        "megaqc_access_token": None,
        "custom_data": {
            "bfq_sequencing_summary": {
                "section_name": "BFQ sequencing summary",
                "plot_type": "html",
                "data": result["summary_html"],
            },
            "bfq_sample_assignment": {
                "section_name": "Expected samples and demultiplexing assignment",
                "plot_type": "table",
                "description": "Missing statistics are unavailable, not zero. Explicit zero counts are retained. SampleSheet identifiers require no submission workbook.",
                "pconfig": {
                    "id": "bfq_sample_assignment",
                    "title": "Sample assignment",
                    "col1_header": "Project / sample",
                },
                "headers": {
                    "project": {"title": "Project"},
                    "sample_id": {"title": "Sample ID"},
                    "expected": {"title": "In SampleSheet"},
                    "assigned_reads": {"title": "Assigned clusters", "format": "{:,.0f}"},
                },
                "data": data,
            },
        },
        "module_order": ["custom_content", "interop", result["demultiplexer"]],
    }


def generate(cfg, tool=None):
    """Generate one report, returning durable JSON-safe notification context.

    Essential malformed demultiplexer/run metadata or MultiQC failures raise with
    diagnostic paths. Optional missing metrics remain explicitly unavailable.
    """
    output = Path(cfg.output_path).resolve()
    report_dir = output / REPORT_DIRECTORY
    inputs = report_dir / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    warnings = []
    samples = _sheet_samples(cfg.run.sample_sheet, warnings)
    module, source = _find_source(output, tool)
    geometry = _run_metadata(cfg, inputs, warnings)
    # Avoid duplicate native/legacy inputs and ensure the BCL Convert module sees
    # RunInfo.xml in exactly the same directory as Demultiplex_Stats.csv.
    for name in (
        "Demultiplex_Stats.csv",
        "Quality_Metrics.csv",
        "Top_Unknown_Barcodes.csv",
        "Stats.json",
        "DemultiplexingStats.xml",
    ):
        target = inputs / name
        target.unlink(missing_ok=True)
        candidate = source.parent / name
        if candidate.is_file() and (
            module == "bclconvert"
            and name.endswith(".csv")
            or module == "bcl2fastq"
            and name.endswith((".json", ".xml"))
        ):
            shutil.copy2(candidate, target)
    try:
        parser = (
            _csv_demux
            if source.suffix == ".csv"
            else _json_demux
            if source.suffix == ".json"
            else _xml_demux
        )
        lanes, unknown = parser(source, samples)
        if not lanes:
            raise ValueError("No lane assignment statistics found")
    except (OSError, ValueError, KeyError, TypeError, ET.ParseError) as error:
        raise RuntimeError(f"Unable to read demultiplexer statistics {source}: {error}") from error
    for sample in samples.values():
        observed = sample.pop("observed_lanes", set())
        sheet_rows = sample.pop("sheet_rows", [])
        expected_lanes = {row["lane"] for row in sheet_rows if row["lane"]}
        if any(not row["lane"] for row in sheet_rows):
            expected_lanes.update(lanes)
        missing_lanes = expected_lanes - observed
        if sample["expected"] and sample["assigned_reads"] is not None and missing_lanes:
            sample["observed_assigned_reads"] = sample["assigned_reads"]
            sample["assigned_reads"] = None
            warnings.append(
                f"{sample['project']} / {sample['sample_id']}: assignment total unavailable; "
                f"statistics missing for expected lanes {', '.join(sorted(missing_lanes))}."
            )
    unmatched = [sample["sample_id"] for sample in samples.values() if not sample["expected"]]
    if unmatched:
        warnings.append(
            "Statistics identifiers not uniquely matched to the SampleSheet: "
            + ", ".join(sorted(set(unmatched)))
            + ". Counts are retained as statistics-only samples."
        )
    if source.suffix == ".xml":
        warnings.append(
            "Stats.json unavailable: XML assignment counts are shown in the BFQ section; native bcl2fastq MultiQC charts and unknown indexes may be unavailable."
        )
    _interop(cfg, inputs, warnings)
    interop_lanes, summary = _interop_metrics(inputs / "interop_summary.csv")
    for lane_id in set(lanes) | set(interop_lanes):
        lane = lanes.setdefault(lane_id, {"total_reads": None, "undetermined_reads": None})
        lane.update(interop_lanes.get(lane_id, {}))
        lane["lane"] = lane_id
        lane["undetermined_percent"] = _percent(lane["undetermined_reads"], lane["total_reads"])
    total = (
        sum(lane["total_reads"] for lane in lanes.values())
        if all(lane["total_reads"] is not None for lane in lanes.values())
        else None
    )
    undetermined = (
        sum(lane["undetermined_reads"] for lane in lanes.values())
        if all(lane["undetermined_reads"] is not None for lane in lanes.values())
        else None
    )
    expected = [sample for sample in samples.values() if sample["expected"]]
    result = {
        "run_id": cfg.run.run_id,
        "projects": sorted({sample["project"] for sample in samples.values()}),
        "report_path": str(report_dir / REPORT_FILENAME),
        "read_geometry": geometry,
        "demultiplexer": module,
        "stats_source": str(source),
        "samples": sorted(
            samples.values(), key=lambda sample: (sample["project"], sample["sample_id"])
        ),
        "lane_metrics": sorted(lanes.values(), key=lambda lane: lane["lane"]),
        "total_reads": total,
        "undetermined_reads": undetermined,
        "undetermined_percent": _percent(undetermined, total),
        "unknown_indexes": unknown,
        "warnings": warnings,
        "zero_read_sample_count": sum(sample["assigned_reads"] == 0 for sample in expected),
        "unavailable_sample_count": sum(sample["assigned_reads"] is None for sample in expected),
        **summary,
    }
    result["summary_text"], result["summary_html"] = _summaries(result)
    config_path = report_dir / "multiqc_config.yaml"
    config_path.write_text(
        yaml.safe_dump(_configuration(result), sort_keys=False, allow_unicode=True)
    )
    report_path = Path(result["report_path"])
    report_path.unlink(missing_ok=True)
    log_path = report_dir / "multiqc.log"
    command = [
        "multiqc",
        *shlex.split(cfg.static.commands.get("multiqc_options", "")),
        "--force",
        "--no-megaqc-upload",
        "--config",
        str(config_path),
        "--outdir",
        str(report_dir),
        "--filename",
        REPORT_FILENAME,
        "-m",
        "custom_content",
        "-m",
        "interop",
        "-m",
        module,
        str(inputs),
    ]
    try:
        with log_path.open("w") as handle:
            subprocess.check_call(command, cwd=report_dir, stdout=handle, stderr=subprocess.STDOUT)
        if not report_path.is_file() or not report_path.stat().st_size:
            raise RuntimeError("MultiQC returned without a nonempty HTML report")
    except (OSError, RuntimeError, subprocess.SubprocessError) as error:
        tail = "\n".join(log_path.read_text(errors="replace").splitlines()[-30:])
        raise RuntimeError(
            f"Sequencing QC generation failed: {error}. Log: {log_path}\n{tail}"
        ) from error
    (report_dir / "summary.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    return result
