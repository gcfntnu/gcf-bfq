"""Small, workflow-independent analysis QC snapshots for durable notifications.

Sample discovery is a normal configmaker output. fastp is optional, so summaries
always disclose report coverage and never treat an absent metric as zero. Capture
these values before finalization; email retries must not require the workdir.
"""

from __future__ import annotations

import json
import math

from html import escape

from bcl2fastq_pipeline import analysis_snapshots


def _count(value):
    if type(value) is not int or value < 0:
        raise ValueError("Expected a non-negative integer count")
    return value


def _discovery(output, project):
    try:
        data = json.loads((output / f"configmaker-analysis-{project}.json").read_text())
        if data.get("schema_version") != 1 or data.get("kind") != "fastq_discovery":
            raise ValueError("Unsupported sample discovery summary")
        count = _count(data["sample_count"])
        missing = data.get("missing_sample_ids", [])
        if not isinstance(missing, list) or not all(isinstance(value, str) for value in missing):
            raise ValueError("Invalid missing sample list")
        return count, missing
    except (OSError, ValueError, KeyError, TypeError, AttributeError):
        return None, None


def _fastp_counts(path):
    data = json.loads(path.read_text())
    if not isinstance(data, dict):
        raise ValueError("Expected a fastp summary object")
    summary = data["summary"]
    before = _count(summary["before_filtering"]["total_reads"])
    after = summary["after_filtering"]
    retained = _count(after["total_reads"])
    if retained > before:
        raise ValueError("Output reads exceed input reads")
    bases = q30_bases = None
    try:
        bases = _count(after["total_bases"])
        if "q30_bases" in after:
            q30_bases = _count(after["q30_bases"])
            if q30_bases > bases:
                raise ValueError("Q30 bases exceed total bases")
        else:
            rate = after["q30_rate"]
            if type(rate) not in (int, float) or not math.isfinite(rate) or not 0 <= rate <= 1:
                raise ValueError("Invalid Q30 rate")
            q30_bases = rate * bases
    except (KeyError, ValueError, TypeError, OverflowError):
        bases = q30_bases = None
    return before, retained, bases, q30_bases


def _fastp(cfg, project):
    result = {
        "reports_found": 0,
        "reports_used": 0,
        "input_reads": None,
        "output_reads": None,
        "retention_pct": None,
        "q30_pct": None,
        "q30_reports_used": 0,
        "warnings": [],
    }
    pipeline = getattr(cfg.run, "pipeline", None)
    if not isinstance(pipeline, str) or not pipeline or "/" in pipeline or pipeline in {".", ".."}:
        result["warnings"].append("fastp summary unavailable: workflow is not recorded.")
        return result
    try:
        workdir = analysis_snapshots.workdir_path(cfg.run.run_id, project)
        marker = json.loads((workdir / analysis_snapshots.MARKER).read_text())
        if marker.get("run_id") != cfg.run.run_id or marker.get("project") != project:
            raise ValueError("Workdir belongs to another run")
    except (OSError, ValueError, KeyError, TypeError, AttributeError, RuntimeError):
        result["warnings"].append(
            "fastp summary unavailable: original analysis workdir cannot be verified."
        )
        return result
    logs = workdir / "data" / "tmp" / pipeline / "bfq" / "logs"
    # Workflow logs are commonly symlinks. Deduplicate aliases of one report.
    try:
        paths = sorted({path.resolve() for path in logs.rglob("*.fastp.json")})
    except (OSError, RuntimeError):
        result["warnings"].append("fastp summary unavailable: report directory cannot be read.")
        return result
    result["reports_found"] = len(paths)
    total_input = total_output = total_bases = q30_bases = 0
    for path in paths:
        try:
            before, after, bases, q30 = _fastp_counts(path)
        except (OSError, ValueError, KeyError, TypeError, AttributeError, OverflowError):
            continue
        result["reports_used"] += 1
        total_input += before
        total_output += after
        if bases is not None:
            result["q30_reports_used"] += 1
            total_bases += bases
            q30_bases += q30
    if result["reports_used"]:
        result["input_reads"] = total_input
        result["output_reads"] = total_output
        if total_input:
            result["retention_pct"] = 100 * total_output / total_input
    if total_bases:
        result["q30_pct"] = 100 * q30_bases / total_bases
    if not paths:
        result["warnings"].append("fastp summary unavailable: no fastp reports for this workflow.")
    elif result["reports_used"] != len(paths):
        result["warnings"].append(
            "Some fastp reports could not be read; totals cover valid reports only."
        )
    if result["q30_reports_used"] < result["reports_used"]:
        result["warnings"].append("Q30 covers only reports with valid after-filtering base counts.")
    return result


def _format_count(value):
    return "Unavailable" if value is None else f"{value:,}"


def _format_pct(value):
    return "Unavailable" if value is None else f"{value:.2f}%"


def _render(projects):
    headings = (
        "Project",
        "FASTQ samples",
        "Planned samples without FASTQs",
        "fastp reports used/found",
        "Input reads",
        "Retained reads",
        "Read retention",
        "After-filter Q30 bases",
    )
    rows = []
    notes = [
        "FASTQ sample counts describe analysis initialization. fastp metrics are optional; "
        "read counts include both mates, not read pairs. Report coverage is shown explicitly "
        "and is not a count of samples passing QC. Q30 is weighted by retained bases."
    ]
    for project in projects:
        fastp = project["fastp"]
        missing = project["missing_sample_ids"]
        rows.append(
            (
                project["project"],
                _format_count(project["sample_count"]),
                _format_count(len(missing) if missing is not None else None),
                f"{fastp['reports_used']}/{fastp['reports_found']}",
                _format_count(fastp["input_reads"]),
                _format_count(fastp["output_reads"]),
                _format_pct(fastp["retention_pct"]),
                _format_pct(fastp["q30_pct"]) + f" ({fastp['q30_reports_used']} reports)",
            )
        )
        if missing:
            notes.append(
                f"{project['project']} planned samples without FASTQs: {', '.join(missing)}."
            )
        if project["sample_count"] is None:
            notes.append(f"{project['project']}: FASTQ discovery summary unavailable.")
        notes.extend(f"{project['project']}: {warning}" for warning in fastp["warnings"])
    html = "<strong>Analysis QC overview</strong><table><thead><tr>"
    html += "".join(f"<th>{escape(value)}</th>" for value in headings) + "</tr></thead><tbody>"
    html += "".join(
        "<tr>" + "".join(f"<td>{escape(value)}</td>" for value in row) + "</tr>" for row in rows
    )
    html += "</tbody></table>" + "".join(f"<p>{escape(note)}</p>" for note in notes)
    lines = ["Analysis QC overview"]
    for row in rows:
        lines.append("; ".join(f"{key}: {value}" for key, value in zip(headings, row, strict=True)))
    return "\n".join([*lines, *notes]), html


def collect(cfg, projects):
    """Return a JSON-compatible snapshot; missing optional outputs are explicit."""
    summaries = []
    for project in sorted(set(projects)):
        if "/" in project or project in {".", ".."}:
            raise ValueError(f"Invalid analysis project: {project!r}")
        count, missing = _discovery(cfg.output_path, project)
        summaries.append(
            {
                "project": project,
                "sample_count": count,
                "missing_sample_ids": missing,
                "fastp": _fastp(cfg, project),
            }
        )
    text, html = _render(summaries)
    return {"schema_version": 1, "projects": summaries, "summary_text": text, "summary_html": html}
