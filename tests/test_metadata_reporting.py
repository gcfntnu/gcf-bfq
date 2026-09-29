"""Real workbook validation in completion/error reporting; no SMTP delivery."""

import json

from types import SimpleNamespace

from bcl2fastq_pipeline.preflight import PreflightValidationError
from configmaker.validation import validate_inputs
from openpyxl import Workbook

from bcl2fastq_pipeline import misc


def metadata_pair(tmp_path):
    sheet = tmp_path / "SampleSheet.csv"
    sheet.write_text(
        "[Data]\nLane,Sample_ID,Sample_Project,index\n"
        "1,001,GCF-2026-043,ACGT\n2,001,GCF-2026-043,ACGT\n"
        "1,002,GCF-2026-043,TGCA\n"
    )
    form = tmp_path / "Sample-Submission-Form.xlsx"
    workbook = Workbook()
    customer = workbook.active
    customer.title = "Sample-Submission-Form"
    for _ in range(14):
        customer.append([])
    customer.append(["Unique Sample ID", "Sample Group"])
    customer.append(["001", "<control>"])
    customer.append(["002", None])
    customer.append(["extra", "treated"])
    lab = workbook.create_sheet("INFO (GCF-lab only)")
    lab.append(["Sample_ID"])
    workbook.save(form)
    return SimpleNamespace(
        output_path=tmp_path,
        run=SimpleNamespace(sample_sheet=sheet, sample_submission_form=form),
    )


def test_planned_email_counts_unique_lane_samples_and_preserves_groups(tmp_path):
    cfg = metadata_pair(tmp_path)
    html = misc.parseSampleSheetMetrics(cfg, projects=["GCF-2026-043"])
    assert "Planned input samples" in html
    assert "2 unique samples in SampleSheet" in html
    assert "across 3 SampleSheet rows" in html
    assert "3 samples in effective submission metadata" in html
    assert "1 additional submission samples allowed" in html
    assert "Sample_Group has 2 unique values" in html
    assert "Missing Sample_Group for 1 samples" in html
    assert "&lt;control&gt;" in html
    assert "<control>" not in html


def test_malformed_workbook_is_reported_and_corrected_bytes_are_revalidated(tmp_path):
    cfg = metadata_pair(tmp_path)
    original = cfg.run.sample_submission_form.read_bytes()
    cfg.run.sample_submission_form.write_bytes(b"not an XLSX workbook")
    html = misc.parseSampleSheetMetrics(cfg)
    assert "Current input validation failed" in html
    assert "Sample-Submission-Form.xlsx" in html
    cfg.run.sample_submission_form.write_bytes(original)
    assert "Current input validation failed" not in misc.parseSampleSheetMetrics(cfg)


def test_error_identity_includes_actual_input_bytes_not_check_time(tmp_path):
    cfg = metadata_pair(tmp_path)
    cfg.run.sample_submission_form.write_bytes(b"invalid workbook one")

    def identity():
        validation = validate_inputs([cfg.run.sample_sheet], [cfg.run.sample_submission_form])
        error = PreflightValidationError(validation)
        return misc.error_failure_signature("demultiplexing", (type(error), error, None), "failed")

    initial = identity()
    assert identity() == initial
    cfg.run.sample_submission_form.write_bytes(b"invalid workbook two")
    assert identity() != initial


def test_error_report_keeps_full_json_after_abbreviated_exception(tmp_path, monkeypatch):
    cfg = metadata_pair(tmp_path)
    cfg.run.run_id = "260918_MN00686_0026_A000HCMFHF"
    cfg.static = SimpleNamespace(paths=SimpleNamespace(report_dir=tmp_path / "reports"))
    cfg.run.sample_submission_form.write_bytes(b"malformed")
    validation = validate_inputs([cfg.run.sample_sheet], [cfg.run.sample_submission_form])
    error = RuntimeError("Input preflight failed (summary)")
    error.validation_result = validation
    monkeypatch.setattr(misc.PipelineConfig, "get", lambda: cfg)
    path = misc.write_error_report((type(error), error, None), "preflight error")
    report = json.loads(path.read_text().split("Complete structured input preflight report:\n")[1])
    assert report["status"] == "failed"
    assert report["inputs"] == validation.to_dict()["inputs"]
    assert report["errors"]


def test_analysis_discovery_is_separate_from_planned_metadata(tmp_path):
    cfg = metadata_pair(tmp_path)
    (tmp_path / "configmaker-analysis-GCF-2026-043.json").write_text(
        json.dumps(
            {
                "schema_version": 1,
                "kind": "fastq_discovery",
                "sample_count": 1,
                "missing_sample_ids": ["002"],
                "samples": [{"sample_id": "001", "flowcells": ["TEST"]}],
            }
        )
    )
    observed = misc.analysisSampleMetrics(cfg, ["GCF-2026-043"])
    assert "1 samples discovered in FASTQs" in observed
    assert "Planned samples without FASTQs: 002" in observed
    assert "2 unique samples in SampleSheet" in misc.parseSampleSheetMetrics(cfg)


def test_legacy_analysis_without_discovery_report_remains_reportable(tmp_path):
    cfg = metadata_pair(tmp_path)
    assert "discovery summary unavailable" in misc.analysisSampleMetrics(cfg, ["GCF-2026-043"])
