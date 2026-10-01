"""Optional read-filter QC is captured before workdir removal; no email delivery."""

import json
import shutil

from types import SimpleNamespace

import pytest

from bcl2fastq_pipeline import analysis_qc, analysis_snapshots

PROJECT = "GCF-2026-043"


@pytest.fixture
def analysis_cfg(tmp_path, monkeypatch):
    monkeypatch.setenv("TMPDIR", str(tmp_path / "work"))
    output = tmp_path / "output"
    output.mkdir()
    cfg = SimpleNamespace(
        output_path=output,
        run=SimpleNamespace(run_id="260930_MN00686_0026_TEST", pipeline="metagenome"),
    )
    (output / f"configmaker-analysis-{PROJECT}.json").write_text(
        json.dumps(
            {
                "schema_version": 1,
                "kind": "fastq_discovery",
                "sample_count": 3,
                "missing_sample_ids": ["missing<&>"],
            }
        )
    )
    workdir = analysis_snapshots.workdir_path(cfg.run.run_id, PROJECT)
    workdir.mkdir(parents=True)
    analysis_snapshots.identify_workdir(workdir, cfg.run.run_id, PROJECT)
    return cfg


def report(cfg, name, before=100, after=80, bases=800, q30=720):  # noqa: PLR0913
    logs = analysis_snapshots.workdir_path(cfg.run.run_id, PROJECT) / "data/tmp/metagenome/bfq/logs"
    path = logs / name / f"{name}.fastp.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(
            {
                "summary": {
                    "before_filtering": {"total_reads": before},
                    "after_filtering": {
                        "total_reads": after,
                        "total_bases": bases,
                        "q30_bases": q30,
                    },
                }
            }
        )
    )
    return path


def test_fastp_aggregation_weights_retention_by_reads_and_q30_by_bases(analysis_cfg):
    report(analysis_cfg, "sample1", before=100, after=80, bases=800, q30=720)
    report(analysis_cfg, "sample2", before=300, after=120, bases=2400, q30=1200)
    snapshot = analysis_qc.collect(analysis_cfg, [PROJECT])
    project = snapshot["projects"][0]
    assert project["sample_count"] == 3
    assert project["missing_sample_ids"] == ["missing<&>"]
    fastp = project["fastp"]
    assert fastp["input_reads"] == 400
    assert fastp["output_reads"] == 200
    assert fastp["retention_pct"] == 50
    assert fastp["q30_pct"] == 60
    assert fastp["reports_used"] == fastp["reports_found"] == 2
    assert "missing&lt;&amp;&gt;" in snapshot["summary_html"]
    assert "missing<&>" not in snapshot["summary_html"]
    assert "read counts include both mates, not read pairs" in snapshot["summary_text"]
    assert json.loads(json.dumps(snapshot)) == snapshot


def test_snapshot_remains_available_after_workdir_removal(analysis_cfg):
    report(analysis_cfg, "sample1")
    snapshot = analysis_qc.collect(analysis_cfg, [PROJECT])
    shutil.rmtree(analysis_snapshots.workdir_path(analysis_cfg.run.run_id, PROJECT))
    assert snapshot["projects"][0]["fastp"]["retention_pct"] == 80
    assert "80.00%" in snapshot["summary_html"]


def test_no_fastp_workflow_explicitly_reports_unavailable(analysis_cfg):
    snapshot = analysis_qc.collect(analysis_cfg, [PROJECT])
    fastp = snapshot["projects"][0]["fastp"]
    assert fastp["input_reads"] is None
    assert fastp["retention_pct"] is None
    assert "no fastp reports for this workflow" in snapshot["summary_text"]
    assert "Unavailable" in snapshot["summary_html"]


def test_missing_legacy_workdir_preserves_discovery_counts(analysis_cfg):
    shutil.rmtree(analysis_snapshots.workdir_path(analysis_cfg.run.run_id, PROJECT))
    snapshot = analysis_qc.collect(analysis_cfg, [PROJECT])
    assert snapshot["projects"][0]["sample_count"] == 3
    assert "workdir cannot be verified" in snapshot["summary_text"]


def test_workdir_for_other_flowcell_is_not_read(analysis_cfg):
    report(analysis_cfg, "sample1")
    workdir = analysis_snapshots.workdir_path(analysis_cfg.run.run_id, PROJECT)
    analysis_snapshots.identify_workdir(workdir, "260930_OTHER_RUN", PROJECT)
    snapshot = analysis_qc.collect(analysis_cfg, [PROJECT])
    assert snapshot["projects"][0]["fastp"]["input_reads"] is None


def test_deduplicates_report_symlinks(analysis_cfg):
    path = report(analysis_cfg, "sample1")
    (path.parent / "alias.fastp.json").symlink_to(path)
    assert (
        analysis_qc.collect(analysis_cfg, [PROJECT])["projects"][0]["fastp"]["reports_found"] == 1
    )


def test_invalid_reports_and_missing_q30_keep_partial_coverage_explicit(analysis_cfg):
    report(analysis_cfg, "good")
    bad = report(analysis_cfg, "bad")
    bad.write_text("malformed JSON")
    partial = report(analysis_cfg, "without_q30")
    data = json.loads(partial.read_text())
    data["summary"]["after_filtering"].pop("q30_bases")
    partial.write_text(json.dumps(data))
    snapshot = analysis_qc.collect(analysis_cfg, [PROJECT])
    fastp = snapshot["projects"][0]["fastp"]
    assert fastp["reports_found"] == 3
    assert fastp["reports_used"] == 2
    assert fastp["q30_reports_used"] == 1
    assert fastp["input_reads"] == 200
    assert fastp["q30_pct"] == 90
    assert "totals cover valid reports only" in snapshot["summary_text"]
    assert "Q30 covers only reports" in snapshot["summary_text"]


def test_q30_rate_fallback_and_zero_reads(analysis_cfg):
    path = report(analysis_cfg, "rate", before=0, after=0, bases=0, q30=0)
    data = json.loads(path.read_text())
    data["summary"]["after_filtering"].pop("q30_bases")
    data["summary"]["after_filtering"]["q30_rate"] = 0.95
    path.write_text(json.dumps(data))
    fastp = analysis_qc.collect(analysis_cfg, [PROJECT])["projects"][0]["fastp"]
    assert fastp["input_reads"] == 0
    assert fastp["output_reads"] == 0
    assert fastp["retention_pct"] is None
    assert fastp["q30_pct"] is None
    report(analysis_cfg, "nonzero", before=100, after=80, bases=800, q30=760)
    fastp = analysis_qc.collect(analysis_cfg, [PROJECT])["projects"][0]["fastp"]
    assert fastp["q30_pct"] == 95


@pytest.mark.parametrize("contents", ["[]", "null", "{}", "broken"])
def test_missing_or_malformed_discovery_is_not_reported_as_zero(analysis_cfg, contents):
    (analysis_cfg.output_path / f"configmaker-analysis-{PROJECT}.json").write_text(contents)
    snapshot = analysis_qc.collect(analysis_cfg, [PROJECT])
    assert snapshot["projects"][0]["sample_count"] is None
    assert snapshot["projects"][0]["missing_sample_ids"] is None
    assert "discovery summary unavailable" in snapshot["summary_text"]


@pytest.mark.parametrize("contents", ["[]", "null", "{}", '{"summary":null}', '{"summary":[]}'])
def test_malformed_optional_fastp_object_does_not_block_reporting(analysis_cfg, contents):
    path = report(analysis_cfg, "broken")
    path.write_text(contents)
    snapshot = analysis_qc.collect(analysis_cfg, [PROJECT])
    assert snapshot["projects"][0]["sample_count"] == 3
    fastp = snapshot["projects"][0]["fastp"]
    assert fastp["reports_found"] == 1
    assert fastp["reports_used"] == 0
    assert fastp["input_reads"] is None
    assert "Some fastp reports could not be read" in snapshot["summary_text"]


@pytest.mark.parametrize(
    "error", [KeyError("TMPDIR"), OSError("missing workdir"), RuntimeError("workdir unavailable")]
)
def test_workdir_resolution_failure_does_not_block_reporting(analysis_cfg, monkeypatch, error):
    def missing_workdir(*args):
        raise error

    monkeypatch.setattr(analysis_snapshots, "workdir_path", missing_workdir)
    snapshot = analysis_qc.collect(analysis_cfg, [PROJECT])
    assert snapshot["projects"][0]["sample_count"] == 3
    assert snapshot["projects"][0]["fastp"]["input_reads"] is None
    assert "workdir cannot be verified" in snapshot["summary_text"]
