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


def report_path(cfg, name, tool):
    logs = (
        analysis_snapshots.workdir_path(cfg.run.run_id, PROJECT)
        / "data/tmp"
        / cfg.run.pipeline
        / "bfq/logs"
    )
    path = logs / name / f"{name}.{tool}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def report(cfg, name, before=100, after=80, bases=800, q30=720):  # noqa: PLR0913
    path = report_path(cfg, name, "fastp")
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


def cutadapt_report(cfg, name, *, before=100, after=40, bases=800, paired=False):  # noqa: PLR0913
    path = report_path(cfg, name, "cutadapt")
    path.write_text(
        json.dumps(
            {
                "tag": "Cutadapt report, merged from the passes in gcf_passes",
                "schema_version": [0, 3],
                "cutadapt_version": "5.2",
                "input": {"path1": f"{name}.fastq", "path2": None, "paired": paired},
                "read_counts": {"input": before, "output": after},
                "basepair_counts": {"output": bases},
                "gcf_passes": [
                    {"pass": "adapter", "file": f"{name}.adapter.json"},
                    {"pass": "filter", "file": f"{name}.filter.json"},
                ],
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


def test_no_filtering_reports_explicitly_reports_unavailable(analysis_cfg):
    snapshot = analysis_qc.collect(analysis_cfg, [PROJECT])
    fastp = snapshot["projects"][0]["fastp"]
    assert fastp["input_reads"] is None
    assert fastp["retention_pct"] is None
    assert "no filtering reports for this workflow" in snapshot["summary_text"]
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
    assert "totals cover usable reports only" in snapshot["summary_text"]
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
    assert "Some filtering reports were invalid or ambiguous" in snapshot["summary_text"]


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


def test_cutadapt_only_smallrna_keeps_native_counts_and_unavailable_q30(analysis_cfg):
    analysis_cfg.run.pipeline = "smallrna"
    cutadapt_report(analysis_cfg, "sample1")
    snapshot = analysis_qc.collect(analysis_cfg, [PROJECT])
    filtering = snapshot["projects"][0]["filtering"]
    assert filtering["input_reads"] == 100
    assert filtering["output_reads"] == 40
    assert filtering["retention_pct"] == 40
    assert filtering["output_bases"] == 800
    assert filtering["bases_reports_used"] == 1
    assert filtering["q30_pct"] is None
    assert filtering["q30_reports_used"] == 0
    assert filtering["reports_found"] == filtering["reports_used"] == 1
    assert snapshot["projects"][0]["fastp"]["input_reads"] is None
    assert snapshot["projects"][0]["fastp"]["reports_found"] == 0
    for body in (snapshot["summary_text"], snapshot["summary_html"]):
        assert "Filtering reports used/found" in body
        assert "cutadapt: 1/1" in body
        assert "Unavailable (0 reports)" in body
        assert "cutadapt does not report Q30" in body
        assert "fastp" not in body
    assert json.loads(json.dumps(snapshot)) == snapshot


def test_mixed_tools_combine_reads_but_only_fastp_bases_weight_q30(analysis_cfg):
    report(analysis_cfg, "fastp1", before=100, after=80, bases=800, q30=720)
    report(analysis_cfg, "fastp2", before=300, after=120, bases=2400, q30=1200)
    cutadapt_report(analysis_cfg, "cutadapt1", before=600, after=300, bases=6000)
    snapshot = analysis_qc.collect(analysis_cfg, [PROJECT])
    filtering = snapshot["projects"][0]["filtering"]
    assert filtering["input_reads"] == 1000
    assert filtering["output_reads"] == 500
    assert filtering["retention_pct"] == 50
    assert filtering["output_bases"] == 9200
    assert filtering["bases_reports_used"] == 3
    assert filtering["q30_pct"] == 60
    assert filtering["q30_reports_used"] == 2
    assert filtering["reports_used"] == filtering["reports_found"] == 3
    assert "fastp: 2/2, cutadapt: 1/1" in snapshot["summary_text"]
    assert "Q30 covers only reports" in snapshot["summary_text"]
    # Schema-v1 consumers still see fastp-only totals, including its own coverage.
    fastp = snapshot["projects"][0]["fastp"]
    assert fastp["input_reads"] == 400
    assert fastp["output_reads"] == 200
    assert fastp["q30_pct"] == 60
    assert fastp["reports_found"] == fastp["reports_used"] == 2


def test_paired_cutadapt_counts_are_normalized_to_individual_reads(analysis_cfg):
    report(analysis_cfg, "fastp1", before=100, after=80)
    cutadapt_report(analysis_cfg, "paired", before=100, after=20, bases=4000, paired=True)
    filtering = analysis_qc.collect(analysis_cfg, [PROJECT])["projects"][0]["filtering"]
    assert filtering["input_reads"] == 300
    assert filtering["output_reads"] == 120
    assert filtering["retention_pct"] == 40
    assert filtering["output_bases"] == 4800
    assert filtering["tools"]["cutadapt"]["input_reads"] == 200


@pytest.mark.parametrize("before, after", [(0, 0), (100, 0)])
def test_cutadapt_zero_reads_remain_zero_and_q30_unavailable(analysis_cfg, before, after):
    cutadapt_report(analysis_cfg, "empty", before=before, after=after, bases=0)
    filtering = analysis_qc.collect(analysis_cfg, [PROJECT])["projects"][0]["filtering"]
    assert filtering["input_reads"] == before
    assert filtering["output_reads"] == filtering["output_bases"] == 0
    assert filtering["retention_pct"] == (0 if before else None)
    assert filtering["q30_pct"] is None
    assert filtering["q30_reports_used"] == 0


@pytest.mark.parametrize(
    "contents",
    [
        "broken",
        "[]",
        "null",
        "{}",
        '{"read_counts":null}',
        '{"read_counts":{"input":100,"output":101},"input":{"paired":false}}',
        '{"read_counts":{"input":-1,"output":0},"input":{"paired":false}}',
        '{"read_counts":{"input":100,"output":true},"input":{"paired":false}}',
        '{"read_counts":{"input":100,"output":40}}',
        '{"read_counts":{"input":100,"output":40},"input":{"paired":"false"}}',
    ],
)
def test_invalid_cutadapt_reports_do_not_hide_valid_metrics(analysis_cfg, contents):
    cutadapt_report(analysis_cfg, "good")
    cutadapt_report(analysis_cfg, "bad").write_text(contents)
    snapshot = analysis_qc.collect(analysis_cfg, [PROJECT])
    filtering = snapshot["projects"][0]["filtering"]
    assert filtering["reports_found"] == 2
    assert filtering["reports_used"] == 1
    assert filtering["input_reads"] == 100
    assert filtering["output_reads"] == 40
    assert "cutadapt: 1/2" in snapshot["summary_text"]
    assert "totals cover usable reports only" in snapshot["summary_text"]


@pytest.mark.parametrize("bases", [None, -1, True, "800"])
def test_missing_cutadapt_bases_preserve_read_metrics(analysis_cfg, bases):
    cutadapt_report(analysis_cfg, "sample1", bases=bases)
    snapshot = analysis_qc.collect(analysis_cfg, [PROJECT])
    filtering = snapshot["projects"][0]["filtering"]
    assert filtering["reports_used"] == 1
    assert filtering["output_reads"] == 40
    assert filtering["output_bases"] is None
    assert filtering["bases_reports_used"] == 0
    assert filtering["q30_pct"] is None
    assert "Retained bases cover only reports" in snapshot["summary_text"]


def test_cutadapt_delivery_symlinks_are_parsed_by_delivery_name(analysis_cfg):
    path = cutadapt_report(analysis_cfg, "sample1")
    # The native target has no .cutadapt.json suffix, as in the workflow delivery rule.
    target = path.parents[2] / "filter" / "sample1.json"
    target.parent.mkdir()
    path.rename(target)
    path.symlink_to(target)
    (path.parent / "alias.cutadapt.json").symlink_to(target)
    provenance = path.parents[2] / "provenance/cutadapt"
    provenance.mkdir(parents=True)
    (provenance / "sample1.adapter.cutadapt.json").write_text(target.read_text())
    filtering = analysis_qc.collect(analysis_cfg, [PROJECT])["projects"][0]["filtering"]
    assert filtering["reports_found"] == filtering["reports_used"] == 1
    assert filtering["input_reads"] == 100


def test_two_tools_for_one_sample_are_excluded_instead_of_double_counted(analysis_cfg):
    cutadapt_report(analysis_cfg, "ambiguous")
    report(analysis_cfg, "ambiguous")
    cutadapt_report(analysis_cfg, "good")
    snapshot = analysis_qc.collect(analysis_cfg, [PROJECT])
    filtering = snapshot["projects"][0]["filtering"]
    assert filtering["reports_found"] == 3
    assert filtering["reports_used"] == 1
    assert filtering["input_reads"] == 100
    assert filtering["output_reads"] == 40
    assert "Multiple filtering reports name the same sample" in snapshot["summary_text"]
