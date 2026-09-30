"""Sequencing report fixtures exercise inputs without analysis or real SMTP."""

import json
import subprocess

from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml

from bcl2fastq_pipeline import sequencing_qc as qc

RUN_INFO = """<RunInfo><Run Id="260925_NB501038_0281_AHL2T7AFXC"><Reads>
<Read Number="1" NumCycles="151" IsIndexedRead="N"/>
<Read Number="2" NumCycles="8" IsIndexedRead="Y"/>
<Read Number="3" NumCycles="8" IsIndexedRead="Y"/>
<Read Number="4" NumCycles="151" IsIndexedRead="N"/>
</Reads></Run></RunInfo>"""
INTEROP = """# Version: v1.5.0
Level,Yield,Projected Yield,Aligned,Error Rate,Intensity C1,%>=Q30
Read 1,15.1,15.1,1.2,0.2,200,94
Read 2 (I),0.8,0.8,0,0,100,90
Read 3 (I),0.8,0.8,0,0,100,90
Read 4,15.1,15.1,1.3,0.2,200,93
Non-indexed,30.2,30.2,1.2,0.2,200,93.5
Total,31.8,31.8,1.2,0.2,200,93.2

Read 1
Lane,Surface,Tiles,Density,Cluster PF,Reads,Aligned,%>=Q30
1,-,10,1000 +/- 20,85.5 +/- 1,100,1.2 +/- 0.1,94
1,1,5,900 +/- 20,85.5 +/- 1,50,1.2 +/- 0.1,20
Read 2 (I)
Lane,Surface,Tiles,Density,Cluster PF,Reads,Aligned,%>=Q30
1,-,10,1000 +/- 20,85.5 +/- 1,100,0,90
Read 3 (I)
Lane,Surface,Tiles,Density,Cluster PF,Reads,Aligned,%>=Q30
1,-,10,1000 +/- 20,85.5 +/- 1,100,0,90
Read 4
Lane,Surface,Tiles,Density,Cluster PF,Reads,Aligned,%>=Q30
1,-,10,1000 +/- 20,85.5 +/- 1,100,1.3 +/- 0.1,93
Extracted: 318
"""


@pytest.fixture
def cfg(tmp_path):
    output, instrument = tmp_path / "output", tmp_path / "instrument"
    output.mkdir()
    instrument.mkdir()
    (instrument / "RunInfo.xml").write_text(RUN_INFO)
    (instrument / "runParameters.xml").write_text("<RunParameters/>")
    (output / "InterOp").mkdir()
    (output / "InterOp" / "IndexMetricsOut.bin").touch()
    sheet = tmp_path / "SampleSheet.csv"
    sheet.write_text(
        "[Header]\nInvestigator Name,Operator\n[Data]\nSample_ID,Sample_Name,Sample_Project,Lane\ns1,first,P1,1\ns2,second,P2,1\ns3,third,P2,1\n"
    )
    workbook = tmp_path / "bad.xlsx"
    workbook.write_text("This workbook must never be read for sequencing QC")
    return SimpleNamespace(
        output_path=output,
        static=SimpleNamespace(commands={"multiqc_options": "-q --interactive"}),
        run=SimpleNamespace(
            run_id="260925_NB501038_0281_AHL2T7AFXC",
            flowcell_path=instrument,
            sample_sheet=sheet,
            sample_submission_form=workbook,
        ),
    )


@pytest.fixture
def report_commands(monkeypatch):
    calls = []

    def interop(command, run_path, output_path, cwd):
        calls.append((command, run_path, output_path, cwd))
        output_path.write_text(INTEROP if command == "interop_summary" else "index fixture")

    def multiqc(command, cwd, stdout, stderr):
        calls.append(command)
        (
            Path(command[command.index("--outdir") + 1]) / command[command.index("--filename") + 1]
        ).write_text("<html>Sequencing QC</html>")

    monkeypatch.setattr(qc, "run_interop_csv", interop)
    monkeypatch.setattr(qc.subprocess, "check_call", multiqc)
    return calls


def write_convert(cfg, *, counts=(1, 0, 9999)):
    reports = cfg.output_path / "Reports"
    reports.mkdir(exist_ok=True)
    (reports / "Demultiplex_Stats.csv").write_text(
        "Lane,SampleID,Index,# Reads,# Perfect Index Reads,# One Mismatch Index Reads\n"
        f"1,s1,AAAA+CCCC,{counts[0]},{counts[0]},0\n1,s2,GGGG+TTTT,{counts[1]},{counts[1]},0\n"
        f"1,Undetermined,,{counts[2]},0,0\n"
    )
    return reports


def write_bcl2fastq(cfg, *, root=None):
    root = root or cfg.output_path / "Stats"
    root.mkdir(parents=True, exist_ok=True)
    (root / "Stats.json").write_text(
        json.dumps(
            {
                "RunId": cfg.run.run_id,
                "ConversionResults": [
                    {
                        "LaneNumber": 1,
                        "TotalClustersPF": 100,
                        "DemuxResults": [
                            {"SampleId": "s1", "SampleName": "first", "NumberReads": 75},
                            {"SampleId": "s2", "SampleName": "second", "NumberReads": 0},
                        ],
                        "Undetermined": {"NumberReads": 25},
                    }
                ],
                "UnknownBarcodes": [{"Lane": 1, "Barcodes": {"ATGC+ATGC": 20}}],
            }
        )
    )
    return root


def test_convert_report_without_analysis_or_fastqs(cfg, report_commands):
    reports = write_convert(cfg)
    (reports / "Top_Unknown_Barcodes.csv").write_text(
        "Lane,index,index2,# Reads\n1,AAAA,GGGG,9000\n"
    )
    result = qc.generate(cfg, tool=["bcl-convert", "4.2.4"])

    assert result["projects"] == ["P1", "P2"]
    assert result["total_reads"] == 10000
    assert result["undetermined_percent"] == 99.99
    assert result["zero_read_sample_count"] == 1
    assert result["unavailable_sample_count"] == 1
    assert [sample["assigned_reads"] for sample in result["samples"]] == [1, 0, None]
    assert result["read_geometry"] == "R151 / I8 / I8 / R151"
    assert result["yield_gb"] == 30.2
    lane = result["lane_metrics"][0]
    assert (lane["density_k_mm2"], lane["pf_percent"], lane["phix_percent"]) == (1000, 85.5, 1.2)
    assert (lane["R1_q30_percent"], lane["R2_q30_percent"]) == (94, 93)
    assert "AAAA+GGGG: 9,000" in result["summary_text"]
    assert not (cfg.output_path / "P1").exists()

    command = report_commands[-1]
    assert command[command.index("--filename") + 1] == "sequencing_qc.html"
    assert command[-2] == "bclconvert"
    inputs = Path(command[-1])
    assert inputs.parent == cfg.output_path / "Stats" / "sequencing_qc"
    assert (inputs / "RunInfo.xml").is_file()
    assert (inputs / "Demultiplex_Stats.csv").is_file()
    config = yaml.safe_load((inputs.parent / "multiqc_config.yaml").read_text())
    assert config["megaqc_upload"] is False
    assert config["megaqc_access_token"] is None
    assert "--no-megaqc-upload" in command
    samples = config["custom_data"]["bfq_sample_assignment"]["data"]
    assert samples["P2 / s2"]["assigned_reads"] == 0
    assert samples["P2 / s3"]["assigned_reads"] == "Unavailable"
    assert json.loads((inputs.parent / "summary.json").read_text()) == result


@pytest.mark.parametrize(
    "tool,nested",
    [
        ("bcl2fastq", False),
        ("cellranger mkfastq", True),
        ("spaceranger mkfastq", True),
        ("cellranger-atac mkfastq", True),
    ],
)
def test_bcl2fastq_and_10x_layouts(cfg, report_commands, tool, nested):
    root = cfg.output_path / "outs" / "fastq_path" / "Stats" if nested else None
    write_bcl2fastq(cfg, root=root)
    result = qc.generate(cfg, tool=tool)
    assert result["demultiplexer"] == "bcl2fastq"
    assert result["total_reads"] == 100
    assert result["undetermined_percent"] == 25
    assert result["unknown_indexes"][0]["reads"] == 20
    assert result["zero_read_sample_count"] == 1
    assert report_commands[-1][-2] == "bcl2fastq"


def test_source_identification_uses_real_outputs_without_environment(
    cfg, report_commands, monkeypatch
):
    write_bcl2fastq(cfg)
    monkeypatch.setenv("FORCE_BCL2FASTQ", "")
    assert qc.generate(cfg)["demultiplexer"] == "bcl2fastq"
    write_convert(cfg)
    assert qc.generate(cfg)["demultiplexer"] == "bclconvert"
    # A known mkfastq/bcl2fastq execution prefers its own native output.
    assert qc.generate(cfg, tool="cellranger mkfastq")["demultiplexer"] == "bcl2fastq"


def test_xml_fallback_preserves_zero_and_avoids_all_barcode_double_count(cfg, report_commands):
    stats = cfg.output_path / "Stats"
    stats.mkdir()
    (stats / "DemultiplexingStats.xml").write_text("""<Stats><Flowcell>
<Project name="P1"><Sample name="s1"><Barcode name="all"><Lane number="1"><BarcodeCount>10</BarcodeCount></Lane></Barcode><Barcode name="AAAA"><Lane number="1"><BarcodeCount>10</BarcodeCount></Lane></Barcode></Sample></Project>
<Project name="P2"><Sample name="s2"><Barcode name="all"><Lane number="1"><BarcodeCount>0</BarcodeCount></Lane></Barcode></Sample></Project>
<Project name="all"><Sample name="all"><Barcode name="all"><Lane number="1"><BarcodeCount>100</BarcodeCount></Lane></Barcode></Sample></Project>
<Project name="default"><Sample name="all"><Barcode name="all"><Lane number="1"><BarcodeCount>90</BarcodeCount></Lane></Barcode></Sample></Project>
</Flowcell></Stats>""")
    result = qc.generate(cfg)
    assert result["samples"][0]["assigned_reads"] == 10
    assert result["undetermined_percent"] == 90
    assert result["zero_read_sample_count"] == 1
    assert result["unknown_indexes"] is None
    assert "XML assignment" in result["summary_text"]


def test_optional_interop_failure_and_missing_counts_are_unavailable(
    cfg, report_commands, monkeypatch
):
    reports = write_convert(cfg)
    (reports / "Demultiplex_Stats.csv").write_text("Lane,SampleID,# Reads\n1,s1,10\n")

    def failed_interop(*args):
        raise RuntimeError("missing tile metrics")

    monkeypatch.setattr(qc, "run_interop_csv", failed_interop)
    result = qc.generate(cfg)
    assert result["undetermined_reads"] is None
    assert result["undetermined_percent"] is None
    assert "Non-index yield: Unavailable" in result["summary_text"]
    assert "missing tile metrics" in result["summary_text"]
    assert "Unknown-index information unavailable" in result["summary_text"]
    json.dumps(result, allow_nan=False)


def test_bclconvert_interop_index_compatibility_link(cfg, report_commands):
    reports = write_convert(cfg)
    (cfg.output_path / "InterOp" / "IndexMetricsOut.bin").unlink()
    (reports / "IndexMetricsOut.bin").write_bytes(b"mock metrics")
    qc.generate(cfg)
    assert (
        cfg.output_path / "InterOp" / "IndexMetricsOut.bin"
    ).resolve() == reports / "IndexMetricsOut.bin"


def test_multilane_counts_deduplicate_planned_samples(cfg, report_commands):
    reports = write_convert(cfg)
    cfg.run.sample_sheet.write_text(
        "[BCLConvert_Data]\nLane,Sample_ID,Sample_Project\n1,s1,P1\n2,s1,P1\n"
    )
    (reports / "Demultiplex_Stats.csv").write_text(
        "Lane,SampleID,# Reads\n1,s1,10\n1,Undetermined,5\n2,s1,20\n2,Undetermined,15\n"
    )
    result = qc.generate(cfg)
    assert len(result["samples"]) == 1
    assert result["samples"][0]["assigned_reads"] == 30
    assert result["undetermined_percent"] == 40
    assert "P1: 1 unique planned samples" in result["summary_text"]


def test_no_samplesheet_uses_project_identifiers_in_demux_stats(cfg, report_commands):
    reports = write_convert(cfg)
    cfg.run.sample_sheet = None
    (reports / "Demultiplex_Stats.csv").write_text(
        "Lane,SampleID,Sample_Project,# Reads\n1,s1,P-from-stats,0\n1,Undetermined,,100\n"
    )
    result = qc.generate(cfg)
    assert result["projects"] == ["P-from-stats"]
    assert result["samples"][0]["assigned_reads"] == 0
    assert "planned samples cannot be enumerated" in result["summary_text"]


def test_plain_10x_sheet_and_html_escaping(cfg, report_commands):
    write_bcl2fastq(cfg)
    cfg.run.sample_sheet.write_text('Lane,Sample,Index,Sample_Project\n1,s1,SI-001,"P<script>"\n')
    result = qc.generate(cfg, tool="cellranger mkfastq")
    assert "P&lt;script&gt;" in result["summary_html"]
    assert "P<script>" not in result["summary_html"]


def test_malformed_demux_is_diagnostic_and_does_not_publish_report(cfg, report_commands):
    reports = write_convert(cfg)
    (reports / "Demultiplex_Stats.csv").write_text("Lane,SampleID,# Reads\n1,s1,not-a-number\n")
    with pytest.raises(
        RuntimeError, match="Unable to read demultiplexer statistics.*Invalid read count"
    ):
        qc.generate(cfg)
    assert not (cfg.output_path / qc.REPORT_DIRECTORY / qc.REPORT_FILENAME).exists()


def test_multiqc_failure_has_log_tail(cfg, report_commands, monkeypatch):
    write_convert(cfg)

    def failed(command, cwd, stdout, stderr):
        stdout.write("Missing required MultiQC module\n")
        raise subprocess.CalledProcessError(2, command)

    monkeypatch.setattr(qc.subprocess, "check_call", failed)
    with pytest.raises(
        RuntimeError,
        match="(?s)Sequencing QC generation failed.*multiqc.log.*Missing required MultiQC module",
    ):
        qc.generate(cfg)


def test_empty_multiqc_success_does_not_accept_stale_html(cfg, report_commands, monkeypatch):
    write_convert(cfg)
    first = qc.generate(cfg)
    monkeypatch.setattr(qc.subprocess, "check_call", lambda *args, **kwargs: None)
    with pytest.raises(RuntimeError, match="without a nonempty HTML report"):
        qc.generate(cfg)
    assert not Path(first["report_path"]).exists()
