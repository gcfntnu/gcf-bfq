"""Project metadata and conversion isolation with synthetic inputs only."""

import csv
import hashlib
import io
import json

from pathlib import Path
from types import SimpleNamespace

import pytest

from bcl2fastq_pipeline.config import ProjectConfig
from test_state_integration import configured_bfq, write_fastq

from bcl2fastq_pipeline import afterFastq, demultiplexing, makeFastq, workflow_config

PROJECTS = ("GCF-2026-101", "GCF-2026-102")


def pair(  # noqa: PLR0913
    tmp_path, project, *, kit="Kit A", settings="", data_settings="", value="", index="AAAAAAAA"
):
    tmp_path.mkdir(parents=True, exist_ok=True)
    sheet = tmp_path / f"SampleSheet_{project}.csv"
    sheet.write_text(
        f"[CustomOptions]\nLibprep,{kit}\nOrganism,organism-{project}\n"
        f"[Settings]\n{settings}[Data]\nSample_ID,Sample_Project,index{data_settings}\n"
        f"sample-{project},{project},{index}{value}\n"
    )
    form = tmp_path / f"Submission-Form_{project}.xlsx"
    form.write_bytes(b"synthetic metadata only")
    return SimpleNamespace(project_id=project, sample_sheet=sheet, submission_form=form)


def test_different_kits_share_one_structured_pass_with_per_sample_adapters(tmp_path):
    first = pair(tmp_path, PROJECTS[0], settings="AdapterRead1,ACGT\n")
    second = pair(tmp_path, PROJECTS[1], kit="Kit B", index="CCCCCCCC")
    original = [p.sample_sheet.read_bytes() for p in (first, second)]
    result = demultiplexing.plan(SimpleNamespace(pairs=(second, first)))
    assert result["passes"] == 1
    assert result["tool"] == "bcl-convert"
    assert result["per_sample_settings"] == ["AdapterRead1"]
    assert "FileFormatVersion,2" in result["sample_sheet_content"]
    rows = list(csv.reader(io.StringIO(result["sample_sheet_content"])))
    header = rows[rows.index(["[BCLConvert_Data]"]) + 1]
    assert [row[-1] for row in rows[-2:]] == ["ACGT", "na"]
    assert header[-1] == "AdapterRead1"
    assert [row["source_row"] for row in result["rows"]] == [8, 7]
    assert [p.sample_sheet.read_bytes() for p in (first, second)] == original
    assert all(len(source["sample_sheet_sha256"]) == 64 for source in result["sources"])


@pytest.mark.parametrize("setting", ["MaskShortReads", "MinimumTrimmedReadLength", "TrimUMI"])
def test_unsupported_differences_fail_before_any_outputs(tmp_path, setting):
    first = pair(tmp_path, PROJECTS[0], settings=f"{setting},10\n")
    second = pair(tmp_path, PROJECTS[1], settings=f"{setting},20\n", index="CCCCCCCC")
    with pytest.raises(demultiplexing.DemultiplexingPlanError, match=f"different {setting}"):
        demultiplexing.plan(SimpleNamespace(pairs=(first, second)))
    assert sorted(p.name for p in tmp_path.iterdir()) == sorted(
        p.name for item in (first, second) for p in (item.sample_sheet, item.submission_form)
    )


def test_row_specific_cycle_settings_remain_bound_to_source_rows(tmp_path):
    first = pair(tmp_path, PROJECTS[0], data_settings=",OverrideCycles", value=",Y75;I8;Y75")
    second = pair(tmp_path, PROJECTS[1], settings="OverrideCycles,Y75;I8;N75\n", index="CCCCCCCC")
    result = demultiplexing.plan(SimpleNamespace(pairs=(first, second)))
    assert [row["effective_parameters"]["OverrideCycles"] for row in result["rows"]] == [
        "Y75;I8;Y75",
        "Y75;I8;N75",
    ]
    assert "OverrideCycles" not in result["global_settings"]
    with pytest.raises(demultiplexing.DemultiplexingPlanError, match="FORCE_BCL2FASTQ"):
        demultiplexing.plan(SimpleNamespace(pairs=(first, second)), force_bcl2fastq=True)


def test_unsupported_backend_is_actionable(tmp_path):
    first = pair(tmp_path, PROJECTS[0], kit="10X Genomics Chromium kit")
    second = pair(tmp_path, PROJECTS[1], index="CCCCCCCC")
    with pytest.raises(demultiplexing.DemultiplexingPlanError, match="BFQ #156"):
        demultiplexing.plan(SimpleNamespace(pairs=(first, second)))


def test_conversion_uses_derived_sheet_and_retains_its_exact_inputs(tmp_path, monkeypatch):
    cfg, source, output = configured_bfq(tmp_path)
    (source / "InterOp").mkdir()
    cfg.run.input_selection = SimpleNamespace(pairs=tuple(pair(output, p) for p in PROJECTS))
    cfg.run.libprep = None
    calls = []
    monkeypatch.setattr(makeFastq, "image_version", lambda _tool: "4.2.4")
    monkeypatch.setattr(makeFastq.subprocess, "check_call", lambda cmd, **kw: calls.append(cmd))
    assert makeFastq.bcl2fq() == ["bcl-convert", "4.2.4"]
    assert calls[0][calls[0].index("--sample-sheet") + 1] == str(cfg.run.sample_sheet)
    provenance = output / "provenance/demultiplexing"
    data = json.loads((provenance / "plan.json").read_text())
    assert data["status"] == "completed"
    assert data["tool_version"] == "4.2.4"
    for p in cfg.run.input_selection.pairs:
        assert (
            provenance / "sources" / p.project_id / p.sample_sheet.name
        ).read_bytes() == p.sample_sheet.read_bytes()


def test_different_project_workflows_receive_independent_pairs_and_snapshot(tmp_path, monkeypatch):
    cfg, _source, output = configured_bfq(tmp_path)
    authoritative = tmp_path / "workflows/libprep.config"
    authoritative.parent.mkdir()
    authoritative.write_text("Kit A PE:\n  workflow: rnaseq\nKit B SE:\n  workflow: default\n")
    monkeypatch.setattr(workflow_config, "AUTHORITATIVE_CONFIG", authoritative)
    monkeypatch.setenv("TMPDIR", str(tmp_path / "analysis"))
    monkeypatch.setenv("SINGULARITY_CACHEDIR", str(tmp_path / "cache"))
    for project, kit, geometry in zip(PROJECTS, ("Kit A", "Kit B"), ((75, 75), (75,)), strict=True):
        selected = pair(output, project, kit=kit)
        cfg.run.project_configs[project] = ProjectConfig(
            project,
            selected.sample_sheet,
            selected.submission_form,
            {"Libprep": kit, "Organism": f"organism-{project}"},
            geometry,
        )
        write_fastq(output / project / f"sample-{project}_R1.fastq.gz")
    commands = []

    def configmaker(cmd, cwd):
        project = cmd[cmd.index("-p") + 1]
        context = cfg.run.project_configs[project]
        assert cmd[cmd.index("--sample-sheet") + 1] == str(context.sample_sheet)
        assert cmd[cmd.index("--sample-submission-form") + 1] == str(context.submission_form)
        assert context.custom["Organism"] == f"organism-{project}"
        assert (
            Path(cmd[cmd.index("--libprep-config") + 1]).read_bytes()
            == cfg.run.libprep_config.content
        )
        commands.append(cmd)
        # An edit to the authoritative config must not affect the next project.
        authoritative.write_text("invalid: [")

    def snakemake(_cmd, cwd, log_path):
        project = cwd.name.rsplit("_", 1)[0]
        result = cwd / "data/tmp" / cfg.run.project_workflows[project] / "bfq"
        result.mkdir(parents=True)
        (result / f"multiqc_{project}.html").write_text(project)
        (result / ".multiqc_config.yaml").write_text("{}")
        (cwd / "data/tmp/sample_info.tsv").write_text(project)

    monkeypatch.setattr(afterFastq.subprocess, "check_call", configmaker)
    monkeypatch.setattr(afterFastq, "run_logged_command", snakemake)
    workdirs = afterFastq.full_align(cfg)
    assert set(workdirs) == set(PROJECTS)
    assert cfg.run.project_workflows == dict(zip(PROJECTS, ("rnaseq", "default"), strict=True))
    assert len(commands) == 2
    assert all((output / f"multiqc_{p}_260918.html").read_text() == p for p in PROJECTS)


def test_scoped_checksums_preserve_unselected_manifest_bytes_and_timestamp(tmp_path):
    cfg, _source, output = configured_bfq(tmp_path)
    for project in PROJECTS:
        write_fastq(output / project / "sample_R1.fastq.gz")
    afterFastq.md5sum_worker(cfg)
    retained = output / f"md5sum_{PROJECTS[0]}_fastq.txt"
    before = retained.read_bytes(), retained.stat().st_mtime_ns
    cfg.run.selected_projects = {PROJECTS[1]}
    afterFastq.md5sum_worker(cfg, force=True)
    assert (retained.read_bytes(), retained.stat().st_mtime_ns) == before
    cfg.run.output_override = output / "provenance/correction/working"
    assert workflow_config.project_workdir(cfg, PROJECTS[1]).is_relative_to(
        cfg.run.output_override.parent
    )


def test_delivery_provenance_contains_only_own_project_rows(tmp_path):
    cfg, _source, output = configured_bfq(tmp_path)
    cfg.run.input_selection = SimpleNamespace(pairs=tuple(pair(output, p) for p in PROJECTS))
    demultiplexing.prepare(cfg)
    directory = demultiplexing.delivery_provenance(cfg, PROJECTS[0])
    data = json.loads((directory / "plan.json").read_text())
    assert data["sample_sheet_sha256"] == cfg.run.demultiplexing_plan["sample_sheet_sha256"]
    assert [row["project"] for row in data["rows"]] == [PROJECTS[0]]
    assert [source["project"] for source in data["sources"]] == [PROJECTS[0]]
    assert PROJECTS[1] not in (directory / "SampleSheet.project.csv").read_text()


def test_delivery_inputs_use_completed_snapshot_and_reject_corruption(tmp_path, monkeypatch):
    cfg, _source, output = configured_bfq(tmp_path)
    monkeypatch.setenv("TMPDIR", str(tmp_path / "analysis"))
    selected = pair(output, PROJECTS[0])
    original = selected.sample_sheet.read_bytes()
    captured = workflow_config.project_workdir(cfg, PROJECTS[0]) / "inputs"
    captured.mkdir(parents=True)
    for source in (selected.sample_sheet, selected.submission_form):
        (captured / source.name).write_bytes(source.read_bytes())
    data = {
        "project": PROJECTS[0],
        "sample_sheet": selected.sample_sheet.name,
        "submission_form": selected.submission_form.name,
        "custom": {"Libprep": "Kit A"},
        "inputs": {
            source.name: hashlib.sha256(source.read_bytes()).hexdigest()
            for source in (selected.sample_sheet, selected.submission_form)
        },
    }
    (captured / "project.json").write_text(json.dumps(data))
    selected.sample_sheet.write_text("later unprocessed edits")
    context = workflow_config.delivery_project_config(cfg, PROJECTS[0])
    assert context.sample_sheet.read_bytes() == original
    context.sample_sheet.write_text("broken snapshot")
    with pytest.raises(ValueError, match="captured hash"):
        workflow_config.delivery_project_config(cfg, PROJECTS[0])


def test_project_deliveries_include_only_their_own_sequencing_report(tmp_path, monkeypatch):
    cfg, _source, output = configured_bfq(tmp_path)
    monkeypatch.setenv("TMPDIR", str(tmp_path / "analysis"))
    for project in PROJECTS:
        selected = pair(output, project)
        cfg.run.project_configs[project] = ProjectConfig(
            project, selected.sample_sheet, selected.submission_form, {"Libprep": "Kit A"}
        )
        cfg.run.project_workflows[project] = "default"
        write_fastq(output / project / "sample_R1.fastq.gz")
        (output / f"sequencer_stats_{project}_260918.html").write_text(project)
    commands = []
    monkeypatch.setattr(
        afterFastq, "run_logged_command", lambda command, **kw: commands.append(command)
    )
    afterFastq.archive_worker(cfg)
    for command in commands:
        if "/QC_" in command[2]:
            continue
        project = next(p for p in PROJECTS if Path(command[2]).name.startswith(p))
        reports = [arg for arg in command if Path(arg).name.startswith("sequencer_stats_")]
        assert reports == [str(output / f"sequencer_stats_{project}_260918.html")]


def test_archive_checksum_selection_matches_complete_project_prefix(tmp_path, monkeypatch):
    from unittest.mock import MagicMock  # noqa: PLC0415

    cfg, _source, output = configured_bfq(tmp_path)
    output.mkdir(parents=True)
    cfg.run.selected_projects = {"GCF-2026-01"}
    names = ("GCF-2026-01_260918.7za", "QC_GCF-2026-01_260918.7za", "GCF-2026-011_260918.7za")
    for name in names:
        (output / name).write_bytes(b"synthetic archive")
    pool = MagicMock()
    pool.__enter__.return_value = pool
    monkeypatch.setattr(afterFastq.mp, "Pool", lambda: pool)
    afterFastq.md5sum_archive_worker(cfg)
    assert {path.name for path in pool.map.call_args.args[1]} == set(names[:2])


def test_downstream_provenance_keeps_project_conversion_until_explicit_refresh(tmp_path):
    cfg, _source, output = configured_bfq(tmp_path)
    cfg.run.input_selection = SimpleNamespace(pairs=tuple(pair(output, p) for p in PROJECTS))
    demultiplexing.prepare(cfg)
    first_hash = cfg.run.demultiplexing_plan["sample_sheet_sha256"]
    own = demultiplexing.delivery_provenance(cfg, PROJECTS[1], refresh=True)
    before = (own / "plan.json").read_bytes()
    # A later shared pass may describe different settings while retained FASTQs
    # still belong to the prior per-project conversion.
    sheet = cfg.run.input_selection.pairs[1].sample_sheet
    sheet.write_text(sheet.read_text().replace("[Settings]\n", "[Settings]\nAdapterRead1,ACGT\n"))
    demultiplexing.prepare(cfg)
    assert cfg.run.demultiplexing_plan["sample_sheet_sha256"] != first_hash
    assert demultiplexing.delivery_provenance(cfg, PROJECTS[1]) == own
    assert (own / "plan.json").read_bytes() == before
    demultiplexing.delivery_provenance(cfg, PROJECTS[1], refresh=True)
    assert json.loads((own / "plan.json").read_text())["sample_sheet_sha256"] != first_hash
    (own / "SampleSheet.project.csv").write_text("corrupted")
    with pytest.raises(demultiplexing.DemultiplexingPlanError, match="derived sheet hash"):
        demultiplexing.delivery_provenance(cfg, PROJECTS[1])
