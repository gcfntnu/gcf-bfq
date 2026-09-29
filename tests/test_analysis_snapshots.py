"""Snapshot contents, atomic publication, and recovery against real state/files."""

import json
import os
import tarfile

from unittest.mock import Mock

import pytest

from bcl2fastq_pipeline.state import FlowcellStateStore, new_state
from test_state_integration import RUN_ID, configured_bfq

from bcl2fastq_pipeline import analysis_snapshots as snapshots

PROJECT = "GCF-2026-001"


def put(root, name, content="original"):
    path = root / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)
    return path


def ready(tmp_path, monkeypatch, projects=(PROJECT,)):
    cfg, source, output = configured_bfq(tmp_path)
    output.mkdir()
    monkeypatch.setenv("TMPDIR", str(tmp_path / "work"))
    store = FlowcellStateStore(cfg.static.paths.manager_dir)
    store.create(new_state(RUN_ID, source, output, origin="new", start_stage="analysis"))
    store.begin_attempt(RUN_ID)
    workdirs = {}
    for project in projects:
        work = snapshots.workdir_path(RUN_ID, project)
        put(work, "config.yaml")
        put(work, "Snakefile")
        workdirs[project] = snapshots.identify_workdir(work, RUN_ID, project)
    store.complete_stage(RUN_ID, "analysis", {"workdirs": workdirs})
    store.start_stage(RUN_ID, "reporting")
    store.complete_stage(RUN_ID, "reporting")
    store.start_stage(RUN_ID, "finalization")
    return cfg, store, output


def commit(store, projects=(PROJECT,)):
    results = snapshots.prepare(store.read(RUN_ID), projects)
    store.complete_run(RUN_ID, list(projects), snapshots=results)
    snapshots.recover(store, RUN_ID)
    return results


def archive_path(output, project=PROJECT):
    return output / "provenance" / f"{project}_analysis.tar.gz"


def read_member(output, member="config.yaml", project=PROJECT):
    with tarfile.open(archive_path(output, project)) as archive:
        return archive.extractfile(member).read().decode()


def reanalyze(store):
    store.set_preparing(RUN_ID, "analysis", reason="test", refresh_inputs=False)
    store.queue(RUN_ID, "analysis")
    store.begin_attempt(RUN_ID)
    work = snapshots.workdir_path(RUN_ID, PROJECT)
    record = snapshots.identify_workdir(work, RUN_ID, PROJECT)
    store.complete_stage(RUN_ID, "analysis", {"workdirs": {PROJECT: record}})
    store.start_stage(RUN_ID, "reporting")
    store.complete_stage(RUN_ID, "reporting")
    store.start_stage(RUN_ID, "finalization")


def test_actual_tree_includes_hidden_dirty_and_unrelated_files_without_data(tmp_path, monkeypatch):
    _cfg, store, output = ready(tmp_path, monkeypatch)
    work = snapshots.workdir_path(RUN_ID, PROJECT)
    put(work, "src/gcf-workflows/libprep.config", "phred: 20 # local edit")
    put(work, "src/gcf-workflows/untracked.rules", "untracked")
    put(work, "src/gcf-workflows/.git/HEAD", "ref: refs/heads/local")
    put(work, ".snakemake/log/run.log", "actual execution")
    put(work, "new-unrelated-file", "included automatically")
    put(work, "nested/data/keep.txt", "nested data is kept")
    put(work, "data/large.fastq", "excluded")
    external = put(tmp_path, "external/secret.fastq", "external bytes")
    (work / "external-dir").symlink_to(external.parent, target_is_directory=True)
    (work / "external-file").symlink_to(external)
    (work / "data-link").symlink_to("data", target_is_directory=True)
    (work / "broken-link").symlink_to("missing")
    os.link(work / "config.yaml", work / "hardlink")

    commit(store)

    with tarfile.open(archive_path(output)) as archive:
        members = {member.name: member for member in archive.getmembers()}
        assert "data" not in members
        assert not any(name.startswith("data/") for name in members)
        for name in ("external-dir", "external-file", "data-link", "broken-link"):
            assert members[name].issym()
            assert not any(key.startswith(name + "/") for key in members)
        assert members["hardlink"].islnk()
        for name in (
            "config.yaml",
            "Snakefile",
            "src/gcf-workflows/untracked.rules",
            "src/gcf-workflows/.git/HEAD",
            ".snakemake/log/run.log",
            "new-unrelated-file",
            "nested/data/keep.txt",
            snapshots.MARKER,
        ):
            assert name in members
        assert (
            archive.extractfile("src/gcf-workflows/libprep.config").read()
            == b"phred: 20 # local edit"
        )
    assert not list((output / "provenance/.staging").glob("*"))
    assert "pending" not in store.read(RUN_ID)["analysis_snapshots"][PROJECT]


@pytest.mark.parametrize("error", [OSError("disk full"), KeyboardInterrupt()])
def test_failed_or_interrupted_write_preserves_previous_snapshot(tmp_path, monkeypatch, error):
    _cfg, store, output = ready(tmp_path, monkeypatch)
    commit(store)
    before = archive_path(output).read_bytes()
    reanalyze(store)
    monkeypatch.setattr(tarfile.TarFile, "add", Mock(side_effect=error))
    with pytest.raises(type(error)):
        snapshots.prepare(store.read(RUN_ID), [PROJECT])
    assert archive_path(output).read_bytes() == before
    assert not list((output / "provenance/.staging").glob("*"))


def test_crash_before_completion_discards_candidate_and_preserves_retained(tmp_path, monkeypatch):
    _cfg, store, output = ready(tmp_path, monkeypatch)
    commit(store)
    before = archive_path(output).read_bytes()
    reanalyze(store)
    put(snapshots.workdir_path(RUN_ID, PROJECT), "config.yaml", "unsuccessful")
    snapshots.prepare(store.read(RUN_ID), [PROJECT])
    snapshots.recover_pending(store)  # Simulated fresh daemon, no completion commit.
    assert archive_path(output).read_bytes() == before
    assert not list((output / "provenance/.staging").glob("*"))


def test_committed_replacement_recovers_after_crash_and_rename_failure(tmp_path, monkeypatch):
    _cfg, store, output = ready(tmp_path, monkeypatch)
    commit(store)
    before = archive_path(output).read_bytes()
    reanalyze(store)
    put(snapshots.workdir_path(RUN_ID, PROJECT), "config.yaml", "new delivered analysis")
    results = snapshots.prepare(store.read(RUN_ID), [PROJECT])
    store.complete_run(RUN_ID, [PROJECT], snapshots=results)
    with monkeypatch.context() as patch:
        patch.setattr(snapshots.os, "replace", Mock(side_effect=OSError("rename failed")))
        with pytest.raises(OSError, match="rename failed"):
            snapshots.recover(store, RUN_ID)
    assert archive_path(output).read_bytes() == before
    assert (output / results[PROJECT]["pending"]).exists()
    snapshots.recover_pending(FlowcellStateStore(store.manager_dir))
    assert read_member(output) == "new delivered analysis"
    assert len(list((output / "provenance").glob("*.tar.gz"))) == 1
    assert not list((output / "provenance/.staging").glob("*"))


def test_crash_after_rename_is_idempotent(tmp_path, monkeypatch):
    _cfg, store, output = ready(tmp_path, monkeypatch)
    results = snapshots.prepare(store.read(RUN_ID), [PROJECT])
    store.complete_run(RUN_ID, [PROJECT], snapshots=results)
    with monkeypatch.context() as patch:
        patch.setattr(store, "mutate", Mock(side_effect=KeyboardInterrupt()))
        with pytest.raises(KeyboardInterrupt):
            snapshots.recover(store, RUN_ID)
    assert read_member(output) == "original"
    before = archive_path(output).stat().st_mtime_ns
    snapshots.recover(store, RUN_ID)
    assert archive_path(output).stat().st_mtime_ns == before
    assert "pending" not in store.read(RUN_ID)["analysis_snapshots"][PROJECT]


def test_all_projects_prepare_before_any_snapshot_is_replaced(tmp_path, monkeypatch):
    projects = (PROJECT, "GCF-2026-002")
    _cfg, store, output = ready(tmp_path, monkeypatch, projects)
    # First project's candidate succeeds; second fails.
    real_write = snapshots._write_archive
    calls = []

    def write(*args):
        calls.append(args)
        if len(calls) == 2:
            raise OSError("second project failed")
        return real_write(*args)

    monkeypatch.setattr(snapshots, "_write_archive", write)
    with pytest.raises(OSError, match="second project"):
        snapshots.prepare(store.read(RUN_ID), projects)
    assert not list((output / "provenance").rglob("*.tar.gz"))


@pytest.mark.parametrize("problem", ["missing", "wrong-owner", "missing-config"])
def test_unavailable_is_explicit_and_does_not_replace_older_snapshot(
    tmp_path, monkeypatch, problem
):
    _cfg, store, output = ready(tmp_path, monkeypatch)
    commit(store)
    before = archive_path(output).read_bytes()
    reanalyze(store)
    work = snapshots.workdir_path(RUN_ID, PROJECT)
    if problem == "missing":
        work.rename(work.with_name("moved"))
    elif problem == "wrong-owner":
        snapshots.identify_workdir(work, "another-flowcell", PROJECT)
    else:
        (work / "config.yaml").unlink()
    results = commit(store)
    assert results[PROJECT]["status"] == "unavailable"
    assert archive_path(output).read_bytes() == before
    assert store.read(RUN_ID)["stages"]["finalization"]["metadata"]["analysis_snapshots"] == results


def test_recorded_workdir_survives_tmpdir_change(tmp_path, monkeypatch):
    _cfg, store, output = ready(tmp_path, monkeypatch)
    monkeypatch.setenv("TMPDIR", str(tmp_path / "different"))
    commit(store)
    assert read_member(output) == "original"


def test_legacy_workdir_can_be_retained_without_reconstructing_it(tmp_path, monkeypatch, caplog):
    _cfg, store, output = ready(tmp_path, monkeypatch)
    work = snapshots.workdir_path(RUN_ID, PROJECT)
    (work / snapshots.MARKER).unlink()
    state = store.read(RUN_ID)
    del state["stages"]["analysis"]["metadata"]["workdirs"]
    store.write(state)
    commit(store)
    assert read_member(output) == "original"
    assert "legacy workdir without an analysis ownership token" in caplog.text


def test_top_level_data_symlink_is_excluded(tmp_path, monkeypatch):
    _cfg, store, output = ready(tmp_path, monkeypatch)
    work = snapshots.workdir_path(RUN_ID, PROJECT)
    external = put(tmp_path, "external/big.fastq")
    (work / "data").symlink_to(external.parent, target_is_directory=True)
    commit(store)
    with tarfile.open(archive_path(output)) as archive:
        assert "data" not in archive.getnames()
    assert json.loads((work / snapshots.MARKER).read_text())["run_id"] == RUN_ID


def test_workdir_reused_during_archive_fails_instead_of_publishing_mixed_tree(
    tmp_path, monkeypatch
):
    _cfg, store, output = ready(tmp_path, monkeypatch)
    real_write = snapshots._write_archive

    def reused(workdir, *args):
        result = real_write(workdir, *args)
        snapshots.identify_workdir(workdir, "other-run", PROJECT)
        return result

    monkeypatch.setattr(snapshots, "_write_archive", reused)
    with pytest.raises(RuntimeError, match="changed while archiving"):
        snapshots.prepare(store.read(RUN_ID), [PROJECT])
    assert not list((output / "provenance").rglob("*.tar.gz"))


def test_multi_project_publication_recovers_after_one_project_was_published(tmp_path, monkeypatch):
    projects = (PROJECT, "GCF-2026-002")
    _cfg, store, output = ready(tmp_path, monkeypatch, projects)
    results = snapshots.prepare(store.read(RUN_ID), projects)
    store.complete_run(RUN_ID, list(projects), snapshots=results)
    replace = snapshots.os.replace

    def fail_second(src, dst):
        if str(dst).endswith("GCF-2026-002_analysis.tar.gz"):
            raise OSError("second publication failed")
        return replace(src, dst)

    with monkeypatch.context() as patch:
        patch.setattr(snapshots.os, "replace", fail_second)
        with pytest.raises(OSError, match="second publication"):
            snapshots.recover(store, RUN_ID)
    assert read_member(output) == "original"
    assert not archive_path(output, projects[1]).exists()
    snapshots.recover(store, RUN_ID)
    assert read_member(output, project=projects[1]) == "original"
    assert len(list((output / "provenance").glob("*.tar.gz"))) == 2


@pytest.mark.parametrize("directory", ["provenance", "provenance/.staging"])
def test_snapshot_storage_rejects_symlink_directory(tmp_path, monkeypatch, directory):
    _cfg, store, output = ready(tmp_path, monkeypatch)
    external = tmp_path / "external"
    external.mkdir()
    link = output / directory
    link.parent.mkdir(exist_ok=True)
    link.symlink_to(external, target_is_directory=True)
    with pytest.raises(RuntimeError, match="must not be a symlink"):
        snapshots.prepare(store.read(RUN_ID), [PROJECT])
    assert not list(external.iterdir())
