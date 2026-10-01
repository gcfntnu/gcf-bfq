import hashlib

from types import SimpleNamespace
from unittest.mock import Mock

import pytest

from bcl2fastq_pipeline.state import apply_cleanup, cleanup_plan

from bcl2fastq_pipeline import afterFastq


@pytest.fixture
def checksum_run(tmp_path):
    project = tmp_path / "GCF-2026-001"
    project.mkdir()
    for name in ["sample_R1.fastq.gz", "sample R2;$.fastq.gz"]:
        (project / name).write_bytes(name.encode())
    return SimpleNamespace(output_path=tmp_path), tmp_path / "md5sum_GCF-2026-001_fastq.txt"


@pytest.mark.parametrize("stage", ["analysis", "reporting", "finalization"])
def test_downstream_restart_reuses_manifest_without_reading_fastqs(
    checksum_run, monkeypatch, stage
):
    cfg, manifest = checksum_run
    afterFastq.md5sum_worker(cfg)
    original = manifest.read_bytes(), manifest.stat().st_mtime_ns
    apply_cleanup(cleanup_plan(cfg.output_path, stage))
    monkeypatch.setattr(afterFastq, "file_md5", Mock(side_effect=AssertionError("FASTQ read")))
    afterFastq.md5sum_worker(cfg)
    assert (manifest.read_bytes(), manifest.stat().st_mtime_ns) == original


@pytest.mark.parametrize(
    "damage", ["missing", "empty", "partial", "duplicate", "bad_digest", "extra", "truncated"]
)
def test_legacy_manifest_is_repaired_once(checksum_run, monkeypatch, damage):
    cfg, manifest = checksum_run
    afterFastq.md5sum_worker(cfg)
    complete = manifest.read_text()
    damaged = {
        "empty": "",
        "partial": complete.splitlines(keepends=True)[0],
        "duplicate": complete + complete,
        "bad_digest": "x" + complete[1:],
        "extra": complete + "0" * 32 + "  GCF-2026-001/missing.fastq.gz\n",
        "truncated": complete.rstrip("\n"),
    }
    if damage == "missing":
        manifest.unlink()
    else:
        manifest.write_text(damaged[damage])
    afterFastq.md5sum_worker(cfg)
    assert manifest.read_text() == complete
    monkeypatch.setattr(afterFastq, "file_md5", Mock(side_effect=AssertionError("FASTQ read")))
    afterFastq.md5sum_worker(cfg)


def test_failed_hashing_never_publishes_partial_manifest(checksum_run, monkeypatch):
    cfg, manifest = checksum_run
    manifest.write_text("incomplete\n")
    monkeypatch.setattr(afterFastq, "file_md5", Mock(side_effect=OSError("read failed")))
    with pytest.raises(RuntimeError, match="FASTQ checksum generation failed"):
        afterFastq.md5sum_worker(cfg)
    assert manifest.read_text() == "incomplete\n"
    assert not list(cfg.output_path.glob(".*.tmp"))
    monkeypatch.undo()
    afterFastq.md5sum_worker(cfg)
    assert len(manifest.read_text().splitlines()) == 2


def test_new_demultiplexing_regenerates_checksums(checksum_run):
    cfg, manifest = checksum_run
    afterFastq.md5sum_worker(cfg)
    fastq = cfg.output_path / "GCF-2026-001/sample_R1.fastq.gz"
    fastq.write_bytes(b"new demultiplexed reads")
    afterFastq.md5sum_worker(cfg, force=True)
    assert hashlib.md5(fastq.read_bytes()).hexdigest() in manifest.read_text()
    assert manifest in cleanup_plan(cfg.output_path, "demultiplexing")


def test_escaped_filenames_and_binary_legacy_entries(checksum_run, monkeypatch):
    cfg, manifest = checksum_run
    (cfg.output_path / "GCF-2026-001/back\\slash\nnewline.fastq.gz").write_bytes(b"reads")
    afterFastq.md5sum_worker(cfg)
    manifest.write_text(manifest.read_text().replace("  GCF-", " *GCF-"))
    monkeypatch.setattr(afterFastq, "file_md5", Mock(side_effect=AssertionError("FASTQ read")))
    afterFastq.md5sum_worker(cfg)
