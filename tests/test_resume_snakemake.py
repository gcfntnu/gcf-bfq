"""Real tiny local DAG: no facility data, scientific tools, containers or mail."""

import subprocess
import sys

from pathlib import Path

import yaml


def test_partial_resume_noop_and_normal_parameter_invalidation(tmp_path, monkeypatch):
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))
    snakefile = tmp_path / "Snakefile"
    snakefile.write_text("""
configfile: "config.yaml"
rule multiqc_report:
    input: "report.txt"
rule expensive:
    output: "expensive.txt"
    shell: "printf expensive > {output}; printf 'expensive\\n' >> executed.log"
rule report:
    input: "expensive.txt"
    output: "report.txt"
    params: repair=config["repair"]
    shell: "printf '%s' '{params.repair}' > {output}; test '{params.repair}' != fail; printf 'report\\n' >> executed.log"
""")
    config = tmp_path / "config.yaml"
    config.write_text(yaml.safe_dump({"repair": "fail"}))
    executable = Path(sys.executable).with_name("snakemake")
    base = [
        str(executable),
        "--cores",
        "1",
        "--scheduler",
        "greedy",
        "--rerun-incomplete",
    ]

    def run(target, expected=0):
        result = subprocess.run(
            [*base, target], check=False, cwd=tmp_path, capture_output=True, text=True
        )
        assert result.returncode == expected, result.stdout + result.stderr
        with (tmp_path / "snakemake-transcripts.txt").open("a") as log:
            log.write(result.stdout + result.stderr)
        return result

    run("multiqc_report", expected=1)
    expensive = tmp_path / "expensive.txt"
    before = expensive.read_bytes(), expensive.stat().st_mtime_ns
    assert not (tmp_path / "report.txt").exists()  # Snakemake owns failed-job cleanup
    run("multiqc_report", expected=1)  # unchanged inputs must not report false success
    config.write_text(yaml.safe_dump({"repair": "repaired"}))
    run("multiqc_report")
    assert (expensive.read_bytes(), expensive.stat().st_mtime_ns) == before
    assert (tmp_path / "executed.log").read_text().splitlines() == ["expensive", "report"]
    report = tmp_path / "report.txt"
    completed = report.read_bytes(), report.stat().st_mtime_ns
    result = run("multiqc_report")
    assert "Nothing to be done" in result.stdout + result.stderr
    assert (report.read_bytes(), report.stat().st_mtime_ns) == completed
    config.write_text(yaml.safe_dump({"repair": "changed-parameter"}))
    run("multiqc_report")
    assert report.read_text() == "changed-parameter"
    assert (expensive.read_bytes(), expensive.stat().st_mtime_ns) == before
    assert (tmp_path / "executed.log").read_text().splitlines() == ["expensive", "report", "report"]
