"""
config.py
=========
Central configuration management for the GCF-BFQ pipeline.

This module replaces the legacy global ConfigParser–based pattern with
a structured, typed configuration model composed of:

    • StaticConfig  – Immutable settings loaded once from bcl2fastq.ini
    • RunContext    – Mutable per-flowcell state
    • PipelineConfig – Singleton container tying them together

It also provides `parse_custom_options()` to read [CustomOptions]
from SampleSheet.csv files produced by the Django sample sheet cleaner.

The design intentionally preserves the “shared mutable state” behavior
of the legacy code: one global configuration object that represents the
current run, but with much cleaner access semantics and type safety.

Example
-------
>>> from bcl2fastq_pipeline.config import PipelineConfig, parse_custom_options
>>> cfg = PipelineConfig.load("/config/bcl2fastq.ini")
>>> cfg.run.begin("/mnt/seq/nova/240415_A01295_0345_BHXXXXXX")
>>> custom_opts, sheet_path = parse_custom_options(cfg.run.flowcell_path / "SampleSheet.csv")
>>> cfg.run.apply_custom(custom_opts, sheet_path)
>>> print(cfg.run.libprep)
'10X Genomics Chromium ...'
>>> print(cfg.static.paths.output_dir)
PosixPath('/mnt/data/output')
"""

from __future__ import annotations

import csv
import re

from configparser import ConfigParser
from dataclasses import dataclass, field
from pathlib import Path
from typing import ClassVar

import yaml


# --------------------------------------------------------------------------- #
# Data classes
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Paths:
    """
    Represents the [Paths] section of bcl2fastq.ini.

    All entries are converted to snake_case attributes automatically.

    Example attributes (depending on .ini):
        ekista_base_dir
        nova_base_dir
        output_dir
        log_dir
        report_dir
        analysis_dir
    """

    ekista_base_dir: Path
    nova_base_dir: Path
    output_dir: Path
    log_dir: Path
    manager_dir: Path
    report_dir: Path | None = None
    analysis_dir: Path | None = None


@dataclass(frozen=True)
class StaticConfig:
    """
    Immutable configuration loaded from bcl2fastq.ini.

    Contains:
        • paths        – Paths()
        • system       – dict of [System] section
        • email        – dict of [Email] section
        • version      – dict of [Version] section
        • commands     – dict of tool command strings (bcl2fastq, cellranger, etc.)
    """

    paths: Paths
    system: dict[str, str] = field(default_factory=dict)
    email: dict[str, str] = field(default_factory=dict)
    version: dict[str, str] = field(default_factory=dict)
    commands: dict[str, str] = field(default_factory=dict)


@dataclass
class RunContext:
    """
    Mutable configuration that describes the currently processed flowcell.

    Attributes
    ----------
    run_id : str
        Flowcell identifier (usually folder name).
    flowcell_path : Path
        Full path to the flowcell directory.
    base_dir : Path
        Parent of flowcell_path.
    instrument_source : str
        Either 'nova' or 'ekista' depending on where the flowcell was found.
    sample_sheet : Optional[Path]
        Path to the SampleSheet.csv used for this run.
    libprep : Optional[str]
        Library prep type from [CustomOptions].
    rerun : bool
        Flag for rerun status.
    custom : dict[str, str]
        Arbitrary key-value pairs from [CustomOptions].
    """

    run_id: str = ""
    flowcell_path: Path | None = None
    base_dir: Path | None = None
    instrument_source: str | None = None
    sample_sheet: Path | None = None
    sample_submission_form: Path | None = None
    libprep: str | None = None
    pipeline: str | None = None
    user: str | None = None
    rerun: bool = False
    sensitive: bool = False
    custom: dict[str, str] = field(default_factory=dict)

    # --- mutators ----------------------------------------------------------- #
    def begin(self, flowcell_path: Path, static_paths: Paths) -> None:
        """
        Initialize run context for a new flowcell.

        Automatically sets:
            • run_id = flowcell_path.name
            • base_dir = flowcell_path.parent
            • instrument_source = 'nova' or 'ekista' if path under known dirs

        Parameters
        ----------
        flowcell_path : Path
            Path to the flowcell directory.
        static_paths : Paths
            From StaticConfig, used to identify instrument source.
        """
        self.flowcell_path = Path(flowcell_path)
        self.run_id = self.flowcell_path.name
        self.base_dir = self.flowcell_path.parent

        # Detect which base_dir the flowcell belongs to
        try:
            f = self.flowcell_path.resolve()
            if str(f).startswith(str(static_paths.nova_base_dir.resolve())):
                self.instrument_source = "nova"
            elif str(f).startswith(str(static_paths.ekista_base_dir.resolve())):
                self.instrument_source = "ekista"
            else:
                self.instrument_source = "unknown"
        except Exception:
            self.instrument_source = "unknown"

    def apply_custom(
        self, custom_opts: dict[str, str], sheet_path: Path, sample_sub_path: Path
    ) -> None:
        """
        Apply [CustomOptions] from a SampleSheet.

        Stores all options in .custom and promotes relevant ones.

        Parameters
        ----------
        custom_opts : dict[str, str]
            Parsed [CustomOptions] entries.
        sheet_path : Path
            Path to the SampleSheet.csv used.
        sample_sub_path : Path
            Path to the SampleSubmissionForm.xlsx used.
        """
        self.sample_sheet = Path(sheet_path)
        self.sample_submission_form = Path(sample_sub_path)
        self.custom = {k.strip(): v.strip() for k, v in custom_opts.items() if k.strip()}

        # Promote well-known keys
        if "Libprep" in self.custom:
            self.libprep = self.custom["Libprep"]
        if "User" in self.custom:
            self.user = self.custom["User"]
        if "Rerun" in self.custom:
            val = self.custom["Rerun"].strip().lower()
            self.rerun = val in ("true", "1", "yes")
        if "SensitiveData" in self.custom:
            val = self.custom["SensitiveData"].strip().lower()
            self.sensitive = val in ("true", "1", "yes")

    def set_pipeline_from_yaml(self, yaml_path: Path) -> None:
        """
        Determine and set the pipeline workflow for this run based on libprep.config.

        The YAML structure is expected to look like:

            Lexogen SENSE mRNA-Seq Library Prep Kit V2 SE:
              workflow: rnaseq
              reads: SE
              adapter: AAAAAAAA
            Lexogen SENSE mRNA-Seq Library Prep Kit V2 PE:
              workflow: rnaseq
              reads: PE
            QIAseq miRNA SE:
              workflow: mirna

        This method:
          • Matches self.libprep (case-insensitive) against YAML keys
          • If not found, also tries "<libprep> SE" and "<libprep> PE"
          • Extracts the inner "workflow" value if present
          • Sets self.pipeline to that workflow name or "UNKNOWN" if not found
        """
        if not self.libprep:
            self.pipeline = None
            return

        yaml_file = Path(yaml_path)
        if not yaml_file.exists():
            print(f"[RunContext] libprep.config not found: {yaml_file}")
            self.pipeline = "UNKNOWN"
            return

        try:
            with open(yaml_file, encoding="utf-8") as fh:
                data = yaml.safe_load(fh) or {}
        except Exception as e:
            print(f"[RunContext] Failed to parse {yaml_file}: {e}")
            self.pipeline = "UNKNOWN"
            return

        libprep_clean = self.libprep.strip().lower()
        candidates = [libprep_clean, f"{libprep_clean} se", f"{libprep_clean} pe"]

        workflow = None

        for key, value in (data or {}).items():
            if not isinstance(value, dict):
                continue  # skip malformed entries

            key_norm = str(key).strip().lower()
            if key_norm in candidates:
                wf = value.get("workflow")
                if wf:
                    workflow = wf
                    break

        self.pipeline = str(workflow).strip() if workflow else "UNKNOWN"

    def reset(self) -> None:
        """Clear run-specific information."""
        self.run_id = ""
        self.flowcell_path = None
        self.base_dir = None
        self.instrument_source = None
        self.sample_sheet = None
        self.sample_submission_form = None
        self.user = None
        self.libprep = None
        self.rerun = False
        self.sensitive = False
        self.pipeline = None
        self.custom.clear()


# --------------------------------------------------------------------------- #
# Singleton container
# --------------------------------------------------------------------------- #
@dataclass
class PipelineConfig:
    """
    Singleton configuration object combining static and run-time settings.

    Use PipelineConfig.load(<path>) once per process to initialize,
    then access anywhere via PipelineConfig.get().

    Example
    -------
    >>> cfg = PipelineConfig.load("/config/bcl2fastq.ini")
    >>> cfg.run.begin("/mnt/seq/nova/240415_A01295_0345_BHXXXXXX")
    >>> cfg = PipelineConfig.get()
    >>> print(cfg.static.paths.output_dir)
    PosixPath('/mnt/data/output')
    """

    static: StaticConfig
    run: RunContext = field(default_factory=RunContext)
    _instance: ClassVar[PipelineConfig | None] = None

    # --- classmethods ------------------------------------------------------- #
    @classmethod
    def load(cls, ini_path: str | Path) -> PipelineConfig:
        """
        Load configuration from bcl2fastq.ini.

        Parameters
        ----------
        ini_path : str or Path
            Path to the bcl2fastq.ini file.

        Returns
        -------
        PipelineConfig
            The singleton instance.
        """
        parser = ConfigParser()
        parser.read(ini_path)

        # Normalize keys to lowercase for convenience
        p = parser["Paths"]

        paths = Paths(
            ekista_base_dir=Path(p.get("ekista_baseDir", "/mnt/seq/ekista")),
            nova_base_dir=Path(p.get("nova_baseDir", "/mnt/seq/nova")),
            output_dir=Path(p.get("outputDir", "/mnt/output")),
            log_dir=Path(p.get("logDir", "/mnt/logs")),
            report_dir=Path(p.get("reportDir", "/mnt/reports")),
            analysis_dir=Path(p.get("analysisDir", "/mnt/analysis")),
            manager_dir=Path(p.get("manager_dir", "/mnt/manager")),
        )

        static = StaticConfig(
            paths=paths,
            system=dict(parser.items("System")) if parser.has_section("System") else {},
            email=dict(parser.items("Email")) if parser.has_section("Email") else {},
            version=dict(parser.items("Version")) if parser.has_section("Version") else {},
            commands={
                k: v
                for k, v in parser.items("Commands")  # example section
            }
            if parser.has_section("Commands")
            else {},
        )

        instance = cls(static=static)
        cls._instance = instance
        return instance

    @classmethod
    def get(cls) -> PipelineConfig:
        """Return the active configuration instance."""
        if cls._instance is None:
            raise RuntimeError("PipelineConfig not initialized. Call PipelineConfig.load() first.")
        return cls._instance

    @property
    def output_path(self) -> Path | None:
        if not self.run.run_id:
            return None
        return self.static.paths.output_dir / self.run.run_id

    def to_file(self, out_path: Path) -> None:
        """
        Write all relevant configuration information to a YAML file.

        The file will include both static (environment) and run-specific
        sections in a human-readable structure:

            static:
              paths:
                ekista_base_dir: /mnt/seq/ekista
                ...
              version:
                gcf_bfq: "1.4.0"
            run:
              run_id: 240415_A01295_0345_BHXXXXXX
              instrument_source: nova
              libprep: Lexogen SENSE mRNA-Seq Library Prep Kit V2
              pipeline: rnaseq
              user: vidar
              custom:
                TrimAdapter: "True"
                ...
        Parameters
        ----------
        out_path : Path
            Path to the output YAML file (created or overwritten).
        """
        # --- convert StaticConfig to simple serializable dict ---
        static_dict = {
            "paths": {k: str(v) for k, v in vars(self.static.paths).items() if v is not None},
            "system": self.static.system,
            "email": self.static.email,
            "version": self.static.version,
            "commands": self.static.commands,
        }

        # --- convert RunContext to dict ---
        run_fields = {
            "run_id": self.run.run_id,
            "flowcell_path": str(self.run.flowcell_path) if self.run.flowcell_path else None,
            "base_dir": str(self.run.base_dir) if self.run.base_dir else None,
            "instrument_source": self.run.instrument_source,
            "sample_sheet": str(self.run.sample_sheet) if self.run.sample_sheet else None,
            "sample_submission_form": str(self.run.sample_submission_form)
            if self.run.sample_submission_form
            else None,
            "libprep": self.run.libprep,
            "pipeline": self.run.pipeline,
            "user": self.run.user,
            "rerun": self.run.rerun,
            "custom": self.run.custom or {},
        }

        cfg_dict = {"static": static_dict, "run": run_fields}

        out_path = Path(out_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)

        with open(out_path, "w", encoding="utf-8") as fh:
            yaml.safe_dump(cfg_dict, fh, sort_keys=False, default_flow_style=False)


# --------------------------------------------------------------------------- #
# SampleSheet parser
# --------------------------------------------------------------------------- #
def parse_custom_options(sample_sheet_path: Path) -> tuple[dict[str, str], Path]:
    """
    Parse the [CustomOptions] section from a SampleSheet.csv.

    Handles messy CSV formatting from Excel edits, trailing commas, etc.

    Parsing strategy:
        - Scan lines until '[CustomOptions]' (case-insensitive)
        - Read subsequent lines until EOF or another '[Section]'
        - Split each line by comma and take only first two fields
        - Ignore empty or malformed rows

    Parameters
    ----------
    sample_sheet_path : Path
        Path to the SampleSheet.csv.

    Returns
    -------
    (dict, Path)
        Dictionary of key-value pairs and the same Path for reference.
    """
    custom_opts: dict[str, str] = {}
    in_custom_section = False
    section_pattern = re.compile(r"^\s*\[.*\]\s*$")

    with open(sample_sheet_path, newline="", encoding="utf-8-sig") as fh:
        reader = csv.reader(fh)
        for row in reader:
            if not row:
                continue
            line = ",".join(row).strip()
            if not in_custom_section:
                if re.match(r"^\s*\[CustomOptions\]\s*$", line, flags=re.IGNORECASE):
                    in_custom_section = True
                continue

            # Exit when encountering another section header
            if section_pattern.match(line) and not line.lower().startswith("[customoptions]"):
                break

            # Parse first two columns if present
            key = (row[0] or "").strip()
            value = (row[1] if len(row) > 1 else "").strip()
            if key:
                custom_opts[key] = value

    return custom_opts, sample_sheet_path
