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
    libprep: str | None = None
    rerun: bool = False
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

    def apply_custom(self, custom_opts: dict[str, str], sheet_path: Path) -> None:
        """
        Apply [CustomOptions] from a SampleSheet.

        Stores all options in .custom and promotes relevant ones.

        Parameters
        ----------
        custom_opts : dict[str, str]
            Parsed [CustomOptions] entries.
        sheet_path : Path
            Path to the SampleSheet.csv used.
        """
        self.sample_sheet = Path(sheet_path)
        self.custom = {k.strip(): v.strip() for k, v in custom_opts.items() if k.strip()}

        # Promote well-known keys
        if "Libprep" in self.custom:
            self.libprep = self.custom["Libprep"]
        if "Rerun" in self.custom:
            val = self.custom["Rerun"].strip().lower()
            self.rerun = val in ("true", "1", "yes")

    def reset(self) -> None:
        """Clear run-specific information."""
        self.run_id = ""
        self.flowcell_path = None
        self.base_dir = None
        self.instrument_source = None
        self.sample_sheet = None
        self.libprep = None
        self.rerun = False
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
