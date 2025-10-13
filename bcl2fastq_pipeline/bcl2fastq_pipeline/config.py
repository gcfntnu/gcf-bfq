from __future__ import annotations

from configparser import ConfigParser
from dataclasses import dataclass, field
from pathlib import Path
from typing import ClassVar


@dataclass
class Paths:
    base_dir: Path
    output_dir: Path
    log_dir: Path | None = None
    report_dir: Path | None = None
    temp_dir: Path | None = None


@dataclass
class Options:
    run_id: str
    lanes: str | None = None
    sample_sheet: str | None = None
    libprep: str | None = None
    rerun: bool = False


@dataclass
class PipelineConfig:
    """Central configuration for a single GCF-BFQ run."""

    paths: Paths
    options: Options
    custom_options: dict[str, str] = field(default_factory=dict)

    # ---- Singleton boilerplate ----
    _instance: ClassVar[PipelineConfig | None] = None

    @classmethod
    def load(cls, ini_path: str | Path) -> PipelineConfig:
        """Initialize from a standard bcl2fastq.ini file."""
        parser = ConfigParser()
        parser.read(ini_path)

        paths = Paths(
            base_dir=Path(parser["Paths"]["baseDir"]),
            output_dir=Path(parser["Paths"]["outputDir"]),
            log_dir=Path(parser["Paths"].get("logDir", parser["Paths"]["outputDir"])) / "logs",
            report_dir=Path(parser["Paths"].get("reportDir", parser["Paths"]["outputDir"]))
            / "reports",
        )
        options = Options(
            run_id=parser["Options"]["runID"],
            lanes=parser["Options"].get("lanes"),
            sample_sheet=parser["Options"].get("sampleSheet"),
            libprep=parser["Options"].get("Libprep"),
            rerun=parser["Options"].getboolean("Rerun", fallback=False),
        )
        instance = cls(paths=paths, options=options)
        cls._instance = instance
        return instance

    @classmethod
    def get(cls) -> PipelineConfig:
        """Return the active configuration (singleton accessor)."""
        if cls._instance is None:
            raise RuntimeError("PipelineConfig not initialized. Call PipelineConfig.load() first.")
        return cls._instance

    def update_from_custom(self, custom_dict: dict[str, str]) -> None:
        """Update values from [CustomOptions] in sample sheet."""
        for key, value in custom_dict.items():
            normalized = key.lower()
            self.custom_options[normalized] = value
            # propagate known attributes into options if relevant
            if hasattr(self.options, normalized):
                setattr(self.options, normalized, value)

    # ---- Path helpers ----
    @property
    def run_path(self) -> Path:
        return self.paths.base_dir / self.options.run_id

    @property
    def output_path(self) -> Path:
        return self.paths.output_dir / self.options.run_id

    @property
    def log_path(self) -> Path:
        return self.paths.log_dir / f"{self.options.run_id}.log"

    @property
    def report_path(self) -> Path:
        return self.paths.report_dir / f"{self.options.run_id}_report.html"
