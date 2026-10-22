"""Workflow-owned 10x conversion profiles, independent of analysis selection."""

import hashlib
import re
import xml.etree.ElementTree as ET

from pathlib import Path

import yaml

from configmaker.indexes import IndexNormalizationError, normalize_rows
from configmaker.project_pairs import _override_geometry

AUTHORITATIVE_PROFILES = Path("/opt/gcf-workflows/demultiplexing.config")


class ProfileError(ValueError):
    """A conversion profile cannot safely represent the requested assay."""


def profile_for(custom):
    kit = custom.get("Libprep", "")
    if "10X Genomics" not in kit:
        backend = custom.get("Demultiplexer", "bcl-convert")
        if backend not in {"bcl-convert", "bcl2fastq"}:
            raise ProfileError(f"Unsupported Demultiplexer {backend!r} for {kit!r}")
        return {"backend": backend}
    try:
        content = AUTHORITATIVE_PROFILES.read_bytes()
        config = yaml.safe_load(content)
        if config["schema_version"] != 1:
            raise ValueError("unsupported schema")
        matches = [
            (name, value)
            for name, value in config["profiles"].items()
            if kit.removesuffix(" PE").removesuffix(" SE") in value["kits"]
        ]
        if len(matches) != 1:
            raise ValueError("no unique kit profile")
        name, profile = matches[0]
    except (OSError, KeyError, TypeError, ValueError, yaml.YAMLError) as error:
        raise ProfileError(
            f"Cannot select a workflow-owned demultiplexing profile for {kit!r}: {error}. "
            "Install the coordinated gcf-workflows #208 configuration; unlisted assays need an explicit reviewed profile (BFQ #156)."
        ) from error
    backend = custom.get("Demultiplexer", profile["default_backend"])
    if backend not in {"bcl-convert", "mkfastq"}:
        raise ProfileError("10x Demultiplexer must be bcl-convert or mkfastq")
    if backend == "bcl-convert" and profile["assay"] == "legacy":
        raise ProfileError(
            f"{name} retains legacy mkfastq: its read roles/barcode retention are unverified for BCL Convert. "
            "Select Demultiplexer,mkfastq; ATAC, joint Multiome and historical v1/v2 do not inherit a GEX profile."
        )
    return {
        **profile,
        "name": name,
        "backend": backend,
        "source": str(AUTHORITATIVE_PROFILES),
        "sha256": hashlib.sha256(content).hexdigest(),
        "bcl_convert_versions": config["bcl_convert_versions"],
    }


def configure(custom, settings, rows, run_info):
    """Return explicit sequences and verified effective BCL Convert settings."""
    profile = profile_for(custom)
    if profile["backend"] != "bcl-convert" or "name" not in profile:
        return profile, settings, rows
    orientation = custom.get("Index2Orientation")
    try:
        expanded = normalize_rows(
            [row for _number, row in rows], i5_orientation=orientation or "forward"
        )
    except IndexNormalizationError as error:
        raise ProfileError(str(error)) from error
    if (
        any(
            row["index2"] and row["_index_provenance"]["representation"] == "catalogue"
            for row in expanded
        )
        and not orientation
    ):
        raise ProfileError(
            "Named dual indexes require Index2Orientation,forward or reverse-complement in [CustomOptions]. "
            "Choose the documented instrument/reagent orientation; BFQ does not guess. Explicit sequence inputs are used as written."
        )
    if not run_info:
        raise ProfileError(
            "10x BCL Convert requires physical RunInfo.xml for cycle/read-role validation"
        )
    try:
        nodes = ET.parse(run_info).findall(".//Reads/Read")
        physical = [(int(node.attrib["NumCycles"]), node.attrib["IsIndexedRead"]) for node in nodes]
        if (
            len(physical) not in {3, 4}
            or [role for _cycles, role in physical] not in [["N", "Y", "N"], ["N", "Y", "Y", "N"]]
            or any(cycles <= 0 for cycles, _role in physical)
        ):
            raise ValueError("expected R1/i7/[i5]/R2 with positive cycles")
    except (OSError, ET.ParseError, ValueError, KeyError) as error:
        raise ProfileError(f"Invalid physical 10x read structure: {error}") from error
    effective = dict(settings)
    if effective.get("FastqCompressionFormat", "gzip").lower() != "gzip":
        raise ProfileError(
            "10x consumers require gzip FASTQs; remove FastqCompressionFormat or set it to gzip"
        )
    effective["FastqCompressionFormat"] = "gzip"
    if effective.get("NoLaneSplitting", "1").lower() not in {"1", "true"}:
        raise ProfileError(
            "NoLaneSplitting conflicts with BFQ's --no-lane-splitting true command; remove it or set it to true"
        )
    effective["NoLaneSplitting"] = "true"
    if "TrimUMI" in effective:
        raise ProfileError(
            "Remove TrimUMI from 10x settings: barcode/UMI cycles are retained as biological Y cycles; BCL Convert 4.2.4 rejects TrimUMI when no U cycles are defined"
        )
    for key, expected in {
        "CreateFastqForIndexReads": "0",
        "MinimumTrimmedReadLength": "8",
        "MaskShortReads": "8",
    }.items():
        if key in effective and effective[key] != expected:
            raise ProfileError(
                f"10x {profile['name']} requires {key},{expected}; correct the source settings explicitly"
            )
        effective[key] = expected
    lengths = set()
    output = []
    for row in expanded:
        for key in ("AdapterRead1", "AdapterRead2"):
            if (row.get(key) or settings.get(key) or "na").lower() not in {"", "na"}:
                raise ProfileError(
                    f"10x {key} trimming can damage barcodes/UMIs; remove it from the source sheet"
                )
            row[key] = "na"
        indexes = [row["index"], row["index2"]]
        if not indexes[0] or any(value and not re.fullmatch("[ACGT]+", value) for value in indexes):
            raise ProfileError("10x profiles require unambiguous explicit A/C/G/T sample indexes")
        if profile.get("require_dual_index") and not indexes[1]:
            raise ProfileError(f"{profile['name']} requires a complete dual-index assignment")
        if len(indexes[0]) not in profile["index_lengths"] or (
            indexes[1] and len(indexes[1]) != len(indexes[0])
        ):
            raise ProfileError(
                f"Unsupported index lengths for {profile['name']}: {[len(value) for value in indexes]}"
            )
        lengths.add(tuple(map(len, indexes)))
        tokens = [f"Y{physical[0][0]}"]
        for position, (cycles, _role) in enumerate(physical[1:-1]):
            used = len(indexes[position])
            if used > cycles:
                raise ProfileError("Sample index exceeds its physical index cycles")
            tokens.append(
                (f"I{used}" if used else "") + (f"N{cycles - used}" if cycles > used else "")
            )
        if len(physical) == 3 and indexes[1]:
            raise ProfileError("Dual-index assignment requires a physical i5 read")
        tokens.append(f"Y{physical[-1][0]}")
        mask = row.get("OverrideCycles") or settings.get("OverrideCycles") or ";".join(tokens)
        try:
            geometry, index_geometry = _override_geometry(
                mask, [(c, r == "Y") for c, r in physical], include_indexes=True
            )
        except ValueError as error:
            raise ProfileError(str(error)) from error
        if (
            "U" in mask.upper()
            or len(geometry) != 2
            or any(
                actual < minimum
                for actual, minimum in zip(
                    geometry, profile["minimum_biological_cycles"], strict=True
                )
            )
        ):
            raise ProfileError(
                f"{profile['name']} must retain barcode/UMI biological reads of at least {profile['minimum_biological_cycles']} cycles"
            )
        if index_geometry != [len(indexes[i]) for i in range(len(index_geometry))]:
            raise ProfileError(
                "OverrideCycles index allocation disagrees with the resolved index sequences"
            )
        row["OverrideCycles"] = mask
        for key in ("BarcodeMismatchesIndex1", "BarcodeMismatchesIndex2"):
            value = row.get(key) or settings.get(key) or "1"
            if value not in {"0", "1", "2"}:
                raise ProfileError(f"{key} must be 0, 1 or 2")
            row[key] = value
        provenance = row.pop("_index_provenance")
        output.append((rows[0][0], row, provenance))
    # Multiple per-library oligos are supported; heterogeneous index geometries
    # within one project need their own explicitly reviewed profile.
    if len(lengths) != 1:
        raise ProfileError("One 10x project must use one index geometry")
    # normalize_rows preserves source order; expand the corresponding source row.
    expanded_rows = []
    position = 0
    for number, original in rows:
        count = len(normalize_rows([original], i5_orientation=orientation or "forward"))
        for _ in range(count):
            _old, row, provenance = output[position]
            row["_index_provenance"] = provenance
            expanded_rows.append((number, row))
            position += 1
    masks = {row["OverrideCycles"] for _number, row in expanded_rows}
    if len(masks) == 1:
        effective["OverrideCycles"] = masks.pop()
    return profile, effective, expanded_rows
