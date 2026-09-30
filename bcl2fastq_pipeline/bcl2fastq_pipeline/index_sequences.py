"""Lossless, in-memory SampleSheet index edits and source-relative orientation.

Only comma-delimited UTF-8 CSV with one [Data] or [BCLConvert_Data] section is
supported. Index cells may use DNA IUPAC symbols in either case. Empty cells are
left alone, but a requested column with no nonempty values is an error. Parsing
records into byte spans, rather than serializing CSV again, preserves unrelated
content, delimiters, quoting, the BOM, and each existing line ending exactly.
"""

from __future__ import annotations

from dataclasses import dataclass

from bcl2fastq_pipeline.state import StateError

INDEX_COLUMNS = {"index1": "index", "index2": "index2"}
_DNA = b"ACGTRYSWKMBDHVNacgtryswkmbdhvn"
_COMPLEMENT = bytes.maketrans(_DNA, b"TGCAYRSWMKVHDBNtgcayrswmkvhdbn")
_DNA_SET = frozenset(_DNA)


class IndexSequenceError(StateError):
    """The SampleSheet cannot safely undergo the requested index edit."""


@dataclass(frozen=True)
class TransformResult:
    content: bytes
    counts: dict[str, int]
    selected_counts: dict[str, int]


@dataclass(frozen=True)
class _Cell:
    value: bytes
    start: int
    end: int


@dataclass(frozen=True)
class _Sheet:
    columns: dict[str, int]
    rows: list[list[_Cell]]


def _records(content: bytes) -> list[list[_Cell]]:
    """Parse strict CSV while retaining spans of field payloads in the input."""
    try:
        content.decode("utf-8-sig")
    except UnicodeDecodeError as error:
        raise IndexSequenceError(
            "SampleSheet must be UTF-8 CSV (an initial BOM is supported)"
        ) from error
    if b"\x00" in content:
        raise IndexSequenceError("SampleSheet contains a NUL byte; expected UTF-8 CSV")
    position = 3 if content.startswith(b"\xef\xbb\xbf") else 0
    size = len(content)
    records = []
    while position < size:
        row = []
        while True:
            if position < size and content[position] == ord('"'):
                position += 1
                start = position
                value = bytearray()
                while True:
                    if position >= size:
                        raise IndexSequenceError("Malformed CSV: unterminated quoted field")
                    character = content[position]
                    if character == ord('"'):
                        if position + 1 < size and content[position + 1] == ord('"'):
                            value.append(character)
                            position += 2
                            continue
                        end = position
                        position += 1
                        break
                    value.append(character)
                    position += 1
                if position < size and content[position] not in b",\r\n":
                    raise IndexSequenceError("Malformed CSV: unexpected text after a closing quote")
                row.append(_Cell(bytes(value), start, end))
            else:
                start = position
                while position < size and content[position] not in b",\r\n":
                    if content[position] == ord('"'):
                        raise IndexSequenceError("Malformed CSV: quote inside an unquoted field")
                    position += 1
                row.append(_Cell(content[start:position], start, position))
            if position == size:
                records.append(row)
                return records
            if content[position] == ord(","):
                position += 1
                continue
            if content[position : position + 2] == b"\r\n":
                position += 2
            else:
                position += 1
            records.append(row)
            break
    return records


def _blank(row: list[_Cell]) -> bool:
    return all(not cell.value.strip() for cell in row)


def _section(row: list[_Cell]) -> bytes | None:
    first = row[0].value.strip()
    if first.startswith(b"[") and first.endswith(b"]") and _blank(row[1:]):
        return first.lower()
    return None


def _read_sheet(content: bytes) -> _Sheet:
    records = _records(content)
    starts = [
        number
        for number, row in enumerate(records)
        if _section(row) in {b"[data]", b"[bclconvert_data]"}
    ]
    if len(starts) != 1:
        raise IndexSequenceError(
            "Expected exactly one [Data] or [BCLConvert_Data] section in comma-delimited CSV"
        )
    data = []
    for row in records[starts[0] + 1 :]:
        if _section(row) is not None:
            break
        if not _blank(row):
            data.append(row)
    if not data:
        raise IndexSequenceError("SampleSheet data section has no header")
    header = data[0]
    columns = {}
    for number, cell in enumerate(header):
        name = cell.value.decode("utf-8").strip().lower()
        if not name:
            continue
        if name in columns:
            raise IndexSequenceError(f"SampleSheet data section has duplicate column {name!r}")
        columns[name] = number
    if len(header) == 1 and any(delimiter in header[0].value for delimiter in (b";", b"\t")):
        raise IndexSequenceError("Only comma-delimited SampleSheet CSV is supported")
    if len(data) == 1:
        raise IndexSequenceError("SampleSheet data section has no sample rows")
    for number, row in enumerate(data[1:], start=1):
        if len(row) > len(header) and not _blank(row[len(header) :]):
            raise IndexSequenceError(
                f"SampleSheet data row {number} has more values than its header"
            )
    return _Sheet(columns, data[1:])


def _index_cells(sheet: _Sheet, index: str) -> list[_Cell]:
    column = INDEX_COLUMNS[index]
    if column not in sheet.columns:
        raise IndexSequenceError(
            f"SampleSheet data section is missing selected column {column!r} ({index})"
        )
    position = sheet.columns[column]
    cells = []
    for number, row in enumerate(sheet.rows, start=1):
        if position >= len(row):
            raise IndexSequenceError(
                f"SampleSheet data row {number} is missing its {column!r} cell"
            )
        cell = row[position]
        if any(character not in _DNA_SET for character in cell.value):
            raise IndexSequenceError(
                f"Invalid {column!r} sequence in data row {number}: expected DNA IUPAC symbols "
                "ACGTRYSWKMBDHVN, without spaces"
            )
        cells.append(cell)
    return cells


def reverse_complement(value: bytes) -> bytes:
    """Reverse-complement a previously validated DNA sequence, retaining case."""
    return value.translate(_COMPLEMENT)[::-1]


def transform(content: bytes, indexes: tuple[str, ...]) -> TransformResult:
    """Toggle selected columns once each; validate every edit before returning.

    ``counts`` counts cells whose bytes change, while ``selected_counts`` counts
    all nonempty cells processed (including palindromes). Repeating a successful
    transform toggles the same cells back to their exact prior byte content.
    """
    selected = tuple(dict.fromkeys(indexes))
    if not selected or any(index not in INDEX_COLUMNS for index in selected):
        raise IndexSequenceError("Select index1, index2, or both for reverse-complementing")
    sheet = _read_sheet(content)
    replacements = []
    counts = {}
    selected_counts = {}
    for index in selected:
        cells = _index_cells(sheet, index)
        selected_counts[index] = sum(bool(cell.value) for cell in cells)
        if not selected_counts[index]:
            raise IndexSequenceError(
                f"Selected column {INDEX_COLUMNS[index]!r} ({index}) has no index sequences"
            )
        counts[index] = 0
        for cell in cells:
            replacement = reverse_complement(cell.value)
            if replacement != cell.value:
                counts[index] += 1
                replacements.append((cell.start, cell.end, replacement))
    parts = []
    previous = 0
    for start, end, replacement in sorted(replacements):
        parts.extend((content[previous:start], replacement))
        previous = end
    parts.append(content[previous:])
    return TransformResult(b"".join(parts), counts, selected_counts)


def _row_keys(sheet: _Sheet, include_project: bool) -> list[tuple[bytes, ...]]:
    sample_column = next(
        (sheet.columns[name] for name in ("sample_id", "sampleid") if name in sheet.columns), None
    )
    if sample_column is None:
        raise IndexSequenceError("Sample_ID is required to compare orientation with the source")
    lane_column = sheet.columns.get("lane")
    project_column = next(
        (
            sheet.columns[name]
            for name in ("sample_project", "sampleproject")
            if name in sheet.columns
        ),
        None,
    )
    columns = [sample_column, lane_column]
    if include_project:
        columns.append(project_column)
    keys = []
    for row in sheet.rows:
        key = tuple(
            row[column].value.strip() if column is not None and column < len(row) else b""
            for column in columns
        )
        if not key[0]:
            raise IndexSequenceError("An empty Sample_ID prevents comparison with the source")
        keys.append(key)
    return keys


def _result(status: str, detail: str) -> dict[str, str]:
    return {"status": status, "detail": detail}


def _compare_index(current: _Sheet, source: _Sheet, index: str) -> dict[str, str]:
    column = INDEX_COLUMNS[index]
    if column not in current.columns:
        return _result("not_present", f"{column} is absent from the effective SampleSheet")
    if column not in source.columns:
        return _result("unavailable", f"{column} is absent from the source SampleSheet")
    try:
        current_cells = _index_cells(current, index)
        source_cells = _index_cells(source, index)
        include_project = any(
            name in current.columns for name in ("sample_project", "sampleproject")
        ) and any(name in source.columns for name in ("sample_project", "sampleproject"))
        current_keys = _row_keys(current, include_project)
        source_keys = _row_keys(source, include_project)
    except IndexSequenceError as error:
        return _result("unavailable", str(error))
    if len(set(current_keys)) != len(current_keys) or len(set(source_keys)) != len(source_keys):
        return _result(
            "ambiguous",
            "Duplicate sample/lane keys prevent an unambiguous comparison with the source",
        )
    if set(current_keys) != set(source_keys):
        return _result("modified", "Sample/lane entries differ from the source SampleSheet")
    current_values = dict(
        zip(current_keys, (cell.value.upper() for cell in current_cells), strict=True)
    )
    source_values = dict(
        zip(source_keys, (cell.value.upper() for cell in source_cells), strict=True)
    )
    return _classify_orientation(current_values, source_values)


def _classify_orientation(
    current_values: dict[tuple[bytes, ...], bytes],
    source_values: dict[tuple[bytes, ...], bytes],
) -> dict[str, str]:
    evidence = set()
    for key, original in source_values.items():
        value = current_values[key]
        reversed_value = reverse_complement(original)
        if value == original == reversed_value:
            continue
        if value == original:
            evidence.add("original")
        elif value == reversed_value:
            evidence.add("reversed")
        else:
            return _result(
                "modified", "Index values include edits other than reverse-complementing the source"
            )
    if evidence == {"original"}:
        return _result("original", "Matches source (original orientation)")
    if evidence == {"reversed"}:
        return _result("reversed", "Reversed relative to source")
    if evidence:
        return _result(
            "mixed", "Some index values match the source and others are reverse-complemented"
        )
    return _result(
        "ambiguous",
        "All values are empty or self reverse-complementary; orientation is indistinguishable",
    )


def orientation(content: bytes, source_content: bytes) -> dict[str, dict[str, str]]:
    """Describe each index relative to source using sample/lane identities.

    Comparison is case-insensitive, ignores row order, and includes project in
    the identity when both sheets provide it. Missing/unparseable comparison
    data produces a status, never an exception that would prevent a valid edit.
    """
    try:
        current = _read_sheet(content)
        source = _read_sheet(source_content)
    except IndexSequenceError as error:
        return {index: _result("unavailable", str(error)) for index in INDEX_COLUMNS}
    return {index: _compare_index(current, source, index) for index in INDEX_COLUMNS}
