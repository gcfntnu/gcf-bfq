import pytest

from bcl2fastq_pipeline.index_sequences import IndexSequenceError, orientation, transform
from bcl2fastq_pipeline.state import StateError


def sheet(rows="s1,AAAC,CGAA\ns2,GGGA,TGTT\n", header="Sample_ID,index,index2"):
    return f"[Data]\n{header}\n{rows}".encode()


def test_edit_preserves_every_unselected_byte_and_repeated_toggle_restores_input():
    original = (
        b'\xef\xbb\xbf[Header],,,\r\nDescription,"a,quoted \'note\' and ""escape""",,\n'
        b"[Settings],,,\rReverseComplementIndexP5,1\r\nReverseComplementIndexP7,0\r\n"
        b"[Data],,,,\r\nLane,Sample_ID,index,index2,Description\r\n"
        b'1,"s1","AaCG",CCCA,"multiline\r\ntext, with ""quotes"""\r\n'
        b'2,s2,TTAA,"gAtC",untouched\n'
        b"[Other]\nindex,GGGT\n"
    )
    result = transform(original, ("index1", "index2"))
    expected = original.replace(b'"AaCG",CCCA', b'"CGtT",TGGG').replace(b'"gAtC"', b'"GaTc"')
    assert result.content == expected
    assert result.counts == {"index1": 1, "index2": 2}
    assert result.selected_counts == {"index1": 2, "index2": 2}
    assert transform(result.content, ("index1", "index2")).content == original


@pytest.mark.parametrize(
    "indexes,expected",
    [
        (("index1",), b"s1,GTTT,CGAA\ns2,TCCC,TGTT\n"),
        (("index2",), b"s1,AAAC,TTCG\ns2,GGGA,AACA\n"),
        (("index1", "index2"), b"s1,GTTT,TTCG\ns2,TCCC,AACA\n"),
        (("index2", "index2"), b"s1,AAAC,TTCG\ns2,GGGA,AACA\n"),
    ],
)
def test_independent_combined_and_duplicate_selection(indexes, expected):
    assert transform(sheet(), indexes).content == b"[Data]\nSample_ID,index,index2\n" + expected


def test_sequential_toggles_equal_combined_and_then_restore():
    original = sheet()
    once = transform(original, ("index1",)).content
    twice = transform(once, ("index2",)).content
    assert twice == transform(original, ("index1", "index2")).content
    assert transform(transform(twice, ("index1",)).content, ("index2",)).content == original


def test_iupac_and_case_are_supported():
    original = sheet("s1,ACGTRYSWKMBDHVNacgtryswkmbdhvn,AAAA\n")
    edited = transform(original, ("index1",)).content
    assert b"nbdhvkmwsryacgtNBDHVKMWSRYACGT" in edited
    assert transform(edited, ("index1",)).content == original


@pytest.mark.parametrize("newline", [b"\n", b"\r\n", b"\r"])
@pytest.mark.parametrize("terminal", [True, False])
def test_line_endings_and_final_newline_preserved(newline, terminal):
    original = newline.join([b"[Data]", b"Sample_ID,index,index2", b"s1,AAAC,CGAA"])
    if terminal:
        original += newline
    assert transform(original, ("index1",)).content == original.replace(b"AAAC", b"GTTT")


def test_empty_selected_cells_remain_unchanged():
    original = sheet('s1,,AAAA\ns2,"",CCCC\ns3,AAAC,GGGG\n')
    result = transform(original, ("index1",))
    assert result.content == original.replace(b"AAAC", b"GTTT")
    assert result.counts == result.selected_counts == {"index1": 1}


def test_palindromic_values_are_valid_even_when_no_bytes_change():
    original = sheet("s1,ACGT,TTAA\n")
    result = transform(original, ("index1",))
    assert result.content == original
    assert result.counts == {"index1": 0}
    assert result.selected_counts == {"index1": 1}


def test_bcl_convert_section_and_trailing_empty_columns():
    original = b"[BCLConvert_Data],,\nSample_ID,index,index2,,,\ns1,AAAC,CGAA\ns2,GGGA,TGTT,,,\n"
    assert b"s1,GTTT,CGAA\n" in transform(original, ("index1",)).content


def test_unselected_invalid_sequences_are_not_changed_or_validated():
    original = sheet("s1,AAAC,an-unsupported-value\n")
    assert transform(original, ("index1",)).content == original.replace(b"AAAC", b"GTTT")


@pytest.mark.parametrize(
    "original,indexes,match",
    [
        (sheet(), (), "Select index1"),
        (sheet(), ("index",), "Select index1"),
        (sheet("s1,AAAC,CG-AA\n"), ("index1", "index2"), "Invalid 'index2'"),
        (sheet("s1,AAAC,CGAA \n"), ("index2",), "Invalid 'index2'"),
        (sheet("s1,AAAC,UCGT\n"), ("index2",), "Invalid 'index2'"),
        (sheet("s1,AAAC,\n"), ("index2",), "no index sequences"),
        (sheet("s1,AAAC\n"), ("index2",), "missing its 'index2'"),
        (sheet("s1,AAAC\n", "Sample_ID,index"), ("index2",), "missing selected column"),
        (sheet("s1,AAAC,CGAA,extra\n"), ("index1",), "more values than"),
        (sheet("", "Sample_ID,index,index2"), ("index1",), "no sample rows"),
        (b"[Data]\n", ("index1",), "no header"),
        (b"Sample_ID,index\ns1,AAAC\n", ("index1",), "exactly one"),
        (sheet() + sheet(), ("index1",), "exactly one"),
        (sheet() + b"[BCLConvert_Data]\nSample_ID,index\ns1,AAAC\n", ("index1",), "exactly one"),
        (sheet("s1,AAAC,CGAA\n", "Sample_ID,index,index"), ("index1",), "duplicate column"),
        (b"[Data]\nSample_ID;index;index2\ns1;AAAC;CGAA\n", ("index1",), "comma-delimited"),
        (b"[Data]\nSample_ID\tindex\tindex2\ns1\tAAAC\tCGAA\n", ("index1",), "comma-delimited"),
        (sheet('s1,"AAAC,CGAA\n'), ("index1",), "unterminated"),
        (sheet('s1,A"AAC,CGAA\n'), ("index1",), "unquoted field"),
        (sheet('s1,"AAAC"extra,CGAA\n'), ("index1",), "after a closing quote"),
        (sheet('s1,"AA""AC",CGAA\n'), ("index1",), "Invalid 'index'"),
        (sheet('s1,"AA\nAC",CGAA\n'), ("index1",), "Invalid 'index'"),
        (sheet() + b"\xff", ("index1",), "UTF-8"),
        (sheet() + b"\x00", ("index1",), "NUL"),
    ],
)
def test_invalid_edits_raise_state_error_without_mutating_input(original, indexes, match):
    before = bytes(original)
    with pytest.raises(IndexSequenceError, match=match) as error:
        transform(original, indexes)
    assert isinstance(error.value, StateError)
    assert original == before


def test_orientation_identifies_original_and_independent_toggles():
    original = sheet()
    assert {index: info["status"] for index, info in orientation(original, original).items()} == {
        "index1": "original",
        "index2": "original",
    }
    result = orientation(transform(original, ("index2",)).content, original)
    assert result["index1"]["status"] == "original"
    assert result["index2"]["status"] == "reversed"
    assert (
        orientation(
            transform(transform(original, ("index2",)).content, ("index2",)).content, original
        )["index2"]["status"]
        == "original"
    )


def test_orientation_matches_sample_lane_project_without_depending_on_row_order():
    header = "Lane,Sample_ID,Sample_Project,index,index2"
    original = sheet("1,s1,P1,AAAC,CGAA\n2,s1,P1,GGGA,TGTT\n1,s1,P2,AAAA,CCCC\n", header)
    reordered = sheet("1,s1,P2,AAAA,CCCC\n2,s1,P1,GGGA,TGTT\n1,s1,P1,AAAC,CGAA\n", header)
    result = orientation(transform(reordered, ("index1",)).content, original)
    assert result["index1"]["status"] == "reversed"
    assert result["index2"]["status"] == "original"


@pytest.mark.parametrize(
    "current,status",
    [
        (sheet(), "original"),
        (sheet("s1,GTTT,CGAA\ns2,TCCC,TGTT\n"), "reversed"),
        (sheet("s1,GTTT,CGAA\ns2,GGGA,TGTT\n"), "mixed"),
        (sheet("s1,TTTT,CGAA\ns2,GGGA,TGTT\n"), "modified"),
        (sheet("s1,AAAC,CGAA\ns3,GGGA,TGTT\n"), "modified"),
        (sheet("s1,AAAC,CGAA\ns1,GGGA,TGTT\n"), "ambiguous"),
        (sheet("s1,aaac,CGAA\ns2,ggga,TGTT\n"), "original"),
        (sheet("s1,AAAC,CGAA\ns2,GG-A,TGTT\n"), "unavailable"),
        (sheet("s1,,CGAA\ns2,GGGA,TGTT\n"), "modified"),
    ],
)
def test_orientation_statuses(current, status):
    assert orientation(current, sheet())["index1"]["status"] == status


def test_palindromes_and_empty_values_do_not_falsely_establish_orientation():
    palindromes = sheet("s1,ACGT,TTAA\ns2,,CCGG\n")
    assert orientation(palindromes, palindromes)["index1"]["status"] == "ambiguous"
    informative = palindromes + b"s3,AAAC,CGAA\n"
    assert orientation(informative, informative)["index1"]["status"] == "original"
    assert (
        orientation(transform(informative, ("index1",)).content, informative)["index1"]["status"]
        == "reversed"
    )


def test_unavailable_and_absent_orientation_data():
    single = sheet("s1,AAAC\n", "Sample_ID,index")
    assert orientation(single, single)["index2"]["status"] == "not_present"
    assert orientation(sheet(), single)["index2"]["status"] == "unavailable"
    assert orientation(sheet(), b"")["index1"]["status"] == "unavailable"
    missing_identity = sheet("AAAC,CGAA\n", "index,index2")
    assert orientation(missing_identity, missing_identity)["index1"]["status"] == "unavailable"


def test_orientation_does_not_silently_ignore_lane_changes():
    original = sheet("1,s1,AAAC,CGAA\n", "Lane,Sample_ID,index,index2")
    current = sheet("2,s1,AAAC,CGAA\n", "Lane,Sample_ID,index,index2")
    assert orientation(current, original)["index1"]["status"] == "modified"


def test_added_project_column_does_not_prevent_matching_source_samples():
    original = sheet("s1,AAAC,CGAA\n")
    current = sheet("s1,P1,AAAC,CGAA\n", "Sample_ID,Sample_Project,index,index2")
    assert orientation(current, original)["index1"]["status"] == "original"
