import pytest
import pandas as pd

from pulsertc.splbam import unwrap_tuple_list
from pulsertc.splbam import get_mismatch_details
from pulsertc.splbam import filter_mismatches


def test_all_mismatches(get_mismatches):

    contig = [
        "1",
        "1",
        "1",
        "1",
        "1",
        "1",
        "2",
        "2",
        "2",
        "3",
        "3",
        "3",
        "3",
        "3",
        "3",
        "3",
    ]
    start = [
        14113,
        14201,
        14201,
        184629,
        184717,
        184717,
        101976764,
        101976764,
        101976852,
        156025613,
        156025618,
        156025633,
        156025665,
        156025618,
        156025633,
        156025665,
    ]
    end = [
        14114,
        14202,
        14202,
        184630,
        184718,
        184718,
        101976765,
        101976765,
        101976853,
        156025614,
        156025619,
        156025634,
        156025666,
        156025619,
        156025634,
        156025666,
    ]
    name = [
        "B",
        "B",
        "B",
        "B",
        "B",
        "B",
        "B",
        "B",
        "B",
        "C",
        "C",
        "C",
        "C",
        "C",
        "C",
        "C",
    ]
    score = [
        True,
        True,
        True,
        True,
        True,
        True,
        True,
        True,
        True,
        True,
        True,
        True,
        True,
        True,
        True,
        True,
    ]
    strand = [
        "+",
        "+",
        "+",
        "+",
        "+",
        "+",
        "-",
        "-",
        "-",
        "+",
        "+",
        "+",
        "+",
        "+",
        "+",
        "+",
    ]
    base = [
        "C",
        "C",
        "C",
        "C",
        "C",
        "C",
        "C",
        "C",
        "C",
        "G",
        "C",
        "T",
        "T",
        "A",
        "T",
        "T",
    ]
    ref = [
        "T",
        "T",
        "T",
        "T",
        "T",
        "T",
        "T",
        "T",
        "T",
        "T",
        "G",
        "C",
        "A",
        "G",
        "C",
        "A",
    ]
    qual = [37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 11, 37, 37, 11, 37, 37]
    pos = [9, 97, 62, 9, 97, 62, 38, 3, 91, 2, 7, 22, 54, 4, 19, 51]
    rlen = [101, 101, 101, 101, 101, 101, 101, 101, 101, 78, 78, 78, 78, 75, 75, 75]
    read1 = [
        False,
        False,
        True,
        False,
        False,
        True,
        True,
        False,
        False,
        True,
        True,
        True,
        True,
        False,
        False,
        False,
    ]
    samflag = [
        163,
        163,
        83,
        419,
        419,
        339,
        355,
        403,
        403,
        83,
        83,
        83,
        83,
        163,
        163,
        163,
    ]

    expected_mismatches = pd.DataFrame(
        zip(
            contig,
            start,
            end,
            name,
            score,
            strand,
            base,
            ref,
            qual,
            pos,
            rlen,
            read1,
            samflag,
        ),
        columns=[
            "contig",
            "start",
            "end",
            "name",
            "score",
            "strand",
            "base",
            "ref",
            "qual",
            "pos",
            "rlen",
            "read1",
            "samflag",
        ],
    )

    all_mismatches = unwrap_tuple_list(get_mismatches, idx=1, concat=True)

    pd.testing.assert_frame_equal(all_mismatches, expected_mismatches, check_exact=True)


def test_mismatch_count(get_mismatches, get_args):

    Category = ["Any", "Any", "Any", "Any", "Any", "Any", "Any", "Any", "Any"]
    Orientation = [
        "First",
        "First",
        "First",
        "First",
        "First",
        "Second",
        "Second",
        "Second",
        "Second",
    ]
    Genomic = ["A", "A", "C", "G", "T", "A", "C", "G", "T"]
    Read = ["C", "G", "G", "A", "A", "T", "T", "A", "C"]
    Coverage = [4, 4, 1, 2, 1, 4, 4, 1, 3]
    Mismatches = [1, 3, 1, 1, 1, 1, 1, 1, 6]

    expected_mismatch_count = pd.DataFrame(
        zip(Category, Orientation, Genomic, Read, Coverage, Mismatches),
        columns=[
            "Category",
            "Orientation",
            "Genomic",
            "Read",
            "Coverage",
            "Mismatches",
        ],
    )

    aligned_pairs = unwrap_tuple_list(get_mismatches)

    # for now, test overal mismatch counts
    # we should also test mismatch details (written to file)
    mismatch_count = get_mismatch_details(
        aligned_pairs, "/dev/null", "/dev/null", get_args.offset
    )

    pd.testing.assert_frame_equal(
        mismatch_count, expected_mismatch_count, check_exact=True
    )


def test_filters(get_mismatches, get_args):

    contig = ["1", "1", "1", "1", "2", "2"]
    start = [14113, 14201, 184629, 184717, 101976764, 101976852]
    end = [14114, 14202, 184630, 184718, 101976765, 101976853]
    name = ["B", "B", "B", "B", "B", "B"]
    score = [True, True, True, True, True, True]
    strand = ["+", "+", "+", "+", "-", "-"]
    base = ["C", "C", "C", "C", "C", "C"]
    ref = ["T", "T", "T", "T", "T", "T"]
    qual = [37, 37, 37, 37, 37, 37]
    pos = [9, 62, 9, 62, 38, 91]
    rlen = [101, 101, 101, 101, 101, 101]
    read1 = [False, True, False, True, True, False]
    samflag = [163, 83, 419, 339, 355, 403]

    expected_mismatches = pd.DataFrame(
        zip(
            contig,
            start,
            end,
            name,
            score,
            strand,
            base,
            ref,
            qual,
            pos,
            rlen,
            read1,
            samflag,
        ),
        columns=[
            "contig",
            "start",
            "end",
            "name",
            "score",
            "strand",
            "base",
            "ref",
            "qual",
            "pos",
            "rlen",
            "read1",
            "samflag",
        ],
    )

    aligned_pairs = unwrap_tuple_list(get_mismatches)
    mismatch_count = get_mismatch_details(
        aligned_pairs, "/dev/null", "/dev/null", get_args.offset
    )

    all_mismatches = unwrap_tuple_list(get_mismatches, idx=1, concat=True)

    # we should also test mismatch_count (written to file)
    all_mismatches = filter_mismatches(
        all_mismatches, mismatch_count, "/dev/null", get_args
    )

    pd.testing.assert_frame_equal(all_mismatches, expected_mismatches, check_exact=True)
