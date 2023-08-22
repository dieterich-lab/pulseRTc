import pytest


def test_get_mismatches(get_bam):

    import pysam as ps

    bam = ps.AlignmentFile(get_bam, "rb")
    SN = (SQ["SN"] for SQ in bam.header["SQ"])
    bam.close()

    assert list(SN) == ["1", "2", "3"]

    # TODO: call utils.apply_parallel_iter
    # create reference dfs for mismatches
    # compare/assert df
    # create another test to test get_mismatch_details
    # compare dfs to ref dfs for counts (write to tmp loc)
    # also add test for filters (for now only mismatches)
    # and compare to remaining read names -> for that we need
    # tu put all filters in one function.
