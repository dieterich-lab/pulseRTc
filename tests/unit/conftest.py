"""
    confest.py
"""

import pytest


class mockParser:
    def parse_args(self):
        self.overwrite = True
        self.do_not_call = False
        self.ref_base = "T"
        self.base_change = "C"
        self.qual = 20
        self.trim5p = 5
        self.trim3p = 5
        self.subtract = False
        self.library_type = "reverse"
        self.offset = 100
        self.num_cpus = 1


class BAMData:

    """\
    Test BAM file.
    """

    def get_header(self):
        header = {
            "HD": {"VN": "1.5"},
            "SQ": [
                {"LN": 248956422, "SN": "1"},
                {"LN": 101991189, "SN": "2"},
                {"LN": 156040895, "SN": "3"},
            ],
        }
        self.header = header

    def get_entries(self):
        query_names = ["A", "A", "B", "B", "B", "B", "B", "B", "C", "C"]
        flags = [419, 339, 163, 83, 419, 339, 355, 403, 83, 163]
        reference_ids = ["1", "1", "1", "1", "1", "1", "2", "2", "3", "3"]
        reference_starts = [
            12006,
            12065,
            14104,
            14139,
            184620,
            184655,
            101976726,
            101976761,
            156025611,
            156025614,
        ]
        mapping_qualities = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        cigars = [
            "99M",
            "101M",
            "101M",
            "101M",
            "101M",
            "101M",
            "101M",
            "101M",
            "78M",
            "75M",
        ]
        next_reference_ids = ["=", "=", "=", "=", "=", "=", "=", "=", "=", "="]
        next_reference_starts = [
            12065,
            12006,
            14139,
            14104,
            184655,
            184620,
            101976761,
            101976726,
            156025614,
            156025611,
        ]
        template_lengths = [160, -160, 136, -136, 136, -136, 136, -136, -75, 75]
        query_sequences = [
            "CAGGTGTCTGACTTCCAGCAACTGCTGGCCTGTGCCAGGGTGCAAGCTGAGCACTGGAGTGGAGTTTTCCTGTGGAGAGGAGCCATGCCTAGAGTGGGA",
            "TGGAGTTTTCCTGTGGAGAGGAGCCATGCCTAGAGTGGGATGGGCCATTGTTCATCTTCTGGCCCCTGTTGTCTGCATGTAACTTAATACCACAACCAGGC",
            "GACATCTTCCACCCCAACACCAGCAATTGTGCCAAGGGCCATTAGGCTCTCAGCATGACTATTTTTAGAGACCCCGTGTCTGTCACTGAAACCTTTTCTGT",
            "GGGCCATTAGGCTCTCAGCATGACTATTTTTAGAGACCCCGTGTCTGTCACTGAAACCTTTTCTGTGGGAGACTATTCCTCCCATCTGCAACAGCTGCCCC",
            "GACATCTTCCACCCCAACACCAGCAATTGTGCCAAGGGCCATTAGGCTCTCAGCATGACTATTTTTAGAGACCCCGTGTCTGTCACTGAAACCTTTTCTGT",
            "GGGCCATTAGGCTCTCAGCATGACTATTTTTAGAGACCCCGTGTCTGTCACTGAAACCTTTTCTGTGGGAGACTATTCCTCCCATCTGCAACAGCTGCCCC",
            "GGGGCAGCTGTTGCAGATGGGAGGAATAGTCTCCCACAGAAAAGGTTTCAGTGACAGACACGGGGTCTCTAAAAATAGTCATGCTGAGAGCCTAATGGCCC",
            "ACAGAAAAGGTTTCAGTGACAGACACGGGGTCTCTAAAAATAGTCATGCTGAGAGCCTAATGGCCCTTGGCACAATTGCTGGTGTTGGGGTGGAAGATGTC",
            "CTGAAGACGCTTCCAGAGAAAATGGCACACCAATCAATAAAGAACTGAGCAGAATCCAACAGTGTGCTTTTAATAAAG",
            "AAGAAGCTTCCAGAGAAAATGGCACACCAATCAATAAAGAACTGAGCAGAATCCAACAGTGTGCTTTTAATAAAG",
        ]
        query_qualities = [
            "FFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFF,FFFFFFFF",
            "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF",
            "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFF",
            "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFF",
            "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFF",
            "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFF",
            "FFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF",
            "FFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF",
            "FFFFFFF,FFFFFFFFFFFF:FFFFFFFF,:FFF,FFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFF::F",
            "FFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF",
        ]
        tags = [
            (("MD:Z", "99"), ("HI:i", 6), ("NM:i", 2), ("nM", 0)),
            (("MD:Z", "101"), ("NH:i", 7), ("HI:i", 6), ("NM:i", 0), ("nM", 0)),
            (("MD:Z", "9T87T3"), ("NH:i", 5), ("HI:i", 1), ("NM:i", 2), ("nM", 3)),
            (("MD:Z", "62T38"), ("NH:i", 5), ("HI:i", 1), ("NM:i", 1), ("nM", 3)),
            (("MD:Z", "9T87T3"), ("NH:i", 5), ("HI:i", 2), ("NM:i", 2), ("nM", 3)),
            (("MD:Z", "62T38"), ("NH:i", 5), ("HI:i", 2), ("NM:i", 1), ("nM", 3)),
            (("MD:Z", "38A62"), ("NH:i", 5), ("HI:i", 3), ("NM:i", 1), ("nM", 3)),
            (("MD:Z", "3A87A9"), ("NH:i", 5), ("HI:i", 3), ("NM:i", 2), ("nM", 3)),
            (
                ("MD:Z", "2T4G14C31A23"),
                ("NH:i", 8),
                ("HI:i", 4),
                ("NM:i", 4),
                ("nM", 7),
            ),
            (("MD:Z", "4G14C31A23"), ("NH:i", 8), ("HI:i", 4), ("NM:i", 3), ("nM", 7)),
        ]
        self.iterator = zip(
            query_names,
            flags,
            reference_ids,
            reference_starts,
            mapping_qualities,
            cigars,
            next_reference_ids,
            next_reference_starts,
            template_lengths,
            query_sequences,
            query_qualities,
            tags,
        )


@pytest.fixture(scope="session")
def get_args():
    args = mockParser()
    args.parse_args()
    return args


@pytest.fixture(scope="session")
def get_bam(tmp_path_factory, get_args):

    """\
    Create test BAM file.

    Parameters
    ----------
    tmp_path_factory
        tmp_path_factory fixture

    Returns
    -------
    :pathlib.Path: base temporary path
    """
    from pathlib import Path
    import pysam as ps

    from pulsertc.utils import index_bam_file

    data = BAMData()
    data.get_header()
    data.get_entries()

    loc = tmp_path_factory.mktemp("data")
    bamf = Path(loc, "pytest.bam")

    def _index(l, cond):
        for i, elem in enumerate(l):
            if cond(elem):
                return i
        return None

    with ps.AlignmentFile(bamf, "wb", header=data.header) as outf:
        for (
            qname,
            flag,
            rname,
            rstart,
            mapq,
            cigar,
            mrnm,
            pnext,
            tlen,
            seq,
            qual,
            tag,
        ) in data.iterator:
            a = ps.AlignedSegment()
            a.query_name = qname
            a.flag = flag
            a.reference_id = _index(data.header["SQ"], lambda x: x["SN"] == rname)
            a.reference_start = rstart
            a.mapping_quality = mapq
            a.cigarstring = cigar
            next_ref_id = rname if mrnm == "=" else mrnm
            a.next_reference_id = _index(
                data.header["SQ"], lambda x: x["SN"] == next_ref_id
            )
            a.next_reference_start = pnext
            a.template_length = tlen
            a.query_sequence = seq
            a.query_qualities = ps.qualitystring_to_array(qual)
            a.tags = tag
            outf.write(a)

    index_bam_file(bamf.as_posix(), get_args)

    return bamf.as_posix()


@pytest.fixture(scope="session")
def get_mismatches(get_bam, get_args):

    """\
    Call get_mismatches

    Parameters
    ----------
    get_bam
        BAM file

    Returns
    -------
    aligned pairs and mismatches
    """
    import pysam as ps

    from pulsertc.utils import apply_parallel_iter
    from pulsertc.splbam import get_mismatches

    num_cpus = get_args.num_cpus
    library_type = get_args.library_type
    offset = get_args.offset

    bam = ps.AlignmentFile(get_bam, "rb")
    SN = (SQ["SN"] for SQ in bam.header["SQ"])
    bam.close()

    mismatch_and_aligned_pairs = apply_parallel_iter(
        SN,
        num_cpus,
        get_mismatches,
        get_bam,
        library_type,
        offset,
        progress_bar=False,
        backend="multiprocessing",
    )

    return mismatch_and_aligned_pairs
