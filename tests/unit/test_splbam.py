class BAMData:

    # def index(l, cond):
    # for i, elem in enumerate(l):
    # if cond(elem):
    # return i
    # return None

    def get_header():
        header = {
            "HD": {"VN": "1.5"},
            "SQ": [
                {"LN": 248956422, "SN": "1"},
                {"LN": 101991189, "SN": "2"},
                {"LN": 156040895, "SN": "3"},
            ],
        }
        return header

    def get_entries():
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
        return zip(
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


# tmpfilename = "/scratch/eboileau/testbam.bam"
# with ps.AlignmentFile(tmpfilename, "wb", header=header) as outf:

# for qname, flag, rname, rstart, mapq, cigar, mrnm, pnext, tlen, seq, qual, tag in entries:

# a = ps.AlignedSegment()
# a.query_name = qname
# a.flag = flag
# a.reference_id = index(
# header["SQ"], lambda x: x["SN"] == rname
# )
# a.reference_start = rstart
# a.mapping_quality = mapq
# a.cigarstring = cigar
# next_ref_id = rname if mrnm == "=" else mrnm
# a.next_reference_id = index(
# header["SQ"], lambda x: x["SN"] == next_ref_id
# )
# a.next_reference_start = pnext
# a.template_length = tlen
# a.query_sequence = seq
# a.query_qualities = ps.qualitystring_to_array(qual)
# a.tags = tag
# outf.write(a)
