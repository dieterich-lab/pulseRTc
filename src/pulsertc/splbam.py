#! /usr/bin/env python3


"""Split a BAM file containing alignments
from a SLAM/TUC/TL-seq experiment into
labelled/unlabelled BAM files, to use as input
for PulseR (counts). Additionally outputs info
on mismatch rates, etc.

NOTE: We use the notation labelled/unlabelled, but
      new/old (RNA) would be less ambiguous, cf.
      biochemical separation.

NOTE: Currently using GRAND-SLAM snpdata - default
      format (P Value = 0.001)
"""

import os
import sys
import logging
import argparse
import shlex
import pandas as pd
import numpy as np
import pysam as ps

from collections import defaultdict

import pulsertc.utils as utils

logger = logging.getLogger(__name__)


default_num_cpus = 1
default_mem = "80G"


def get_mismatches(contig, bam_file, lib_type, offset):
    """
    Find mismatches and report all aligned pairs
    """

    bam = ps.AlignmentFile(bam_file, "rb").fetch(contig=contig)
    bed_entries = []
    aligned_pairs_dict = defaultdict(int)

    for r in bam:
        rna_template = "+"
        if lib_type == "stranded" and (
            (r.is_read1 and r.is_reverse) or (r.is_read2 and r.is_forward)
        ):
            rna_template = "-"
        if lib_type == "reverse" and (
            (r.is_read1 and r.is_forward) or (r.is_read2 and r.is_reverse)
        ):
            rna_template = "-"

        # only for estimation of conversion rates
        # read filtering done during abundance estimation (featureCounts and/or Salmon options)
        used = utils.is_used(r)

        # only matched bases are returned - no None on either side
        aligned_pairs = r.get_aligned_pairs(with_seq=True, matches_only=True)
        for aligned_pair in aligned_pairs:
            mpos = aligned_pair[0]
            ref = aligned_pair[2].upper()
            base = r.query_sequence[aligned_pair[0]]
            # mismatches in lowercase
            if aligned_pair[2].islower():
                entry = {
                    "contig": contig,
                    "start": aligned_pair[1],
                    "end": aligned_pair[1] + 1,
                    "name": r.query_name,
                    "score": used,
                    "strand": rna_template,
                    "base": base
                    if rna_template == "+"
                    else base.translate(str.maketrans("ACGTacgtNnXx", "TGCAtgcaNnXx")),
                    "ref": ref
                    if rna_template == "+"
                    else ref.translate(str.maketrans("ACGTacgtNnXx", "TGCAtgcaNnXx")),
                    "base_qual": r.query_qualities[aligned_pair[0]],
                    "m_pos": mpos,
                    "rlen": r.query_length,
                    "read1": r.is_read1,
                    "samflag": r.flag,
                }
                bed_entries.append(entry)

            # report aligned pairs in a way that is consistent
            # with the read orientation
            if used:
                if r.is_reverse:
                    ref = ref.translate(str.maketrans("ACGTacgtNnXx", "TGCAtgcaNnXx"))
                    base = base.translate(str.maketrans("ACGTacgtNnXx", "TGCAtgcaNnXx"))
                # shift position of read 2
                if r.is_read2:
                    mpos += offset
                key = "{}:{}:{}".format(ref, base, mpos)
                aligned_pairs_dict[key] += 1

    bed = pd.DataFrame(bed_entries)

    return (aligned_pairs_dict, bed)


def get_mismatch_details(details, filename1, filename2, offset):
    """
    Wrangle all aligned pairs per read for quality control
    """

    # first sort mismatch details list of dict into one dict
    mismatch_details = defaultdict(int)
    for d in details:
        for key, val in d.items():
            mismatch_details[key] += val

    # get coverage
    coverage = defaultdict(int)
    for key, val in mismatch_details.items():
        genomic, read, pos = key.split(":")
        ckey = "{}:{}".format(genomic, pos)
        if genomic == "N" or read == "N":
            continue
        coverage[ckey] += val

    # and convert to final format
    sl = []
    for key, val in mismatch_details.items():
        genomic, read, pos = key.split(":")
        if genomic == read or genomic == "N" or read == "N":
            continue
        cov = coverage["{}:{}".format(genomic, pos)]
        e = {
            "Category": "Any",
            "Genomic": genomic,
            "Read": read,
            "Position": int(pos),
            "Coverage": cov,
            "Mismatches": val,
        }
        sl.append(e)
    mismatch_details_df = pd.DataFrame(sl)
    mismatch_details_df = mismatch_details_df[
        ["Category", "Genomic", "Read", "Position", "Coverage", "Mismatches"]
    ]
    mismatch_details_df.sort_values(by=["Genomic", "Read", "Position"], inplace=True)
    mismatch_details_df.to_csv(filename1, sep="\t", index=False, compression="gzip")

    # now sum everything to get overall` mismatches, but split by read
    mismatch_details_df.loc[
        mismatch_details_df["Position"] < offset, "Orientation"
    ] = "First"
    mismatch_details_df.loc[
        mismatch_details_df["Position"] >= offset, "Orientation"
    ] = "Second"

    cov = (
        mismatch_details_df[["Genomic", "Orientation", "Coverage"]]
        .drop_duplicates()
        .groupby(["Genomic", "Orientation"])
        .sum()
    )
    cov = pd.DataFrame(cov.to_records())

    mis = (
        mismatch_details_df[["Genomic", "Read", "Orientation", "Mismatches"]]
        .groupby(["Genomic", "Read", "Orientation"])
        .sum()
    )
    mis = pd.DataFrame(mis.to_records())

    mismatches_counts = mis.merge(cov, how="left", on=["Genomic", "Orientation"])
    mismatches_counts["Category"] = "Any"
    mismatches_counts = mismatches_counts[
        ["Category", "Orientation", "Genomic", "Read", "Coverage", "Mismatches"]
    ]
    mismatches_counts.sort_values(by=["Orientation", "Genomic", "Read"], inplace=True)

    mismatches_counts.to_csv(filename2, sep="\t", index=False, compression="gzip")

    return mismatches_counts


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Split a BAM file containing alignments from
                                     a SLAM/TUC/TL-seq experiment into labelled/unlabelled BAM files.
                                     Requires the MD tag. All reads are kept incl. unmapped,
                                     secondary/supplementary, etc. these can be filtered out when
                                     counting (featureCounts/Salmon). Optionally, SNPs can be subtracted.
                                     This scripts also does a lot of wrangling to estimate conversion rates.""",
    )

    parser.add_argument("bam", help="""The input BAM file (full path).""")

    parser.add_argument("outdir_bam", help="""The output directory (BAM files).""")

    parser.add_argument(
        "outdir_mm", help="""The output directory (mismatch information)."""
    )

    parser.add_argument("name", help="""The output base name without extension.""")

    parser.add_argument(
        "-l",
        "--library-type",
        help="""Library type: stranded,
                        reverse-stranded, or infer. Unstranded libraries are not
                        handled. GTF annotation file required if infer. Ignore inward/
                        outward orientation. Library type is only inferred from
                        counts to one or the other, based on selected flags.""",
        type=str,
        choices=["stranded", "reverse", "infer"],
        default="infer",
    )

    parser.add_argument(
        "-a",
        "--gtf",
        help="""Annotation GTF file, required if
                        [--library-type infer]""",
        type=str,
    )

    parser.add_argument(
        "-ssize",
        "--sample-size",
        help="""Sample size if
                        [--library-type infer]""",
        type=int,
        default=1000,
    )

    parser.add_argument(
        "-isize",
        "--insert-size",
        help="""For quality control only, insert
                        size is interpreted as sum of length of read 1 and read 2, ignoring
                        overlap, inner distance, etc. If None, this will be inferred from
                        the average read length.""",
        type=int,
        default=None,
    )

    parser.add_argument(
        "-s",
        "--subtract",
        help="""SNPs to be subtracted (GS default format)""",
        type=str,
    )

    parser.add_argument(
        "--vcf", help="""Use this flag if SNPs are in VCF format""", action="store_true"
    )

    parser.add_argument(
        "-ref",
        "--ref-base",
        help="""Conversion reference base.""",
        choices=["A", "C", "G", "T"],
        default="T",
    )

    parser.add_argument(
        "-bc",
        "--base-change",
        help="""Conversion base (substitution/mismatch).""",
        choices=["A", "C", "G", "T"],
        default="C",
    )

    parser.add_argument(
        "-q",
        "--base-qual",
        help="The minimum base quality for any given mismatch (default: 20).",
        type=int,
        default=20,
    )

    parser.add_argument(
        "--trim5p",
        help="The number bases to trim at the 5' ends of reads (default: 0).",
        type=int,
        default=0,
    )

    parser.add_argument(
        "--trim3p",
        help="The number bases to trim at the 3' ends of reads (default: 0).",
        type=int,
        default=0,
    )

    parser.add_argument(
        "--overwrite",
        help="""If this flag is present, then existing files
        will be overwritten.""",
        action="store_true",
    )

    parser.add_argument(
        "-t",
        "--tmp",
        help="""Optional argument: where to write
        temporary files. If not specified, programs-specific tmp will be used.""",
        default=None,
    )

    utils.add_sbatch_options(parser, num_cpus=default_num_cpus, mem=default_mem)
    utils.add_logging_options(parser)
    args = parser.parse_args()
    utils.update_logging(args)

    msg = "[splbam]: {}".format(" ".join(sys.argv))
    logger.info(msg)

    # if using slurm, submit the script
    if args.use_slurm:
        cmd = "{}".format(" ".join(shlex.quote(s) for s in sys.argv))
        utils.check_sbatch(cmd, args=args)
        return

    # check output path
    exist = utils.check_files_exist(
        [args.outdir_bam, args.outdir_mm], raise_on_error=True, logger=logger
    )

    # check that all files exist
    input_files = [args.bam]
    if args.subtract:
        input_files.append(args.subtract)

    if args.library_type == "infer":

        if args.gtf is None:
            msg = "Annotation file (GTF format) is required to infer library type. Terminating!"
            logger.error(msg)
            return
        input_files.append(args.gtf)

    exist = utils.check_files_exist(input_files, raise_on_error=True, logger=logger)

    # create the output files - BAM
    labelled_filename = "{}.labelled.unsrt.bam".format(args.name)
    labelled_filename = os.path.join(args.outdir_bam, labelled_filename)

    sorted_labelled_filename = "{}.labelled.bam".format(args.name)
    sorted_labelled_filename = os.path.join(args.outdir_bam, sorted_labelled_filename)

    unlabelled_filename = "{}.unlabelled.unsrt.bam".format(args.name)
    unlabelled_filename = os.path.join(args.outdir_bam, unlabelled_filename)

    sorted_unlabelled_filename = "{}.unlabelled.bam".format(args.name)
    sorted_unlabelled_filename = os.path.join(
        args.outdir_bam, sorted_unlabelled_filename
    )

    # create the output files - mismatches
    mismatch_details_filename = "{}.mismatchDetails.tab.gz".format(args.name)
    mismatch_details_filename = os.path.join(args.outdir_mm, mismatch_details_filename)

    mismatch_filename = "{}.mismatches.tab.gz".format(args.name)
    mismatch_filename = os.path.join(args.outdir_mm, mismatch_filename)

    mismatch_filename_final = "{}.mismatches-used.tab.gz".format(args.name)
    mismatch_filename_final = os.path.join(args.outdir_mm, mismatch_filename_final)

    # and check if exist
    out_files = [
        labelled_filename,
        sorted_labelled_filename,
        unlabelled_filename,
        sorted_unlabelled_filename,
        mismatch_details_filename,
        mismatch_filename,
        mismatch_filename_final,
    ]
    all_out_exists = all([os.path.exists(of) for of in out_files])
    if args.overwrite or not all_out_exists:
        pass
    else:
        msg = "All output files {} already exist. Skipping call.".format(out_files)
        logger.warning(msg)
        return

    # infer library type
    if args.library_type == "infer":
        msg = (
            f"Library type inferred using {args.bam}, with first {args.sample_size} gene records "
            f"for each strand from {args.gtf}.\n"
        )

        args.library_type = utils.infer_library_type(
            args.bam, args.gtf, sample_size=args.sample_size
        )
        if args.library_type is None:
            return

    # infer read length - for quality control
    if args.insert_size is None:
        args.insert_size = (
            utils.infer_read_length(args.bam, sample_size=args.sample_size) * 2
        )

    msg = "Getting all mismatches"
    logger.info(msg)

    # first get all mismatches
    bam = ps.AlignmentFile(args.bam, "rb")
    SN = (SQ["SN"] for SQ in bam.header["SQ"])
    bam.close()

    mismatch_and_aligned_pairs = utils.apply_parallel_iter(
        SN,
        args.num_cpus,
        get_mismatches,
        args.bam,
        args.library_type,
        args.insert_size / 2,
        progress_bar=False,
        backend="multiprocessing",
    )

    # use all aligned pairs to estimate conversion rates/mismatch ratios
    # this should be rewritten...
    aligned_pairs = [a for a, b in mismatch_and_aligned_pairs]
    mismatch_count = get_mismatch_details(
        aligned_pairs, mismatch_details_filename, mismatch_filename
    )

    all_mismatches = [b for a, b in mismatch_and_aligned_pairs]
    all_mismatches = pd.concat(all_mismatches)

    # get the conversion of interest
    m_ref = all_mismatches["ref"] == args.ref_base
    m_bc = all_mismatches["base"] == args.base_change
    all_mismatches = all_mismatches[m_ref & m_bc]

    # filter base quality - no offset of 33 needs to be subtracted
    m_qual = all_mismatches["base_qual"] >= args.base_qual

    # below we keep track of discarded mismatches to adjust rates
    discarded = all_mismatches[~m_qual].copy()
    m = discarded["read1"] & discarded["score"]
    discarded_first = discarded[m].shape[0]
    m = ~discarded["read1"] & discarded["score"]
    discarded_second = discarded[m].shape[0]

    all_mismatches = all_mismatches[m_qual]

    # discard mismatches found at read ends
    m_trim5p = all_mismatches["m_pos"] < args.trim5p
    m_trim3p = all_mismatches["m_pos"] > (all_mismatches["rlen"] - args.trim3p)
    discarded = all_mismatches[m_trim5p | m_trim3p].copy()
    m = discarded["read1"] & discarded["score"]
    discarded_first += discarded[m].shape[0]
    m = ~discarded["read1"] & discarded["score"]
    discarded_second += discarded[m].shape[0]

    all_mismatches = all_mismatches[~m_trim5p & ~m_trim3p]

    # remove all SNPs from what remains
    if args.subtract:
        # currently GRAND-SLAM snpdata default format
        if args.vcf:
            snps = utils.fmt_convert(args.subtract)
        else:
            snps = pd.read_csv(args.subtract, sep="\t")
        snps = snps.Location.unique()
        # add field
        all_mismatches["Location"] = all_mismatches[["contig", "start"]].apply(
            lambda x: ":".join([str(s) for s in x]), axis=1
        )

        discarded = all_mismatches[all_mismatches.Location.isin(snps)].copy()
        m = discarded["read1"] & discarded["score"]
        discarded_first += discarded[m].shape[0]
        m = ~discarded["read1"] & discarded["score"]
        discarded_second += discarded[m].shape[0]

        all_mismatches = all_mismatches[~all_mismatches.Location.isin(snps)]

    # adjust final mismatch counts
    # if stranded, read 1 is in the same direction as RNA template, and vice versa
    m_read1 = (mismatch_count.Genomic == args.ref_base) & (
        mismatch_count.Read == args.base_change
    )
    m_read2 = (
        mismatch_count.Genomic
        == args.ref_base.translate(str.maketrans("ACGTacgtNnXx", "TGCAtgcaNnXx"))
    ) & (
        mismatch_count.Read
        == args.base_change.translate(str.maketrans("ACGTacgtNnXx", "TGCAtgcaNnXx"))
    )
    if args.library_type == "reverse":
        m_read2 = (mismatch_count.Genomic == args.ref_base) & (
            mismatch_count.Read == args.base_change
        )
        m_read1 = (
            mismatch_count.Genomic
            == args.ref_base.translate(str.maketrans("ACGTacgtNnXx", "TGCAtgcaNnXx"))
        ) & (
            mismatch_count.Read
            == args.base_change.translate(str.maketrans("ACGTacgtNnXx", "TGCAtgcaNnXx"))
        )

    m = (mismatch_count.Orientation == "First") & m_read1
    mismatch_count.loc[m, "Mismatches"] = (
        mismatch_count.loc[m, "Mismatches"] - discarded_first
    )
    n = (mismatch_count.Orientation == "Second") & m_read2
    mismatch_count.loc[n, "Mismatches"] = (
        mismatch_count.loc[n, "Mismatches"] - discarded_second
    )
    mismatch_count = mismatch_count[m | n]
    mismatch_count.to_csv(
        mismatch_filename_final, sep="\t", index=False, compression="gzip"
    )

    # what remains are true conversions, other reads are classified as unlabelled
    # NOTE: we keep all query_name for which at least one read has a mismatch
    # this include read pairs, but also multi-mapping reads...
    # in practice, we filter them at the abundance estimation step...
    true_conversions = all_mismatches.name.unique()

    # now split the BAM file
    msg = "Reading the alignments and splitting the input BAM file"
    logger.info(msg)

    # requires a lot of memory however...
    bam = ps.AlignmentFile(args.bam, "rb")
    qname_index = ps.IndexedReads(bam)
    qname_index.build()

    # we first "split" the query names, sort by query name and write each file in turn
    # we don't sort the lists, this will not be faster, the index is just fine
    true_conversions = set(true_conversions)
    all_qnames = set([a.query_name for a in bam.fetch(until_eof=True)])
    all_qnames = all_qnames - true_conversions

    labelled = ps.AlignmentFile(labelled_filename, "wb", template=bam)
    unlabelled = ps.AlignmentFile(unlabelled_filename, "wb", template=bam)

    # labelled/new
    for qname in true_conversions:
        alignments = qname_index.find(qname)
        for a in alignments:
            labelled.write(a)
    labelled.close()

    # unlabelled/old
    for qname in all_qnames:
        alignments = qname_index.find(qname)
        for a in alignments:
            unlabelled.write(a)
    unlabelled.close()

    bam.close()

    # create the bamtools index if it does not already exists
    args.num_cpus = 6  # limit... otherwise this is problematic?!
    args.keep_intermediate_files = False  # delete unsorted bam file

    utils.sort_bam_file(labelled_filename, sorted_labelled_filename, args)
    utils.index_bam_file(sorted_labelled_filename, args)

    utils.sort_bam_file(unlabelled_filename, sorted_unlabelled_filename, args)
    utils.index_bam_file(sorted_unlabelled_filename, args)


if __name__ == "__main__":
    main()
