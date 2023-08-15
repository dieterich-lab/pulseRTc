#! /usr/bin/env python3

"""Wrapper function for splbam. 
"""

import argparse
import logging
import shlex
import yaml

from pathlib import Path

import pulsertc.utils as utils

logger = logging.getLogger(__name__)


default_num_cpus = 1
default_mem = '80G'


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="""Wrapper function for splbam""")

    parser.add_argument('config', help="The yaml configuration file (full path).")
    
    parser.add_argument('-l', '--library-type', help="""Library type: stranded, 
                        reverse-stranded, or infer. Unstranded libraries are not 
                        handled. GTF annotation file required if infer. Ignore inward/
                        outward orientation. Library type is only inferred from 
                        counts to one or the other, based on selected flags.""", type="str", 
                        choices=["stranded", "reverse", "infer"], default="infer")
    
    parser.add_argument('--gtf', help="""Annotation GTF file, required if 
                        [--library-type infer]""", type="str")
    
    parser.add_argument('-s', '--sample-size', help="""Sample size if 
                        [--library-type infer]""", type=int, default=1000)
    
    parser.add_argument('-isize', '--insert-size', help="""For quality control only, insert 
                        size is interpreted as sum of length of read 1 and read 2, ignoring
                        overlap, inner distance, etc. If None, this will be inferred from
                        the average read length.""", type=int, default=None)

    parser.add_argument('--overwrite', help="""If this flag is present, existing files 
                        will be overwritten.""", action='store_true')

    utils.add_sbatch_options(parser, num_cpus=default_num_cpus, mem=default_mem)
    utils.add_logging_options(parser)
    args = parser.parse_args()
    utils.update_logging(args)

    config = yaml.load(open(args.config), Loader=yaml.FullLoader)

    # check that all of the necessary programs are callable
    programs = ["splbam"]
    utils.check_programs_exist(programs)
    
    required_keys = ["samples",
                     "parent"]
    utils.check_keys_exist(config, required_keys)
    
    # create output directory structure
    Path(config["parent"], "mapping").mkdir(parents=True, exist_ok=True)
    Path(config["parent"], "mismatches").mkdir(parents=True, exist_ok=True)
    # if running with Salmon, this directory is created later
    Path(config["parent"], "tables", "featureCounts").mkdir(parents=True, exist_ok=True)
    
    # handle all option strings to call splbam
    # use defaults for most options (not passed)
    logging_str = utils.get_logging_options_string(args)

    all_opt_str = ""
    if 'base_qual' in config.keys():
        all_opt_str = f"-q {config['base_qual']}"
    if 'trim5p' in config.keys():
        all_opt_str = f"{all_opt_str} --trim5p {config['trim5p']}"
    if 'trim3p' in config.keys():
        all_opt_str = f"{all_opt_str} --trim3p {config['trim3p']}"
    if 'ref_base' in config.keys():
        all_opt_str = f"{all_opt_str} -ref {config['ref_base']}"
    if 'base_change' in config.keys():
        all_opt_str = f"{all_opt_str} -bc {config['base_change']}"
    
    vcf_str = ""
    if 'vcf' in config.keys() and config['vcf'] == True:
        vcf_str = "--vcf"
    if 'vcf' in config.keys() and not 'snpdata' in config.keys():
        msg = "[snpdata] is empty, ignoring [vcf]!"
        logger.warning(msg)
    snpdata_str = ""
    if 'snpdata' in config.keys():
        snpdata_str = f"-s {config['snpdata']} {vcf_str}"
    
    # handle do_not_call so that we do call splbam, but that it does not run anything
    call = not args.do_not_call
    do_not_call_str = ""
    if not call:
        do_not_call_str = "--do-not-call"
    args.do_not_call = False

    overwrite_str = ""
    if args.overwrite:
        overwrite_str = "--overwrite"
        
    mem_str = "--mem {}".format(shlex.quote(args.mem))
    
    if args.library_type == "infer":
        
        if args.gtf is None:
            msg = "Annotation file (GTF format) is required to infer library type. Terminating!"
            logger.error(msg)
            return
        
        bamf = Path(config["bamloc"], list(config['samples'].values())[0]).as_posix()
        exist = utils.check_files_exist([bamf, args.gtf], raise_on_error=True, logger=logger)
        
        msg = f"Library type inferred using {bamf}, with first {sample_size} gene records " \
              f"for each strand from {args.gtf}.\n"
        
        lib_counts = infer_library_type(bamf, gtff, sample_size=args.sample_size)
        stranded = lib_counts["stranded"]
        reverse = lib_counts["reverse"]
        total = stranded + reverse
        proportions = f"stranded {stranded/total*100}, reverse-stranded {reverse/total*100}"
        margin = 0.3*max(stranded, reverse)
        if abs(stranded-reverse) < margin:
            msg = f"Unable to infer library type: {proportions}.\n" \
                  f"Increase [--sample-size], or specify libray type. Terminating!"
            logger.error(msg)
            return
        
        args.library_type = "stranded"
        if reverse > stranded:
            args.library_type = "reverse"
            
        msg = f"Inferred library type is {args.library_type}, based on: {proportions}."
        logger.warning(msg)
    
    # TODO: pass argument lib type
    # only stranded and reverse-stranded - docs must match Salmon, featureCounts, etc.
    for name, bam in config["samples"].items():
        
        bam_path = Path(config["bamloc"], bam).as_posix()
        cmd = f"splbam {bam_path} {Path(config['parent'], 'mapping').as_posix()} " \
              f"{Path(config['parent'], 'mismatches').as_posix()} {name} " \
              f"{all_opt_str} {snpdata_str} --num-cpus {args.num_cpus} {mem_str} " \
              f"{logging_str} {overwrite_str} {do_not_call_str}"
        utils.check_sbatch(cmd, args=args)
    
    
if __name__ == '__main__':
    main()   
