
pulseRTc
========

Starting from alignments (BAM files), **pulseRTc** can be used to infer genome-wide kinetics of RNA abundance from 4sU-tagging metabolic labeling methods, including biochemical enrichment after thiol-specific biotinylation, and recent approaches such as SLAM-seq, TimeLapse-seq or TUC-seq that rely on bioinformatic enrichment of newly transcribed RNAs.

We use [pulseR](https://dieterich-lab.github.io/pulseR/index.html), a kinetic and statistical modelling framework, to estimate RNA decay rates (half-lives). 
To handle data arising from nucleotide conversion protocols, alignment files are split into old/pre-existing (unlabeled) and new (metabolically labeled, T -> C) reads. 

**Note:** Kinetic models for biochemical separation (with or without spike-ins) and nucleotide conversion protocols are pre-defined for a *pulse labeling experiment*. These models can be modified to reflect differences in experimental procedures ( *e.g.* chase experiment), but this requires a bit of R scripting, see under [pulser](pulser/models.R). Currently, the pre-processing workflow for nucleotide conversion protocols splits input BAM files into "labeled" (reads with T -> C) and "unlabeled" (reads without T -> C), according to the presence or not of characteristic mismatches ( after quality filtering, SNPs removal, *etc.* ).

**Note:** [GRAND-SLAM](https://github.com/erhard-lab/gedi/wiki/GRAND-SLAM) can be used as another method to estimate old and new RNAs (for RT *e.g.* SLAM-seq or direct conversion *e.g.* TUC- or TimeLapse-seq). It is used by default for SNP calling, and requires **Java 1.8**. As **GRAND-SLAM** is currently only available under license agreement, if you want to use it, you need to install it yourself! We do not know which version is currently available, but according to the latest documention **Java 1.8** is still a requirement. Alternatively, for SNP calling **BCFtools** can be used (installed as a dependency by default). 


Getting started
===============

Clone the git repository

```
git clone https://github.com/dieterich-lab/pulseRTc.git
cd pulseRTc
```

Create a conda environment and activate it

```
conda env create --name pulsertc --file environment.yml
conda activate pulsertc
```

Install additional R packages

```
# use "all" to install additional packages for plotting, DGE, GSEA, etc.
# otherwise only pulseR is installed 
Rscript scripts/R/install_extra_r_pkgs.R [all]
```

Install the package. To install the local VCS project in development mode, use the `--editable` or `-e` option, otherwise
this flag can be ignored

```
pip --verbose install [-e] . 2<&1 | tee install.log
```

> **Note**
> This is a note

> **Warning**
> This is a warning


<div class="panel panel-info">
**Note**
{: .panel-heading}
<div class="panel-body">

NOTE DESCRIPTION

**Note:** Dependencies include [Subread](http://subread.sourceforge.net/) (featureCounts) for read quantification, [Samtools](http://www.htslib.org/), and [BCFtools](http://samtools.github.io/bcftools/howtos/index.html). **Samtools** can be used to pre-process BAM files and/or for downstream analyses, but in many cases **pysam** will be sufficient, assuming your BAM files are ready to be processed. If installing both, make sure to install compatible versions (see [environment.yml](environment.yml)). 

</div>
</div>


**Note:** Dependencies include [Subread](http://subread.sourceforge.net/) (featureCounts) for read quantification, [Samtools](http://www.htslib.org/), and [BCFtools](http://samtools.github.io/bcftools/howtos/index.html). **Samtools** can be used to pre-process BAM files and/or for downstream analyses, but in many cases **pysam** will be sufficient, assuming your BAM files are ready to be processed. If installing both, make sure to install compatible versions (see [environment.yml](environment.yml)). 

**Note:** [Salmon](https://salmon.readthedocs.io/en/latest/) is used for transcript abundance quantification, however it is NOT currently installable via `conda`. This has been a recurring issue and seems to be related to conflicting libraries, see *e.g.* [issue 147](https://github.com/COMBINE-lab/salmon/issues/147), or [issue 594](https://github.com/COMBINE-lab/salmon/issues/594). We currently use **Salmon 1.9.0**. In addition, we use [GffRead](http://ccb.jhu.edu/software/stringtie/gff.shtml) to extract transcript sequences (not installed by default). Alternatively, you can specify paths to binaries for both **Salmon** and **GffRead**. For gene read count quantification, **featureCounts** is used instead, and thus **Salmon** and **GffRead** do not need to be installed! 

**Note:** Some pre-processing scripts use [Picard Tools](https://broadinstitute.github.io/picard/), but it is not installed by default, as only older versions are available on `conda` (we use version 2.5.0). In most cases, if your BAM files are ready to be processed, you won't need it. 

**Note:** Extra **R ** packages are optional, and include standard packages for DGE, *etc*. The `all` option is currently broken, as packages such as **DESeq2** fail to install in the conda environment. If you already have a system-wide (or managed via environment modules) **R ** installation with extra packages, the search path can be updated for an interactive session (or added to scripts) using *e.g.*

```
.libPaths( c( .libPaths(), "/beegfs/homes/eboileau/R/x86_64-pc-linux-gnu-library/4.2", "/beegfs/biosw/R/4.2.1_deb11/lib/R/library") )
```

or `R_LIBS_USER` can be modified using *Renviron* in the source directory of the environment.


****************


Salmon options


- The alignment-based mode of Salmon does not require indexing. Rather, you can simply provide Salmon with a FASTA file of the transcripts and a SAM/BAM file containing the alignments you wish to use for quantification.

you can provide salmon with your aligned reads. If you do this, be sure that the reads are not sorted by position, and that all alignments for the same read appear consecutively in the alignment file.

Salmon, like eXpress 1, uses a streaming inference method to perform transcript-level quantification. One of the fundamental assumptions of such inference methods is that observations (i.e. reads or alignments) are made “at random”. This means, for example, that alignments should not be sorted by target or position. 

Does this also applies to quasi-mapping mode, i.e. reads in FASTQ?

If the fastqs come right from the sequencer, you are fine. The thing is that you can transform BAM back to fastq for realignments/requantification, and as BAMs are often coordiate-sorted, the resulting fastq would not be randomly ordered. Therefore the recommendation is to shuffle fastq prior to quantification.

- The (quasi-mapping) mapping-based mode of Salmon runs in two phases; indexing and quantification.
The indexing step is independent of the reads, and only needs to be run once for a particular set of reference transcripts.


SO FAR, IMPLEMENT THIS ONE ONLY, THSI AVOIDS TO MAP BACK FILES, INSTALL e.g. STAR, etc.



Preparing transcriptome indices (mapping-based mode)


Assume that transcripts.fa contains the set of transcripts you wish to quantify. We generally recommend that you build a decoy-aware transcriptome file.

The first is to compute a set of decoy sequences by mapping the annotated transcripts you wish to index against a hard-masked version of the organism’s genome. This can be done with e.g. MashMap2, and we provide some simple scripts to greatly simplify this whole process. Specifically, you can use the generateDecoyTranscriptome.sh script, whose instructions you can find in this README.

https://github.com/COMBINE-lab/SalmonTools/blob/master/scripts/generateDecoyTranscriptome.sh
https://github.com/COMBINE-lab/SalmonTools/blob/master/README.md

The second is to use the entire genome of the organism as the decoy sequence. This can be done by concatenating the genome to the end of the transcriptome you want to index and populating the decoys.txt file with the chromosome names. Detailed instructions on how to prepare this type of decoy sequence is available here. This scheme provides a more comprehensive set of decoys, but, obviously, requires considerably more memory to build the index.

> NEED TO INSTALL GFFREAD IF USING SALMON, CAN BE DONE VIA CONDA, BUT CURRENTLY NOT DONE (USE MODULES!)






Salmon does not currently have built-in support for interleaved FASTQ files (i.e., paired-end files where both pairs are stored in the same file).

ONLY PAIRED-END!




    Counting number of reads per gene. With --quantMode GeneCounts option STAR will count number reads per gene while mapping. A read is counted if it overlaps (1nt or more) one and only one gene. Both ends of the pairedend read are checked for overlaps. The counts coincide with those produced by htseq-count with default parameters. This option requires annotations (GTF or GFF with –sjdbGTFfile option) used at the genome generation step, or at the mapping step. STAR outputs read counts per gene into ReadsPerGene.out.tab file with 4 columns which correspond to different strandedness options:

    column 1: gene ID

    column 2: counts for unstranded RNA-seq

    column 3: counts for the 1st read strand aligned with RNA (htseq-count option -s yes)

    column 4: counts for the 2nd read strand aligned with RNA (htseq-count option -s reverse)

    Select the output according to the strandedness of your data. Note, that if you have stranded data and choose one of the columns 3 or 4, the other column (4 or 3) will give you the count of antisense reads. With --quantMode TranscriptomeSAM GeneCounts, and get both the Aligned.toTranscriptome.out.bam and ReadsPerGene.out.tab outputs.

    
Then when implemented, we could use directly transcriptome if unsorted...

