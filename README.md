# pulseRTc

Starting from alignments (BAM files), **pulseRTc** can be used to infer genome-wide kinetics of RNA abundance from 4sU-tagging metabolic labeling methods, including biochemical enrichment after thiol-specific biotinylation, and recent approaches such as [SLAM-seq](https://www.nature.com/articles/nmeth.4435), [TimeLapse-seq](https://www.nature.com/articles/nmeth.4582) or [TUC-seq](https://onlinelibrary.wiley.com/doi/10.1002/anie.201707465) that rely on bioinformatic enrichment of newly transcribed RNAs.

We use [pulseR](https://dieterich-lab.github.io/pulseR/index.html), a kinetic and statistical modelling framework, to estimate RNA decay rates (half-lives).

> **Warning**\
> Kinetic models for biochemical separation (with or without spike-ins) and nucleotide conversion protocols are pre-defined for a _pulse labeling experiment_. These models can be modified to reflect differences in experimental procedures ( _e.g._ chase experiment), but this requires a bit of R scripting, see [models](scripts/R/pulser/models.R). Currently, the pre-processing workflow for nucleotide conversion protocols splits input alignment (BAM) files into "labelled" (reads with T -> C) and "unlabelled" (old/pre-existing reads without T -> C), according to the presence or not of characteristic mismatches ( after quality filtering, SNPs removal, _etc._ ).

# Getting started

Clone the git repository

```
git clone https://github.com/dieterich-lab/pulseRTc.git
cd pulseRTc
```

Create a conda environment and activate it

```
# or use mamba...
conda env create --name pulsertc --file environment.yml
conda activate pulsertc
```

> **Note**\
> General dependencies include [Samtools](http://www.htslib.org/), [BCFtools](http://samtools.github.io/bcftools/howtos/index.html), and [Picard Tools](https://broadinstitute.github.io/picard/). **Picard Tools** is optional. For gene abundance quantification, [Subread](http://subread.sourceforge.net/) (featureCounts) is also installed. See below for transcript abundance quantification. Comment out corresponding lines in [environment.yml](environment.yaml), before creating the environment, to skip installation of selected dependencies.

> **Warning**\
> For transcript abundance quantification, [GffRead](http://ccb.jhu.edu/software/stringtie/gff.shtml#gffread) and [Salmon](https://salmon.readthedocs.io/en/latest/) are required, however **Salmon** is currently NOT installable via `conda`, see [here](https://github.com/dieterich-lab/pulseRTc/issues/2). If you have existing installations, you can specify paths to binaries for both **Salmon** and **GffRead** in [makefile.vars](scripts/makefile.vars).

Install additional R packages

```
# install pulseR
Rscript scripts/R/install_extra_r_pkgs.R
```

Install **pulseRTc**. To install the local VCS project in development mode, use the `--editable` or `-e` option, otherwise
this flag can be ignored

```
pip --verbose install [-e] . 2>&1 | tee install.log
```

# Running **pulseRTc**

The workflow is made up of different parts which may or may not all be run sequentially. To get started, `cd scripts`, and set all paths, variables, and program options in [makefile.vars](scripts/makefile.vars). None of these paths is required to be located in this repository, _i.e._ they can be anywhere on your system, as long as they are accessible from the calling context. Most of them are self-explanatory, _e.g_ `INDEXLOC` is the path to your `FASTA` and `GTF` files, _etc._ You also need to define which commands `MKCMD` are used to call programs, _i.e._ whether you use Slurm, or simply bash, _etc._

We also provide examples of [SAMPLELIST.TXT](data/SAMPLELIST.TXT) (optional) and [config.yaml](data/config.yaml) (required) that are defined using variables `SAMPLELIST` and `YML`, respectively. `DATALOC` is the location where these files reside.

> **Note**\
> Make sure you have permission for all scripts! This can be done _e.g._ by `chmod +x filename`

### Pre-processing

If your BAM files are ready to be processed, the first step is to run **BCFtools** to get SNPs, see next step. In general, removing duplicate ( _e.g_ due to PCR artefacts) might help to reduce the number of false positives for variant calling.

> **Warning**\
> BAM files need to be sorted and indexed. They require the MD tag! If they do not, you can try `samtools calmd file.sorted.bam reference.fasta > file.md.sorted.bam` where the reference fasta must be the same as the one used for mapping.

### Estimate SNPs

If your BAM files are ready to be processed, then you want to get SNPs. We recommend to use **BCFtools**, as it is installed in the environment. You need to make sure `BCFLOC` and `OUTPUT_BCF` are set, and `make run-bcftools`. The final location of your BAM files is defined in the [config.yaml](data/config.yaml) using the `bamloc` key. This is automatically picked-up by [makefile.vars](scripts/makefile.vars).

### Run splbam (split BAM files into "old" and "new")

Once this is done, you need to update the `snpdata` key in the [config.yaml](data/config.yaml) with the output of **BCFtools** (or **GRAND-SLAM**), and set `vcf: True` if using the **BCFtools** output. You need to set `parent` and any program options, and `make run-workflow`. After completion, you can get some summary statistics using `make plot-mm`.

### Abundance estimation - gene

We are now ready to proceed to read abundance quantification, as **pulseR** needs read counts as input. If we are interested in gene read counts, then we can use **featureCounts**, and `make count-fc`. For transcript quantification, see further below.

### pulseR

To prepare **pulseR** input from **featureCounts** tables, we `make prep-fc`, then fit the models using `make pulse-all-fc-counts`. Before fitting the models, selected time points have to be defined in [time_pts.R](scripts/R/pulser/time_pts.R).

> **Note**\
> `make pulse-all-fc-counts` will fit all time points defined in [time_pts.R](scripts/R/pulser/time_pts.R) using _allSets_. After these results are available, you can fit subsets of time points using _timeSets_ and `make pulse-set-fc-counts`.

> **Warning**\
> The script [prep_fc.R](scripts/R/pulser/prep_fc.R) called using `make prep-fc` is NOT all-purpose, and may need to be adjusted depending on your data. If you already have count tables and associated metadata (R objects), they can be put under `pulsedir/featureCounts/data`, where `pulsedir` is given in the [config.yaml](data/config.yaml). The count tables must be named `counts.rds` (gene by samples), and the metadata `samples.rds`, with columns `sample fraction time rep`, where `sample` must match the columns names of `counts.rds`.

If you have matching results from **GRAND-SLAM**, you can combine all estimates from **pulseR** and **GRAND-SLAM** together by `make match-gs`.

### Abundance estimation - transcript

We use **Salmon** for transcript abundance quantification (and **GffRead** to extract transcript sequences). Variables `INDEX_SALMON`, `SALMON_PATH`, `SALMON_OPTS`, and `GFFREAD_PATH` need to be set in [makefile.vars](scripts/makefile.vars). We use **Salmon** in (quasi-mapping) mapping-based mode with a decoy-aware transcriptome (using the entire genome): the indexing step is independent of the reads, and only needs to be run once for a particular set of reference transcripts. This can be done with `get-salmon-index`. This generally results in larger index files, see [here](https://github.com/COMBINE-lab/salmon) for other options. If you already have an index, you can skip this part, but make sure to set `INDEX_SALMON` to the parent directory of your index files, which must be located in a sub-directory called **index**.

Since we already have _split_ BAM files, we only need to convert them back to FASTQ. Salmon does not currently have built-in support for interleaved FASTQ files, so READ1 and READ1 (paired-end sequencing reads) flags are directed to different files. Since our BAM files were coordinate sorted, we shuffle and groups reads beforehand. This is done with `make run-bam2fastq`, followed by `make run-salmon-qm`.

### Downstrean analyses

Some scripts are available _e.g._ to compare **pulseR** and **GRAND-SLAM** results, or perform DGE to check gene expression at varying labeling times _vs._ 0 (unlabelled), however many of these currently contain hard coded parameters/options.

# More explanations

### Notes

We use the terminology _labelled_ (new) and _unlabelled_ (old), although this can be inadequate, _e.g._ a read can be new without being labelled (few Ts present, low incorporation, _etc._).

### SNPs estimation

For SNPs estimation, we recommend to use all samples, including 0h time points.

### Read mapping and processing (splbam)

Mapped reads must have the MD attribute (_e.g._ for STAR `--outSAMattributes MD`). This is necessary to find mismatches. Reads are processed using the entire query sequence, including soft-clipped bases. You can use `--trim5p` and `--trimp3p` to remove mismatches at the ends of the reads, but note that this is done based on the mapped sequence.

Reads are sorted into first and second in pair for quality control. In general, one read is sense, the other one is antisense. We define a read pair to be sense, if the first read is sense, and antisense if the first read is antisense. For _stranded_ libraries, sense pairs map to the + strand, and vice versa. For _reverse-stranded_ libraries, antisense pairs map to the + strand, and vice versa.

To sort reads, all mismatches are reported with respect to the + strand, _e.g._ a T->C on the + RNA template is reported as T->C, and an A->G on the - RNA template is also reported as T->C. However, for quality control, aligned pairs are reported consistently with the read orientation/protocol at sequencing.
