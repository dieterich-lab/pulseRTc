
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


**Note:** Dependencies include [Subread](http://subread.sourceforge.net/) (featureCounts) for read quantification, [Samtools](http://www.htslib.org/), and [BCFtools](http://samtools.github.io/bcftools/howtos/index.html). **Samtools** can be used to pre-process BAM files and/or for downstream analyses, but in many cases **pysam** will be sufficient, assuming your BAM files are ready to be processed. If installing both, make sure to install compatible versions (see [environment.yml](environment.yml)). 

**Note:** [Salmon](https://salmon.readthedocs.io/en/latest/) is used for transcript abundance quantification, however it is NOT currently installable via `conda`. This has been a recurring issue and seems to be related to conflicting libraries, see *e.g.* [issue 147](https://github.com/COMBINE-lab/salmon/issues/147), or [issue 594](https://github.com/COMBINE-lab/salmon/issues/594). We currently use **Salmon 1.9.0** installed via environment modules. For standard gene read count quantification, **featureCounts** is used instead, and thus **Salmon** does not need to be installed!

**Note:** Some pre-processing scripts use [Picard Tools](https://broadinstitute.github.io/picard/), but it is not installed by default, as only older versions are available on `conda` (we use version 2.5.0). In most cases, if your BAM files are ready to be processed, you won't need it. 

**Note:** Extra **R ** packages are optional, and include standard packages for DGE, *etc*. The `all` option is currently broken, as packages such as **DESeq2** fail to install in the conda environment. If you already have a system-wide (or managed via environment modules) **R ** installation with extra packages, the search path can be updated for an interactive session (or added to scripts) using *e.g.*

```
.libPaths( c( .libPaths(), "/beegfs/homes/eboileau/R/x86_64-pc-linux-gnu-library/4.2", "/beegfs/biosw/R/4.2.1_deb11/lib/R/library") )
```

or `R_LIBS_USER` can be modified using *Renviron* in the source directory of the environment.



