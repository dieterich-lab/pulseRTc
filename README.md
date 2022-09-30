



git clone --recursive https://github.com/dieterich-lab/single-cell-nanopore.git
cd single-cell-nanopore

conda env create --name scNapBar --file environment.yaml
conda activate scNapBar


conda env create --prefix /prj/rpbp-dev/working-envs/install-conda -f environment.yml
conda activate /prj/rpbp-dev/working-envs/install-conda


TODO:

- install base env with environment.yml
- activate, and if wan to do additonal ananlyses, DGE, etc, then
conda install --yes --file r_requirements.yml
no activate and just run script



For strict version control, we could set --no-update-deps. Conda attempts to install the newest versions of the requested packages. To accomplish this, it may update some packages that are already installed, or install additional packages. To prevent existing packages from updating, use the --freeze-installed option. This may force conda to install older versions of the requested packages, and it does not prevent additional dependency packages from being installed.





For GS

java version "1.8.0_181"
Java(TM) SE Runtime Environment (build 1.8.0_181-b13)
Java HotSpot(TM) 64-Bit Server VM (build 25.181-b13, mixed mode)

do we install java?



Save environment settings
conda env export -p ${path_to_env} > ${out_path}/conda_environment_dge_BioC3.8.full.yml


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










First create a virtual environment:
 
```
python3 -m venv /path/to/virtual/environment
```

For information about Python virtual environments, see the [venv](https://docs.python.org/3/library/venv.html) documentation.
To activate the new virtual environment and install `pulseRTc`:

```
# Activate the new virtual environment.
source /path/to/virtual/environment/bin/activate

# If necessary, upgrade pip and wheel or additional packages (such as setuptools if installing in editable mode).
pip install --upgrade pip setuptools wheel



```

### Running

To run the workflow, call the main wrapper `run` with a configuration file. See *makefile* for examples.

Scripts and example configuration files are under *run*. Intermediate results are under *workflow*.
The pulseR scripts are under *pulser*, and the final results under *results*.

The pulseR code and analysis scripts are in R. Minimal requirements are


```{r}
install.packages("devtools")
library(devtools)

install_github("dieterich-lab/pulseR", subdir="pkg")
install_cran("tidyverse")

``` 

**Note:** The directory *workflow/mapping* does not contain the alignment files. Mapping is left to the user.
Some large intermediate files are missing from *workflow/mismatches*.




*************************


example worlfow pipeline set-by-step from makefile, and scritps explanation below
explain also files, etc.

whole pre-processing is not necessairy applicable, in general you just need to make sure
you have a list of BAM files ready for processing, sorted, dedup, etc. then these only need
to be listed in confiog.yml for the workflow, etc.
So you can skip this part (list makefile commands and scripts), or you can adjust them to your need.

if using picard AND GS, make sure java version are not conlfucting....


SAMPLELIST.TXT

*************

scripts

sortnrename
custome script assume a directory structure and naming convention, script might nend to be adjusted!

markdups


