#!/usr/bin/env Rscript

# Install required packages (pulseR) and additional
# packages for DGE, GSEA, etc. 

args = commandArgs(trailingOnly = TRUE)

# set mirror
options(repos = c(CRAN = "https://cran.uni-muenster.de/"))

# always install pulseR
devtools::install_github("dieterich-lab/pulseR", subdir = "pkg")

if (length(args)>0 & args[1]=="all") {
    cran_pkgs <- c() # add any extra packages here
    cran_pkgs <- cran_pkgs[!cran_pkgs %in% installed.packages()[,1]]
    # DESeq2 fails to install in the newly created conda environment...
    # we temporarily rely on locally installed packages ()
    bioc_pkgs <- c(
        #"DESeq2",
        #"IHW",
        #"Glimma",
        #"topGO",
        #"org.Hs.eg.db",
        #"biomaRt"
    )
    bioc_pkgs <- bioc_pkgs[!bioc_pkgs %in% installed.packages()[,1]]
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    devtools::install_github("dieterich-lab/CellPlot", build_vignettes = FALSE)
} else {
  cran_pkgs <- c()
  bioc_pkgs <- c()
}

if (length(cran_pkgs) > 0)
    install.packages(cran_pkgs)
    
if (length(bioc_pkgs) > 0)
    BiocManager::install(bioc_pkgs)
