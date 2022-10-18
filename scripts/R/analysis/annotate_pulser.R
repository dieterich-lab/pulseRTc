#! /usr/bin/env Rscript

# Annotate pulseR estimates (all time points) and
# split estimates between low, intermediate and high decay rates
# based on the estimated rates for GAPDH (low), PDLIM5 (intermediate), MYC (high).

# Hard coded: uses org.Hs.eg.db and ensembl identifiers, etc.

# Usage: ./annotate_pulser.R [RESLOC_PULSER] [SRC] <0/1>
# 1: [RESLOC_PULSER] Results directory (pulseR)
# 2. [SRC] pulseR script source directory
# 3 <0/1> 0: only annotate (default), 1: split estimates 

# WARNING Hard coded: update path for missing libraries...
.libPaths( c( .libPaths(), "/beegfs/homes/eboileau/R/x86_64-pc-linux-gnu-library/4.2", "/beegfs/biosw/R/4.2.1_deb11/lib/R/library") )

library(dplyr)
library(tibble)
library(purrr)

library(org.Hs.eg.db)
library(biomaRt)


args <- commandArgs(trailingOnly=TRUE) 
split <- FALSE
if (length(args)<2) {
  stop("./annotate_pulser.R [RESLOC_PULSER] [SRC] <0/1>\n", call.=FALSE)
} else {
  if (length(args)>2 & as.integer(args[3])==1) {
    split <- TRUE
  }
}

src <- args[2]
source(file.path(src, "time_pts.R", fsep=.Platform$file.sep)) 
timeSets <- lapply(allSets, function(ts) { paste(ts, collapse="-") })

# pulseR results 
pulseDir <- file.path(args[1], "featureCounts", fsep=.Platform$file.sep)
pulseFit <- "pulsefit"
pulseCis <- "pulsecis"
# combined output
tabDir <- file.path(pulseDir, "analysis", "tables", fsep=.Platform$file.sep)

# genes
# actually we use only low and high...
MYC <- "ENSG00000136997"
PDLIM5 <- "ENSG00000163110"
GAPDH <- "ENSG00000111640"

# data 
fit <- readRDS(file.path(pulseDir, paste(pulseFit, "-", timeSets, ".rds", sep=""), fsep=.Platform$file.sep))
rates <- fit$fit$d
names(rates) <- rownames(fit$pd$counts)

mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host="apr2019.archive.ensembl.org")
resMArt <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol", "gene_biotype", "description"), mart=mart)
    
if (isTRUE(split)) {
    high <- rates[names(rates)==MYC]
    int <- rates[names(rates)==PDLIM5]
    low <- rates[names(rates)==GAPDH]

    highr <- rates[rates>=high]
    intr <- rates[(rates>=low)&(rates<high)]
    lowr <- rates[rates<low]

    lowr <- as.data.frame(lowr)
    lowr$ensembl_gene_id <- rownames(lowr)
    colnames(lowr) <- c("rate", 'ensembl_gene_id')
    lowr$class <- 'low'
    lowr <- merge(lowr, resMArt, by='ensembl_gene_id', all.x=T)

    intr <- as.data.frame(intr)
    intr$ensembl_gene_id <- rownames(intr)
    colnames(intr) <- c("rate", 'ensembl_gene_id')
    intr$class <- 'intermediate'
    intr <- merge(intr, resMArt, by='ensembl_gene_id', all.x=T)

    highr <- as.data.frame(highr)
    highr$ensembl_gene_id <- rownames(highr)
    colnames(highr) <- c("rate", 'ensembl_gene_id')
    highr$class <- 'high'
    highr <- merge(highr, resMArt, by='ensembl_gene_id', all.x=T)

    annotated <- do.call("rbind", list(lowr, intr, highr))
} else {
    rates <- as.data.frame(rates)
    rates$ensembl_gene_id <- rownames(rates)
    colnames(rates) <- c("rate", 'ensembl_gene_id')
    annotated <- merge(rates, resMArt, by='ensembl_gene_id', all.x=T)
}
  
library(openxlsx)

workBook <- createWorkbook()
addWorksheet(workBook, sheetName='annotated')
writeDataTable(workBook, sheet=1, x=annotated, rowNames = FALSE)
saveWorkbook(workBook, file.path(tabDir, paste(pulseFit, "-", timeSets, "-annotated.xlsx", sep="")), overwrite=TRUE)
