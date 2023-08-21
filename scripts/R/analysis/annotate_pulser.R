#! /usr/bin/env Rscript

# Annotate pulseR estimates (all time points) and
# split estimates between low, intermediate and high decay rates
# based on the estimated rates for GAPDH (low), PDLIM5 (intermediate), MYC (high).

# Hard coded: uses org.Hs.eg.db and ensembl identifiers, etc.

# Usage: ./annotate_pulser.R [RESLOC_PULSER] [SRC] [0/1] <0/1>
# 1: [RESLOC_PULSER] Results directory (pulseR)
# 2: [SRC] pulseR script source directory
# 3: [0/1] 0: gene, 1: transcript annotation
# 4: <0/1> 0: only annotate (default), 1: split estimates

# WARNING Hard coded: update path for missing libraries...
.libPaths( c( .libPaths(), "/beegfs/homes/eboileau/R/x86_64-pc-linux-gnu-library/4.2", "/beegfs/biosw/R/4.2.1_deb11/lib/R/library") )

# WARNING Hard coded annotation org.Hs.eg.db, etc.

library(dplyr)
library(tibble)
library(purrr)

library(org.Hs.eg.db)
library(biomaRt)


args <- commandArgs(trailingOnly=TRUE)
split <- FALSE
if (length(args)<3) {
  stop("./annotate_pulser.R [RESLOC_PULSER] [SRC] [0/1] <0/1>\n", call.=FALSE)
} else {
  if (length(args)>3 & as.integer(args[4])==1) {
    split <- TRUE
  }
}

src <- args[2]
source(file.path(src, "time_pts.R", fsep=.Platform$file.sep))
timeSets <- lapply(allSets, function(ts) { paste(ts, collapse="-") })

# Ensembl 102
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host="nov2020.archive.ensembl.org")
# resMArt.prot <- getBM(attributes=c("ensembl_gene_id", "ensembl_transcript_id", "uniprotswissprot", "uniprotsptrembl"), mart=mart)
if (args[3] == 0) {
    loc <- "featureCounts"
    key <- "ensembl_gene_id"
    resMArt <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol", "gene_biotype", "description"), mart=mart)
} else {
    loc <- "Salmon"
    key <- "ensembl_transcript_id"
    resMArt <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol", "ensembl_transcript_id", "transcript_biotype"), mart=mart)
}

# pulseR results
pulseDir <- file.path(args[1], loc, fsep=.Platform$file.sep)
pulseFit <- "pulsefit"
# combined output
tabDir <- file.path(pulseDir, "analysis", "tables", fsep=.Platform$file.sep)
if (!dir.exists(tabDir)) {dir.create(tabDir, recursive=TRUE)}

# genes
# actually we use only low and high...
MYC <- "ENSG00000136997"
PDLIM5 <- "ENSG00000163110"
GAPDH <- "ENSG00000111640"

# data
fit <- readRDS(file.path(pulseDir, paste(pulseFit, "-", timeSets, ".rds", sep=""), fsep=.Platform$file.sep))
rates <- fit$fit$d
names(rates) <- rownames(fit$pd$counts)

# split currently hard coded to deal with gene ids only...
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
    rates[[key]] <- rownames(rates)
    colnames(rates) <- c("rate", key)
    annotated <- merge(rates, resMArt, by=key, all.x=T)
}

# add half-lives
annotated$half_life <- log(2)/annotated$rate
annotated <- annotated[,c(1,3,4,5,6,2,7)]

library(openxlsx)

workBook <- createWorkbook()
addWorksheet(workBook, sheetName='annotated')
writeDataTable(workBook, sheet=1, x=annotated, rowNames = FALSE)
saveWorkbook(workBook, file.path(tabDir, paste(pulseFit, "-", timeSets, "-annotated.xlsx", sep="")), overwrite=TRUE)
