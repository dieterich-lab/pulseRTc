#! /usr/bin/env Rscript

##' Prepare count tables and metadata for analysis from Salmon counts
##' Note: this function is not all purpose, we know the exact format of the data!
##' No filtering based on expression, etc.
##' ** Uses ENSG identifiers

# Usage: ./prep_salmon.R [LOC] [OUTPUT] [NAME] <DATA>
# 1: [LOC] Input directory (count table)
# 2. [OUTPUT] Output directory - created if does not exists
# 3. [NAME] Name of input file
# 4. <DATA> std: biochemical separation (anything else defaults to output from splbam).
#           If given, then name will be used to construct input file name (default fc-).

library(purrr)
library(tibble)
library(dplyr)

# get params and options
args <- commandArgs(trailingOnly=TRUE)
data <- "fc"
if (length(args)<3) {
  stop("USAGE: plot_mismatches.R [LOC] [OUTPUT] [NAME] <DATA>\n", call.=FALSE)
} else {
  if (length(args)>3) {
    data <- args[4]
  }
}

getData <- function(data, dataDir, dataFile, outDir){
    loc <- file.path(dataDir, "tables", "featureCounts", fsep=.Platform$file.sep)
    cts <- read.table(file.path(loc, dataFile, fsep=.Platform$file.sep), row.names=1, sep = "\t", header=TRUE, check.names=FALSE)
    cts <- cts %>% dplyr::select(starts_with(file.path(dataDir, "mapping", fsep=.Platform$file.sep)))
    depth <- length(strsplit(colnames(cts)[1], "/", fixed=T)[[1]])
    conditions <- data.frame(vapply(strsplit(colnames(cts), "/", fixed=T), "[", "", depth), stringsAsFactors=FALSE)
    colnames(conditions) <- "sample"
    if (data == 'std') {
        conditions$fraction <- vapply(strsplit(conditions$sample, "_", fixed=T), "[", "", 2)
        conditions$fraction[grepl("flow", conditions$fraction)] <- "flow_through"
        conditions$fraction[grepl("enriched", conditions$fraction)] <- "pull_down"
        conditions$time <- as.numeric(gsub("h", "", gsub("^.*_", "", vapply(strsplit(conditions$sample, ".", fixed=T), "[", "", 1))))
        conditions$rep <- gsub("^[1-9][0-9][0-9][0-9][0-9][0-9]", "", vapply(strsplit(conditions$sample, "_", fixed=T), "[", "", 1))
        conditions$sample <- gsub("[A-Z]$", "", vapply(strsplit(conditions$sample, "_", fixed=T), "[", "", 1))
        # redefine columns
        colnames(cts) <- conditions$sample
        # remove spike-in mix, but keep ERCC
        cts <- cts[grepl("^ENSG|^ERCC", rownames(cts)),]
        # if we use fractional counts (multimappers)
        cts <- round(cts)
    } else { # slam, tls, or tuc
        conditions$fraction <- vapply(strsplit(vapply(strsplit(conditions$sample, "_", fixed=T), "[", "", 3), ".", fixed=T), "[", "", 2)
        conditions$time <- as.numeric(gsub("h", "", vapply(strsplit(vapply(strsplit(conditions$sample, "_", fixed=T), "[", "", 3), ".", fixed=T), "[", "", 1)))
        conditions$rep <- gsub("WT", "", vapply(strsplit(conditions$sample, "_", fixed=T), "[", "", 2))
        conditions$sample <- gsub("[A-Z]$", "", vapply(strsplit(conditions$sample, "_", fixed=T), "[", "", 1))
        conditions$sample <- paste(conditions$sample, "L", sep="")
        conditions$sample[grepl("unlabelled", conditions$fraction)] <- gsub("L", "U", conditions$sample[grepl("unlabelled", conditions$fraction)])
        # redefine columns
        colnames(cts) <- conditions$sample
        # remove spike-in mix and ERCC (unused)
        cts <- cts[grepl("^ENSG", rownames(cts)),]
        # if we use fractional counts (multimappers)
        cts <- round(cts)
    }
    saveRDS(cts, file.path(outDir, "counts.rds", fsep=.Platform$file.sep))
    saveRDS(conditions, file.path(outDir, "samples.rds", fsep=.Platform$file.sep))
}


## Call

outDir <- file.path(args[2], "featureCounts", "data", fsep=.Platform$file.sep)
if (!dir.exists(outDir)) {dir.create(outDir, recursive=TRUE)}

tsv <- paste(data, args[3], sep="")
print(paste("Processing ", file.path(args[1], "tables", tsv, fsep=.Platform$file.sep), " ...", sep=""))
print(paste("Writing to ", outDir, " ...", sep=""))

getData(data, args[1], tsv, outDir)
