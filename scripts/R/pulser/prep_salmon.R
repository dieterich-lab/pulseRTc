#! /usr/bin/env Rscript

##' Prepare count tables and metadata for analysis from Salmon counts
##' Note: this function is not all purpose, we know the exact format of the data!
##' No filtering based on expression, etc.
##' Currently dealing with SLAM-like data preparation only!

# Usage: ./prep_salmon.R [LOC] [OUTPUT] <DATA>
# 1: [LOC] Input directory (Salmon directories)
# 2. [OUTPUT] Output directory - created if does not exists

library(purrr)
library(tibble)
library(dplyr)

# get params and options
args <- commandArgs(trailingOnly=TRUE)
data <- "salmon"
if (length(args)<2) {
  stop("USAGE: plot_mismatches.R [LOC] [OUTPUT] <DATA>\n", call.=FALSE)
} else {
  if (length(args)>2) {
    data <- args[3]
  }
}

# NumReads
colClasses <- c("character", "NULL", "NULL", "NULL", "numeric")

getData <- function(data, dataDir, dataFiles, outDir){
    squants <- map(dataFiles, function(.id) {
        df <- read.table(file.path(dataDir, .id, "quant.sf", fsep=.Platform$file.sep), 
                         header = T, 
                         row.names = 1, 
                         colClasses = colClasses)
        colnames(df) <- .id
        df
    })
    cts <- bind_cols(squants)
    conditions <- data.frame(colnames(cts), stringsAsFactors=FALSE)
    colnames(conditions) <- "sample"
    if (data == 'std') {
        stop("OPTION <std> NOT IMPLEMENTED! TERMINATING\n", call.=FALSE)
    } else { # slam, tls, or tuc
        conditions$fraction <- vapply(strsplit(vapply(strsplit(conditions$sample, "_", fixed=T), "[", "", 3), ".", fixed=T), "[", "", 2)
        conditions$time <- as.numeric(gsub("h", "", vapply(strsplit(vapply(strsplit(conditions$sample, "_", fixed=T), "[", "", 3), ".", fixed=T), "[", "", 1)))
        conditions$rep <- gsub("WT", "", vapply(strsplit(conditions$sample, "_", fixed=T), "[", "", 2))
        conditions$sample <- gsub("[A-Z]$", "", vapply(strsplit(conditions$sample, "_", fixed=T), "[", "", 1))
        conditions$sample <- paste(conditions$sample, "L", sep="")
        conditions$sample[grepl("unlabelled", conditions$fraction)] <- gsub("L", "U", conditions$sample[grepl("unlabelled", conditions$fraction)])
        # redefine columns
        colnames(cts) <- conditions$sample
        # if we use fractional counts (multimappers)
        cts <- round(cts)
    }
    saveRDS(cts, file.path(outDir, "counts.rds", fsep=.Platform$file.sep))
    saveRDS(conditions, file.path(outDir, "samples.rds", fsep=.Platform$file.sep))
}


## Call

outDir <- file.path(args[2], "Salmon", "data", fsep=.Platform$file.sep)
if (!dir.exists(outDir)) {dir.create(outDir, recursive=TRUE)}

loc <- file.path(args[1], "tables", "Salmon", fsep=.Platform$file.sep)
sfiles <- list.files(loc)
print(paste("Processing ", sfiles, " ...", sep=""))
print(paste("Writing to ", outDir, " ...", sep=""))

getData(data, loc, sfiles, outDir)
