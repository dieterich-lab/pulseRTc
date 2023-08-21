#! /usr/bin/env Rscript

# Uses output of featureCounts on "un-split" BAM files
# to check gene expression at varying labeling times vs. 0 (unlabelled).

# Usage: ./run_deseq.R [LOC] [RESLOC] [NAME] [SRC]
# 1: [LOC] Location of BAM files
# 2: [RESLOC] Results directory (featureCounts for DESeq)
# 3. [NAME] Name of input file
# 4. [SRC] pulseR script source directory

# WARNING Hard coded: update path for missing libraries...
.libPaths( c( .libPaths(), "/beegfs/homes/eboileau/R/x86_64-pc-linux-gnu-library/4.2", "/beegfs/biosw/R/4.2.1_deb11/lib/R/library") )

# WARNING Hard coded annotation org.Hs.eg.db, etc. and FRD/LFC cut-off

library(DESeq2)
library(IHW)

library("Glimma")
library("genefilter")
library("org.Hs.eg.db")

library(dplyr)
library(purrr)
library(tibble)

library(openxlsx)


# get params and options
args <- commandArgs(trailingOnly=TRUE)
data <- "dge"
if (length(args)<4) { stop("USAGE: plot_mismatches.R [LOC] [RESLOC] [NAME] [SRC] \n", call.=FALSE) }
source(file.path(args[4], "utils.R", fsep=.Platform$file.sep))
source(file.path(args[4], "time_pts.R", fsep=.Platform$file.sep))

# ---------------------------------------------------------

## Default LFC and FDR threshold

lfcThreshold.set <- log2(1.2)
altHypothesis.set <- "greaterAbs"
alpha.set <- 0.05

# ---------------------------------------------------------

## functions

# filter out lowly expressed genes using geometric mean
whichHigh <- function(x, level) {
  apply(x, 1, function(y) exp(mean(log(y)))) > level
}

# DESeq results
get_results <- function (contrast, dds) {

    num <- as.character(contrast[1])
    denom <- as.character(contrast[2]) # reference

    print(paste("Contrast: time_", num, "_vs_", denom, sep=""))

    res <- results(dds,
                   contrast=c("time", num, denom),
                   lfcThreshold=lfcThreshold.set,
                   altHypothesis=altHypothesis.set,
                   alpha=alpha.set,
                   filterFun=ihw)
    res$padj[is.na(res$padj)] <- 1

    res$symbol <- mapIds(org.Hs.eg.db,
                         keys=rownames(res),
                         column="SYMBOL",
                         keytype="ENSEMBL",
                         multiVals="first")

    res.shrunken <- lfcShrink(dds,
                              coef=resultsNames(dds)[grepl(paste('time_', num, '_vs_0', sep=''), resultsNames(dds))],
                              res=res,
                              type="apeglm",
                              lfcThreshold=lfcThreshold.set)
    res.shrunken$pvalue <- res$pvalue
    res.shrunken$padj <- res$padj

    res.shrunken$symbol <- mapIds(org.Hs.eg.db,
                                  keys=rownames(res.shrunken),
                                  column="SYMBOL",
                                  keytype="ENSEMBL",
                                  multiVals="first")

    is.de <- as.numeric(res.shrunken$padj < alpha.set & abs(res.shrunken$log2FoldChange) > lfcThreshold.set)
    anno <- data.frame(GeneID=rownames(res.shrunken), symbol=res.shrunken$symbol)
    glMDPlot(res.shrunken,
             counts=counts(dds ,normalized=TRUE),
             anno,
             dds$time,
             samples=colnames(dds),
             status=is.de,
             transform = FALSE,
             xlab = "logMeanExpr",
             ylab = "log2FoldChange",
             side.ylab = "NormalizedCount",
             path=outDir,
             folder=paste("glimma-plots", num, "_vs_", denom, sep=""),
             launch=FALSE)

    res.tib <- res %>%
        data.frame() %>%
        rownames_to_column(var="gene") %>%
        as_tibble()

    res.shrunken.tib <- res.shrunken %>%
        data.frame() %>%
        rownames_to_column(var="gene") %>%
        as_tibble()

    # write to disk
    wb <- createWorkbook()

    addWorksheet(wb, sheetName=paste0(num, "_vs_", denom, sep=""))
    writeDataTable(wb, sheet=1, x=res.tib)

    addWorksheet(wb, sheetName=paste0(num, "_vs_", denom, "_shrunken", sep=""))
    writeDataTable(wb, sheet=2, x=res.shrunken.tib)

    filen <- paste0("time_", num, "_vs_", denom, ".xlsx", sep="")
    filen <- file.path(outDir, filen, fsep=.Platform$file.sep)
    saveWorkbook(wb, filen, overwrite=TRUE)

}


# ---------------------------------------------------------

## call

outDir <- file.path(args[2], "DESeq", "results", fsep=.Platform$file.sep)
if (!dir.exists(outDir)) {dir.create(outDir, recursive=TRUE)}

loc <- file.path(args[2], "DESeq", fsep=.Platform$file.sep)
tsv <- paste(data, args[3], sep="")

print(paste("Processing ", file.path(loc, "tables", tsv, fsep=.Platform$file.sep), " ...", sep=""))
print(paste("Writing to ", outDir, " ...", sep=""))

# read counts and prep coldata
cts <- read.table(file.path(loc, "tables", tsv, fsep=.Platform$file.sep), row.names=1, sep = "\t", header=TRUE, check.names=FALSE)
cts <- cts %>% dplyr::select(starts_with(file.path(args[1], fsep=.Platform$file.sep)))
depth <- length(strsplit(colnames(cts)[1], "/", fixed=T)[[1]])
conditions <- data.frame(vapply(strsplit(colnames(cts), "/", fixed=T), "[", "", depth), stringsAsFactors=FALSE)
colnames(conditions) <- "sample"
conditions$time <- as.numeric(gsub("h", "", vapply(strsplit(vapply(strsplit(conditions$sample, "_", fixed=T), "[", "", 3), ".", fixed=T), "[", "", 1)))
# here we're not using time as a continuous variable e.g. to identify genes that change as a function of time
# but rather we want to compare labelling time points in terms of expression signatures, so we use them as factors/contrasts
conditions$time <- as.factor(conditions$time)
conditions$rep <- gsub("WT", "", vapply(strsplit(conditions$sample, "_", fixed=T), "[", "", 2))
conditions$rep <- as.factor(conditions$rep)
conditions$sample <- gsub("[A-Z]$", "", vapply(strsplit(conditions$sample, "_", fixed=T), "[", "", 1))
rownames(conditions) <- conditions$sample
# redefine columns
colnames(cts) <- conditions$sample
# remove spike-in mix and ERCC (unused) if any
cts <- cts[grepl("^ENSG", rownames(cts)),]
# if we use fractional counts (multimappers)
cts <- round(cts)

# expression threshold - same as for pulseR
# note however that the `results` function performs independent filtering by default using the
# mean of normalized counts as a filter statistic
expressionThreshold <- 50
highExpr <- whichHigh(1 + cts, expressionThreshold)
cts <- cts[highExpr,]

stopifnot(all(rownames(conditions) == colnames(cts)))

# contrasts to test - based on all sets (all time points)
allSets <- as.character(unlist(allSets))
allSets <- allSets[!allSets==0]
contrasts <- data.frame(allSets, rep("0", length(allSets)))
colnames(contrasts) <- c('cond', 'ref')

# construct DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData=cts,
                              colData=conditions,
                              design=~time)

# fit all at once
dds$time <- relevel(dds$time, ref='0')
dds <- DESeq(dds)

apply(contrasts, 1, get_results, dds=dds)
