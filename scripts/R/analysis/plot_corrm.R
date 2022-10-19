#! /usr/bin/env Rscript

# Compare pulseR (featureCounts) with GRAND-SLAM

# Usage: ./plot_corrm.R [RESLOC_PULSER] [SRC]
# 1: [RESLOC_PULSER] Results directory (pulseR)
# 2. [SRC] pulseR script source directory


library(dplyr)
library(tibble)
library(purrr)

library(ggplot2)
library(RColorBrewer)

library(cowplot)
library(grid)
library(gridExtra)


## input/params

# pulseR, GS c('#3182BD', '#31A354')
pal <- c(rev(scales::brewer_pal(pal="Blues")(6))[2],  
         rev(scales::brewer_pal(pal="Greens")(6))[2])
         

args <- commandArgs(trailingOnly=TRUE)
if (length(args)<2) { stop("./plot_corrm.R [RESLOC_PULSER] [SRC]\n", call.=FALSE) }

src <- args[2]
source(file.path(src, "utils.R", fsep=.Platform$file.sep)) 

# pulseR/GRAND-SLAM combined results 
pulseDir <- file.path(args[1], "featureCounts", fsep=.Platform$file.sep)
pulseFit <- "pulsefit"
pulseCis <- "pulsecis"
# combined output
tabDir <- file.path(pulseDir, "analysis", "tables", fsep=.Platform$file.sep)
# combined output figures - created
figDir <- file.path(pulseDir, "analysis", "figures", fsep=.Platform$file.sep)
if (!dir.exists(figDir)) {dir.create(figDir, recursive=TRUE)}


## plotting

shadesOfGrey <- colorRampPalette(c("white", "grey100", "grey90", "grey80", "grey70", "grey60", 
                                   "grey50", "grey40", "grey30", "grey20", "grey10", "grey0"))
                                   
                                   
myHeatmap <- function(cormat, labels, llim=0.6) {
    # melt fo plotting...
	melted_cormat <- reshape2::melt(cormat, na.rm = TRUE)
	# manually set column order 
    col_order <- rev(levels(melted_cormat$Var2))
    melted_cormat$Var1 <- factor(melted_cormat$Var1, levels = col_order)
    melted_cormat$Var2 <- factor(melted_cormat$Var2, levels = col_order)
	# and create a ggheatmap
	ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value)) +
        geom_tile(color = "white") +
        scale_fill_gradient(low = "#ece7f2", high = "#3182BD", limit = c(llim,1), space = "Lab", 
            name="Pearson\nCorrelation",na.value='white',guide=F) +
        theme_minimal() + 
        theme(text=element_text(size=12),
              axis.text.y = element_text(size=12),
              axis.text.x = element_text(angle = 90, vjust = 0.5, size = 12, hjust = 1)) +
		scale_x_discrete(labels=labels) + scale_y_discrete(labels=labels) +
        coord_fixed()
	# add text (corr values)
	ggheatmap + 
        geom_text(aes(Var2, Var1, label = value), color = "black", size = 2) +
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              panel.grid.major = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
	          axis.ticks = element_blank())
}                                   
                                   
                                   
myLine <- function(x, y, ...){
    par(new = TRUE)
    smoothScatter(x, y, 
                  colramp = shadesOfGrey,
                  xaxt='n', yaxt = 'n',
                  ...)
    abline(a = 0,b = 1, ...)
}


# Viewport function
vplayout <- function(x,y){
  viewport(layout.pos.row = x, layout.pos.col = y)
}


log10.axis <- function(side, at, ...) {
    at.minor <- log10(outer(1:9, 10^(min(at):max(at))))
    axis(side=side, at=at.minor, labels=NA, tcl=par("tcl")*0.5, ...)
}


lower.panel <- function(x, y, i, j, ...){
    # if comparing to ERCC, use residual SE from linear fit
    # ERCC always need to be on first column!
    # else use RMSD
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    if (j==1) {
        f <- lm(y ~ x)
        rse <- round(sqrt(deviance(f)/f$df.residual), 2)
        rse <- format(rse, digits=2, nsmall=2)
        txt <- paste0("RSE\n", rse)
    } else {
        rmsd <- round(sqrt(mean((y - x)^2)), 2)
        rmsd <- format(rmsd, digits=2, nsmall=2)
        txt <- paste0("RMSD\n", rmsd)
    }
    cex.cor <- .8/strwidth(txt)
    text(grconvertX(0.5,"npc"), grconvertY(0.5, "npc"), txt, cex = cex.cor)
}


# passing indices to lower panels
# adapted from R function pairs
myPairs <-

function (x, labels, panel = points, ...,

          horInd = 1:nc, verInd = 1:nc,

          lower.panel = panel, upper.panel = panel,

          diag.panel = NULL, text.panel = textPanel,

          label.pos = 0.5 + has.diag/3, line.main = 3,

          cex.labels = NULL, font.labels = 1,

          row1attop = TRUE, gap = 1, log = "",

          horOdd = !row1attop, verOdd = !row1attop)

{

    if(doText <- missing(text.panel) || is.function(text.panel))

	textPanel <-

	    function(x = 0.5, y = 0.5, txt, cex, font)

		text(x, y, txt, cex = cex, font = font)

 

    localAxis <- function(side, x, y, xpd, bg, col=NULL, main, oma, ...) {

      ## Explicitly ignore any color argument passed in as

      ## it was most likely meant for the data points and

      ## not for the axis.

        xpd <- NA

        if(side %% 2L == 1L && xl[j]) xpd <- FALSE

        if(side %% 2L == 0L && yl[i]) xpd <- FALSE

        if(side %% 2L == 1L) Axis(x, side = side, xpd = xpd, ...)

        else Axis(y, side = side, xpd = xpd, ...)

    }

 

    localPlot <- function(..., main, oma, font.main, cex.main) plot(...)

    localLowerPanel <- function(..., main, oma, font.main, cex.main)

        lower.panel(...)

    localUpperPanel <- function(..., main, oma, font.main, cex.main)

        upper.panel(...)

 

    localDiagPanel <- function(..., main, oma, font.main, cex.main)

        diag.panel(...)

 

    dots <- list(...); nmdots <- names(dots)

    if (!is.matrix(x)) {

        x <- as.data.frame(x)

        for(i in seq_along(names(x))) {

            if(is.factor(x[[i]]) || is.logical(x[[i]]))

               x[[i]] <- as.numeric(x[[i]])

            if(!is.numeric(unclass(x[[i]])))

                stop("non-numeric argument to 'pairs'")

        }

    } else if (!is.numeric(x)) stop("non-numeric argument to 'pairs'")

    panel <- match.fun(panel)

    if((has.lower <- !is.null(lower.panel)) && !missing(lower.panel))

        lower.panel <- match.fun(lower.panel)

    if((has.upper <- !is.null(upper.panel)) && !missing(upper.panel))

        upper.panel <- match.fun(upper.panel)

    if((has.diag  <- !is.null( diag.panel)) && !missing( diag.panel))

        diag.panel <- match.fun( diag.panel)

 

    if(row1attop) {

        tmp <- lower.panel; lower.panel <- upper.panel; upper.panel <- tmp

        tmp <- has.lower; has.lower <- has.upper; has.upper <- tmp

    }

 

    nc <- ncol(x)

    if (nc < 2L) stop("only one column in the argument to 'pairs'")

    if(!all(1L <= horInd & horInd <= nc))

        stop("invalid argument 'horInd'")

    if(!all(1L <= verInd & verInd <= nc))

        stop("invalid argument 'verInd'")

    if(doText) {

	if (missing(labels)) {

	    labels <- colnames(x)

	    if (is.null(labels)) labels <- paste("var", 1L:nc)

	}

	else if(is.null(labels)) doText <- FALSE

    }

    oma  <- if("oma"  %in% nmdots) dots$oma

    main <- if("main" %in% nmdots) dots$main

    if (is.null(oma))

	oma <- c(4, 4, if(!is.null(main)) 6 else 4, 4)

    opar <- par(mfcol = c(length(horInd), length(verInd)),

                mar = rep.int(gap/2, 4), oma = oma)

    on.exit(par(opar))

    dev.hold(); on.exit(dev.flush(), add = TRUE)

 

    xl <- yl <- logical(nc)

    if (is.numeric(log)) xl[log] <- yl[log] <- TRUE

    else {xl[] <- grepl("x", log); yl[] <- grepl("y", log)}

    ni <- length(iSet <- if(row1attop) horInd else rev(horInd))

    nj <- length(jSet <- verInd)

    for(j in jSet)

        for(i in iSet) {

            l <- paste0(if(xl[j]) "x" else "",

                        if(yl[i]) "y" else "")

            localPlot(x[, j], x[, i], xlab = "", ylab = "",

                      axes = FALSE, type = "n", ..., log = l)

            if(i == j || (i < j && has.lower) || (i > j && has.upper) ) {

                box()

                j.odd <- (match(j, jSet) + horOdd) %% 2L

                i.odd <- (match(i, iSet) + verOdd) %% 2L

                if(i == iSet[1L] && (!j.odd || !has.upper || !has.lower))

                    localAxis(3L, x[, j], x[, i], ...)

                if(i == iSet[ni] && ( j.odd || !has.upper || !has.lower))

                    localAxis(1L, x[, j], x[, i], ...)

                if(j == jSet[1L] && (!i.odd || !has.upper || !has.lower))

                    localAxis(2L, x[, j], x[, i], ...)

                if(j == jSet[nj] && ( i.odd || !has.upper || !has.lower))

                    localAxis(4L, x[, j], x[, i], ...)

                mfg <- par("mfg")

                if(i == j) {

                    if (has.diag) localDiagPanel(as.vector(x[, i]), ...)

		    if (doText) {

                        par(usr = c(0, 1, 0, 1))

                        if(is.null(cex.labels)) {

                            l.wid <- strwidth(labels, "user")

                            cex.labels <- max(0.8, min(2, .9 / max(l.wid)))

                        }

                        xlp <- if(xl[i]) 10^0.5 else 0.5

                        ylp <- if(yl[j]) 10^label.pos else label.pos

                        text.panel(xlp, ylp, labels[i],

                                   cex = cex.labels, font = font.labels)

                    }

                } else if(i < j)

                    localLowerPanel(as.vector(x[, j]), as.vector(x[, i]), ...)

                else

                    localUpperPanel(as.vector(x[, j]), as.vector(x[, i]), i, j, ...)

                if (any(par("mfg") != mfg))

                    stop("the 'panel' function made a new plot")

            }

            else par(new = FALSE)

        }

    if (!is.null(main)) {

        font.main <- if("font.main" %in% nmdots) dots$font.main else par("font.main")

        cex.main  <- if("cex.main"  %in% nmdots) dots$cex.main  else par("cex.main")

        mtext(main, 3, line.main, outer=TRUE, at = 0.5, cex = cex.main, font = font.main)

    }

    invisible(NULL)

}

corrLabels <- function(n) {
    labels <- c()
    for (c in gsub("\n", " ", n)) {
        components <- unlist(strsplit(c, " "))
        labels <- c(labels, bquote(.(components[1])[.(components[2])]))
    }
    names(labels) <- n
    labels
}


corrLabelsFull <- function(n) {
    labels <- c()
    for (c in gsub("\n", " ", n)) {
        components <- unlist(strsplit(c, " "))
        labels <- c(labels, bquote(.(components[1])[.(components[2])]^.(paste(components[-1:-2], collapse=","))))
    }
    names(labels) <- n
    labels
}


plotCorrd <- function(r, l, n, m=0.6, hjust=1.5, vjust=1.5) {
    title <- substitute(paste("Decay rate ", delta, " ", h^-1, " (n=",n, ")"), list(n = n))
    myHeatmap(round(cor(as.matrix(r), use='p'), 2), l, llim = m) + ggtitle(title) +
    theme(plot.title = element_text(hjust = hjust, vjust = vjust))
}


## call

res <- map(dir(tabDir, "tbl-0"), ~readRDS(file.path(tabDir, .x)))
timeSets <- extractTimesFromNames(dir(tabDir, "tbl-0"))
names(res) <- timeSets

# adjust all names
res <- map(res,
           function(df) {
           d <- df[, grep("d$", colnames(df))]
           n <- gsub(" d$", "", colnames(d))
           n <- gsub("pulseR ", "pulseR\n", n)
           n <- gsub("GS ", "GS\n", n)
           colnames(d) <- n
           d
})

# corr heatmap
supp_heatmaps <- map(res, function(r) {
        labels <- corrLabels(colnames(r))
        plotCorrd(r, labels, dim(r)[1])
    })
names(supp_heatmaps) <- names(res)

# pairs
# we need to write them to disk first, does not work as function...
map2(res, names(res), function(x, y) {
    filen <- file.path(figDir, paste("pair-", y, ".pdf", sep=""), fsep=.Platform$file.sep)
    pdf(filen, width=16, height=16)
    myPairs(x, 
        lower.panel = lower.panel, 
        upper.panel = myLine,
        cex.labels = 4, cex.axis = 2.5,
        main = "")
    dev.off()    
})

map2(supp_heatmaps, names(supp_heatmaps), function(x, y) {
    filen <- file.path(figDir, paste("pair-", y, ".pdf", sep=""), fsep=.Platform$file.sep)
    p2 <- ggdraw() + draw_image(filen)
    g <- plot_grid(x, p2, labels = "auto", ncol = 2, nrow = 1, align = 'h')
    save_plot(file.path(figDir, paste("figure3_d-", y, ".pdf", sep=""), fsep=.Platform$file.sep), g)
})

# all together
names(res) <- NULL
all_res <- do.call(cbind, res)
labels <- corrLabelsFull(colnames(all_res))
heatmap <- plotCorrd(all_res, labels, dim(all_res)[1], m=0.5, hjust=0.5)
ggsave(file.path(figDir, paste("figure3_d-all.pdf", sep=""), fsep=.Platform$file.sep), plot=heatmap)
