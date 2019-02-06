#' @include AllClasses.R
#' @include AllGenerics.R
#' @include methods-codonFreq.R
NULL



##------------------------------------------------------------------------------
## Mean codon usage frequency difference (MCUFD)
##------------------------------------------------------------------------------
#' Calculate mean codon usage similarity.
#'
#' Compare codon usage between sequences in a \code{codonFreq} object with those
#'    in a reference database. Mean codon usage frequency difference
#'    (\code{MCUFD}) is calculated as described in
#'    \href{https://goo.gl/7sGZeE}{Stedman et al. (2013)}.
#'
#' @param cFobj An object of class \code{codonFreq}.
#' @param exclude A character vector of codons to be excluded from comparisons.
#' @param minlen Numeric, the minimum length of sequence (in codons) to be
#'    included in the analysis. Default = 500.
#' @param norm Logical, should the data be normalised? If TRUE, a normalisation
#'    step will be performed. Default = FALSE.
#'
#' @return A list of data.frames containing values which indicate the degree of
#'    codon usage similarity between the sequences used to construct the
#'    \code{codonFreq} object and those used to construct the reference
#'    database. The reference database entries are in rows and the
#'    \code{codonFreq} sequences are in columns.
#'
#' @examples
#' ## To be added
#'
#' @name MCUFD
#' @rdname MCUFD
#'
#' @export
setMethod("MCUFD",
    signature(cFobj = "codonFreq"),
    function(cFobj, exclude, minlen, norm) {
        if (isTRUE(norm)) cFobj <- normalise(cFobj)
        seqidx <- 1:nrow(cFobj@freq)
        if (length(exclude) > 0) {
            keepcF <- which(!(colnames(cFobj@freq) %in% exclude))
            cFobj <- cFobj[seqidx, keepcF]
            keepRef <- which(!(colnames(gbnorm) %in% exclude))
            gbnorm <- gbnorm[, keepRef]
        }
        if (minlen > 0) {
            codidx <- 1:ncol(cFobj@freq)
            keepSeq <- which(cFobj@ncod > minlen)
            cFobj <- cFobj[keepSeq, codidx]
            gbnorm <- subset(gbnorm, X..Codons >= minlen)
        }
        resList <- vector("list", length = nseq(cFobj))
        names(resList) <- cFobj@seqID
        sapply(1:nrow(cFobj@freq), function(i) {
            cFrow <- cFobj@freq[i,]
            gbnorm[, 7:ncol(gbnorm)] <- abs(
                sweep(gbnorm[, 7:ncol(gbnorm)], 2, cFrow)
            )
            gbnorm$MCUFD <- rowMeans(gbnorm[, 7:ncol(gbnorm)], na.rm = TRUE)
            gbnorm$seqid <- cFobj@seqID[i]
            gbnorm <- gbnorm[order(gbnorm$MCUFD),]
            rownames(gbnorm) <- NULL
            resList[[i]] <<- gbnorm
        })
        resList
})



##------------------------------------------------------------------------------
## Plot functions
##------------------------------------------------------------------------------
#' MCUFD plot
#'
#' Plot the distribution of MCUFD values per taxon
#'
#' @param cFres List of data frames containing MCUFD results.
#' @param type Character, the type of plot to make. Options are "bar" (default)
#'    and "line".
#' @param n Numeric, the number of top-ranked reference taxa to be plotted
#'    per input sequence. Default = 100.
#' @param rank Character, taxonomic rank to be used for categorisation of the
#'    CUFD hit sequences. Options are "Domain", "Kingdom", and "Phylum".
#'    Default = "Phylum".
#' @param fname Character, name of figure generated.
#' @param units Numeric, units to be used for defining the plot size.
#'    Options are "in" (default), "cm", and "mm".
#' @param width Numeric, width of the figure (in \code{units}).
#' @param height Numeric, height of the figure (in \code{units}).
#' @param dpi Numeric, resolution of the figure (default = 600).
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#'    MCUFD_plot(
#'        MCUFD_tmp, type = "bar", fname = "MCUFD_bar_kingdom", n = 100,
#'        rank = "Kingdom"
#'    )
#'
#' @rdname plots
#'
#' @export
setMethod("MCUFD_plot",
    signature(cFres = "list"),
    function(cFres, type, n, rank, fname, units, width, height, dpi) {
        keepCols <- switch(rank,
            "Domain" = c(4, 68, 69),
            "Kingdom" = c(5, 68, 69),
            "Phylum" = c(6, 68, 69)
        )
        if (!is.na(n)) {
            cFres <- lapply(cFres, "[", 1:n, keepCols, drop = FALSE)
        } else {
            cFres <- lapply(cFres, "[", , keepCols, drop = FALSE)
        }
        cFres <- dplyr::bind_rows(cFres)
        cFres$seqid <- as.factor(str_sub(cFres$seqid, 1, 20))
        nseq <- nlevels(cFres$seqid)
        ## get the five most common taxon ranks:
        keepTax <- sort(
            as.data.frame(
                cFres %>% count_(rank, sort = TRUE) %>% top_n(5, n)
            )[,1]
        )
        cFtop <- subset(cFres, cFres[,1] %in% keepTax)
        cFres$retax <- ifelse(
            cFres[, rank] %in% keepTax, cFres[, rank], "Other"
        )
        cFres$retax <- factor(
            cFres$retax, levels = c(keepTax, "Other")
        )
        cc1 <- 12
        brewer_pallette <- brewer.pal(8, "Set1")
        if (type == "bar") {
            plt <- ggplot(cFres, aes_string(x = "seqid", fill = "retax")) +
                geom_bar(position = "fill") +
                theme_classic() +
                scale_fill_manual(
                    values = c(
                        brewer_pallette[1], brewer_pallette[5],
                        brewer_pallette[2], brewer_pallette[3],
                        brewer_pallette[4], "grey"
                    ),
                    name = rank,
                    breaks = levels(cFres$retax),
                    labels = levels(cFres$retax)
                ) +
                theme(
                    text = element_text(size = cc1),
                    axis.text.x = element_text(
                        angle = 45, hjust = 1, margin = margin(t = 5)
                    ),
                    plot.margin = unit(c(1,1,1,1), "cm")
                ) +
                xlab("Sequence") +
                ylab(paste0("Proportion of top ", n, " taxa"))+
                scale_x_discrete(expand = c(0.03, 0)) +
                scale_y_continuous(expand = c(0.02, 0))
        } else if (type == "histogram") {
            breaks <- pretty(
                range(cFres$MCUFD), n = (nclass.FD(cFres$MCUFD))*5, min.n = 1
            )
            bwidth <- breaks[2] - breaks[1]
            plt <- ggplot(cFres,
                aes_string(
                    x = "MCUFD", fill = rank, group = rank
                )) +
                geom_histogram(binwidth = bwidth) +
                theme_classic() +
                theme(
                    text = element_text(size = cc1)
                ) +
                xlab("MCUFD") + ylab("Density")
        } else {
            stop('Plot type must be either "bar" or "histogram".')
        }
        ggsave(
            plt,
            file = paste0(fname, ".png"),
            device = "png",
            units = units,
            width = width,
            height = height,
            dpi = dpi
        )
        plt
    }
)




##------------------------------------------------------------------------------
## Enrichment testing
##------------------------------------------------------------------------------
.enrichTest <- function(df, n, rank) {
    df[, rank] <- as.character(df[, rank])
    fullCounts <- data.frame(df %>% group_by_(rank) %>% tally())
    subCounts <- data.frame(df[1:n,] %>% group_by_(rank) %>% tally())
    fullTax <- fullCounts[, rank]
    ntax <- length(na.omit(fullTax))
    resDF <- data.frame(
        Taxon = character(ntax),
        foldEnrich = numeric(ntax),
        pvalue = numeric(ntax),
        padj = numeric(ntax),
        stringsAsFactors = FALSE
    )
    sapply(1:length(na.omit(fullTax)), function(x) {
        taxon <- na.omit(fullTax)[x]
        taxCountFull <- fullCounts[,2][grepl(taxon, fullCounts[,1])]
        taxCountSub <- subCounts[,2][grepl(taxon, subCounts[,1])]
        if (length(taxCountSub) == 0) taxCountSub <- 0
        otherCountFull <- sum(fullCounts[,2]) - taxCountFull
        otherCountSub <- sum(subCounts[,2]) - taxCountSub
        ratioFull <- taxCountFull/(otherCountFull+taxCountFull)
        ratioSub <- taxCountSub/(otherCountSub+taxCountSub)
        foldEnrich <- ratioSub - ratioFull
        contTab <- matrix(
            c(taxCountFull, taxCountSub, otherCountFull, otherCountSub),
            nrow = 2, ncol = 2
        )
        fisherp <- fisher.test(contTab)$p.value
        resDF[x, 1] <<- taxon
        resDF[x, 2] <<- foldEnrich
        resDF[x, 3] <<- fisherp
    })
    resDF[, 4] <- p.adjust(resDF[, 3], method = "BH")
    resDF <- resDF[order(resDF$padj),]
    resDF
}


#' MCUFD enrichment testing
#'
#' Compare taxon proportions in the top n hits with those in the full set of
#'    taxa.
#'
#' @param cFres List of data frames containing MCUFD results.
#' @param n Numeric, the number of top-ranked reference taxa to be plotted
#'    per input sequence. Default = 100.
#' @param rank Character, taxonomic rank to be used for categorisation of the
#'    CUFD hit sequences. Options are "Domain", "Kingdom", and "Phylum".
#'    Default = "Phylum".
#' @param plot Logical, should the enrichment results be plotted?
#'    Default = FALSE.
#' @param pthresh Numeric, adjusted p-value threshold to be used for subsetting
#'    the data if plot == TRUE.
#' @param fname Character, name of the output file (used if plot == TRUE).
#' @param units Numeric, units to be used for defining the plot size.
#'    Options are "in" (default), "cm", and "mm".
#' @param width Numeric, width of the figure (in \code{units}).
#' @param height Numeric, height of the figure (in \code{units}).
#' @param dpi Numeric, resolution of the figure (default = 600).
#'
#' @return A list of data frames containing enrichment results.
#'
#' @examples
#'    enrich_tmp <- MCUFD_enrich(
#'        MCUFD_tmp, n = 100, rank = "Phylum", plot = TRUE, pthresh = 0.05,
#'        fname = "enrich_tmp2", height = 5, width = 7
#'    )
#'
#' @export
setMethod("MCUFD_enrich",
    signature(cFres = "list"),
    function(cFres, n, rank, plot, pthresh, fname, units, width, height, dpi) {
        resList <- vector("list", length(cFres))
        sapply(seq_along(cFres), function(i) {
            resList[[i]] <<- .enrichTest(cFres[[i]], n, rank)
        })
        if (isTRUE(plot)) {
            resMerge <- dplyr::bind_rows(resList)
            if (!is.na(pthresh)) {
                resMerge <- subset(resMerge, padj < pthresh)
            }
            resMelt <- melt(
                resMerge, id.vars = c("Taxon"), measure.vars = "foldEnrich"
            )
            resMelt$Taxon <- factor(
                resMelt$Taxon, levels = rev(levels(factor(resMelt$Taxon)))
            )
            cc1 <- 12
            plt <- ggplot(
                    resMelt, aes(x = Taxon, y = value)
                ) +
                geom_point(
                    shape = 21, size = 4, fill = "red", alpha = 0.5,
                    show.legend = FALSE
                ) +
                theme_classic() +
                theme(text = element_text(size = cc1)) +
                coord_flip() +
                geom_hline(
                    yintercept = 0, linetype = "dashed",
                    color = "grey", size = 1
                ) +
                xlab("Taxon") +
                ylab(paste0("Fold enrichment in top ", n, " taxa"))
            ggsave(
                plt,
                file = paste0(fname, ".png"),
                device = "png",
                units = units,
                width = width,
                height = height,
                dpi = dpi
            )
        }
        resList
    }
)
