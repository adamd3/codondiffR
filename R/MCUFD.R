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
#'    virusSet <- readSeq(example = TRUE)
#'    virusCF <- codonFreq(virusSet)
#'    exclCod <- c("ATT", "TGT")
#'    MCUFD_tmp <- MCUFD(virusCF, exclude = exclCod, norm = TRUE)
#'    range(MCUFD_tmp[[1]]$MCUFD)
#'    table(MCUFD_tmp[[1]]$Kingdom, useNA = "always")
#'    MCUFD_tmp[[1]]$Species[1:10]
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
#' @param save Logical, should the figure be saved to file? Default = FALSE.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#'    virusSet <- readSeq(example = TRUE)
#'    virusCF <- codonFreq(virusSet)
#'    exclCod <- c("ATT", "TGT")
#'    MCUFD_tmp <- MCUFD(virusCF, exclude = exclCod, norm = TRUE)
#'    MCUFD_plot(
#'        cFres = MCUFD_tmp, type = "bar", save = TRUE,
#'        fname = "MCUFD_bar_kingdom", n = 100, rank = "Kingdom"
#'    )
#'
#' @rdname plots
#'
#' @export
setMethod("MCUFD_plot",
    signature(cFres = "list"),
    function(cFres, type, n, rank, fname, units, width, height, dpi, save) {
        keepCols <- switch(rank,
            "Domain" = c(4, ncol(cFres[[1]])-1, ncol(cFres[[1]])),
            "Kingdom" = c(5,  ncol(cFres[[1]])-1, ncol(cFres[[1]])),
            "Phylum" = c(6,  ncol(cFres[[1]])-1, ncol(cFres[[1]]))
        )
        if (is.na(n)) {
            cFres <- lapply(cFres, "[", , keepCols, drop = FALSE)
        } else {
            cFres <- lapply(cFres, "[", 1:n, keepCols, drop = FALSE)
        }
        cFres <- dplyr::bind_rows(cFres)
        cFres$seqid <- factor(
            str_sub(cFres$seqid, 1, 25),
            levels = unique(str_sub(cFres$seqid, 1, 25))
        )
        nseq <- nlevels(cFres$seqid)
        ## get the x most common taxon ranks:
        x <- 5
        keepTax <- sort(
            as.data.frame(
                cFres %>% count_(rank, sort = TRUE) %>% top_n(n = x, wt = n)
            )[1:x,1]
        )
        cFtop <- subset(cFres, cFres[,1] %in% keepTax)
        cFres$retax <- ifelse(
            cFres[, rank] %in% keepTax, cFres[, rank], "Other"
        )
        cFres$retax <- factor(
            cFres$retax, levels = c(keepTax, "Other")
        )
        cc1 <- 12
        brewer_pallette <- c(
            brewer.pal(8, "Set1")[1], brewer.pal(8, "Set1")[5],
            brewer.pal(8, "Set1")[2:4], brewer.pal(8, "Set1")[6:9]
        )
        if (type == "bar") {
            plt <- ggplot(cFres, aes(x = "seqid", fill = "retax")) +
                geom_bar(position = "fill") +
                theme_classic() +
                scale_fill_manual(
                    values = c(brewer_pallette[1:x], "grey"),
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
                ylab(paste0("Proportion of top ", n, " taxa")) +
                scale_x_discrete(expand = c(0.03, 0)) +
                scale_y_continuous(expand = c(0.02, 0))
        } else if (type == "histogram") {
            breaks <- pretty(
                range(cFres$MCUFD), n = (nclass.FD(cFres$MCUFD))*5, min.n = 1
            )
            bwidth <- breaks[2] - breaks[1]
            plt <- ggplot(cFres,
                aes(
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
        if (isTRUE(save)) {
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
        foldEnrich <- (ratioSub - ratioFull)
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
#' @param save Logical, should the enrichment plot be saved to file?
#'    Default = FALSE.
#' @param pthresh Numeric, adjusted p-value threshold to be used for subsetting
#'    the data if plot == TRUE. Default = NA.
#' @param fname Character, name of the output file (used if plot == TRUE).
#' @param units Numeric, units to be used for defining the plot size.
#'    Options are "in" (default), "cm", and "mm".
#' @param width Numeric, width of the figure (in \code{units}).
#' @param height Numeric, height of the figure (in \code{units}).
#' @param dpi Numeric, resolution of the figure (default = 600).
#' @param ptype Character, type of plot to make (if plot == TRUE).  Options are
#'    "heatmap" (default) and "dotplot".
#' @param outtab Character (optional); name of file to which enrichment test
#'     results will be saved. Format is tab-delimited.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#'    virusSet <- readSeq(example = TRUE)
#'    virusCF <- codonFreq(virusSet)
#'    exclCod <- c("ATT", "TGT")
#'    MCUFD_tmp <- MCUFD(virusCF, exclude = exclCod, norm = TRUE)
#'    enrich_tmp <- MCUFD_enrich(
#'        MCUFD_tmp, n = 100, rank = "Phylum", save = FALSE, pthresh = 0.01
#'    )
#'
#' @export
setMethod("MCUFD_enrich",
    signature(cFres = "list"),
    function(
        cFres, n, rank, save, pthresh, fname,
        units, width, height, dpi, ptype, outtab
    ) {
        resList <- vector("list", length(cFres))
        sapply(seq_along(cFres), function(i) {
            resList[[i]] <<- .enrichTest(cFres[[i]], n, rank)
        })
        resMerge <- dplyr::bind_rows(resList, .id = "seqid")
        if (!is.na(pthresh)) {
            # resMerge <- subset(resMerge, padj < pthresh)
            keepTax <- unique(resMerge$Taxon[resMerge$padj < pthresh])
            resMerge <- subset(resMerge, Taxon %in% keepTax)
            resMerge$padj <- ifelse(
                resMerge$padj < pthresh, resMerge$padj, 0
            )
        }
        if (!is.null(outtab)) {
            write.table(
                resMerge, file = paste0(outtab, ".tsv"),
                quote = FALSE, sep = "\t", row.names = FALSE
             )
        }
        resMelt <- melt(
            resMerge, id.vars = c("Taxon", "seqid"),
            measure.vars = "foldEnrich"
        )
        subMelt <- melt(
            resMerge, id.vars = c("Taxon"), measure.vars = "foldEnrich"
        )
        resMelt$Taxon <- factor(
            resMelt$Taxon, levels = rev(levels(factor(resMelt$Taxon)))
        )
        resMelt$seqid <- factor(
            resMelt$seqid, levels = c(1:length(resMelt$seqid))
        )
        subMelt$Taxon <- factor(
            subMelt$Taxon, levels = rev(levels(factor(subMelt$Taxon)))
        )
        cc1 <- 12
        heatpal <- brewer.pal(3, "Set1")[c(3,1)]
        heatlim <- c(
            floor(min(resMelt$value)), ceiling(max(resMelt$value))
        )
        plt <- switch(ptype,
            "dotplot" = ggplot(subMelt, aes(x = Taxon, y = value)) +
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
                ylab(paste0("Fold enrichment in top ", n, " taxa")),
            "heatmap" = ggplot(
                resMelt, aes(x = Taxon, y = seqid, fill = value)
                ) +
                geom_tile(colour = "black") +
                scale_fill_gradient2(
                    paste0("Fold enrichment\nin top ", n, " taxa"),
                    low = heatpal[1],
                    mid = "white", high = heatpal[2],
                    midpoint = 0, limits = heatlim
                ) +
                scale_y_discrete(
                    name = "Sequence",
                    labels = str_sub(names(cFres), 1, 20)
                ) +
                labs(x = rank) +
                theme_bw() +
                theme(
                    text = element_text(size = cc1),
                    # panel.grid.major = element_blank(),
                    # panel.grid.minor = element_blank(),
                    axis.text.x = element_text(
                        angle = 45, hjust = 1, size = cc1*1.5,
                        #hjust = -0.5, #vjust = 0.5,
                        margin = margin(2,0,0,0)
                    ),
                    axis.text.y = element_text(
                        size=cc1, margin = margin(0,2,0,0)
                    ),
                    axis.title.x = element_text(
                        colour = "black", #size=cc1,
                        #hjust = -0.5, #vjust = 0.5,
                        margin = margin(2,0,0,0)
                    ),
                    legend.position = "top",
                ) +
                ylab(rank)
        )
        if (isTRUE(save)) {
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
        plt
    }
)
