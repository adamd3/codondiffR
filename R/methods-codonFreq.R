#' @include AllClasses.R
#' @include AllGenerics.R
#' @import ggplot2
#' @import reshape2
#' @import MASS
#' @import caret
#' @import e1071
#' @import RColorBrewer
#' @import dplyr
#' @import stringr
#' @import scales
#' @import ggbiplot
#' @importFrom Biostrings trinucleotideFrequency
#' @importFrom Biostrings readDNAStringSet
#' @importFrom Biostrings width
NULL


##------------------------------------------------------------------------------
## Validity check
##------------------------------------------------------------------------------
setMethod(.validity, "codonFreq", function(object) {
    errors <- character()
    if (!is.numeric(object@freq)) {
        msg <- paste("Codon frequencies must be numeric values. \n")
        errors <- c(errors, msg)
    }
    if (length(errors) == 0) TRUE else errors
})



##------------------------------------------------------------------------------
## Read in sequences
##------------------------------------------------------------------------------
#' Read nucleotide sequences from a file.
#'
#' Reads a fasta-format file containing one or more protein-coding
#'     DNA sequences.
#'
#' @param file Character, path to a file containing one or more protein-coding
#'     DNA sequences in fasta format.
#' @param ... Other arguments passed to readDNAStringSet.
#'
#' @return Returns a \code{DNAStringSet} object.
#'
#' @examples
#'     mySeq <- readSeq(file = "example.fasta")
#'
#' @export
readSeq <- function(file = character(), ...) {
    SSobj <- readDNAStringSet(file, format="fasta")
    SSobj
}



##------------------------------------------------------------------------------
## codonFreq constructor
##------------------------------------------------------------------------------
#' Create new objects of class \code{codonFreq}.
#'
#' @name codonFreq
#' @rdname codonFreq
#'
#' @param object An object of class \code{DNAStringSet}.
#'
#' @return a \code{codonFreq} object.
setMethod("codonFreq", "DNAStringSet", function(object) {
    if (!all(width(object) %% 3 == 0)) {
        stop(paste0(
            "All sequences lengths must be a multiple of 3.\n",
            "Partial codons and non-coding regions are not allowed."
        ))
    }
    emptyidx <- which(width(object) == 0)
    if (length(emptyidx) > 0) {
      warning(
          paste0(
              "\nThe following sequences are empty and will be ignored: \n",
              names(object)[emptyidx],
              "\n"
          )
      )
      object <- object[-emptyidx]
    }
    freqmat <- Biostrings::trinucleotideFrequency(
        object, step = 3, as.prob = TRUE
    )
    freqmat <- freqmat[,order(colnames(freqmat))]
    new(
        "codonFreq",
        seqID = names(object),
        freq = freqmat,
        ncod = width(object)/3
    )
})



##------------------------------------------------------------------------------
## Show method
##------------------------------------------------------------------------------
#' Display the data in a \code{codonFreq} object.
#'
#' @rdname show-codonFreq
#'
#' @param object A \code{codonFreq} object.
#'
#' @export
setMethod("show", "codonFreq", function(object) {
    sc <- length(object@seqID)
    cat(
        "An object of class codonFreq, with data for",
        sc,
        "sequences.\n"
    )
    cdf <- as.data.frame(object@freq)
    rownames(cdf) <- object@seqID
    if (sc < 7) {
        print(cdf)
    } else {
        outdf <- rbind(cdf[1:3,], "..." = " ", cdf[(sc-2):sc,])
        print(outdf)
    }
})



##------------------------------------------------------------------------------
## Codon usage plots
##------------------------------------------------------------------------------
#' Plot the relative frequencies of codons in a \code{codonFreq} object.
#'
#' @rdname plot-codonFreq
#'
#' @param object A \code{codonFreq} object.
#' @param fname Character, name of figure generated.
#' @param units Numeric, units to be used for defining the plot size.
#'    Options are "in" (default), "cm", and "mm".
#' @param width Numeric, width of the figure (in \code{units}).
#' @param height Numeric, height of the figure (in \code{units}).
#' @param dpi Numeric, resolution of the figure (default = 600).
#' @param groups Character, optional vector giving groups for each of the
#'    sequences in the \code{codonFreq} object. If not supplied, the sequences
#'    will be treated as a single group.
#' @param ptype Character, type of plot to make.  Options are "boxplot"
#'    (default) and "jitter".
#' @param order Character, determines how codons are order. The default
#'    (median) sorts codons by median usage in the full dataset. An alternative
#'    option is "alphabetical". A vector of (64) codons can also be passed,
#'    which will be used for ordering, exactly as in the vector.
#' @param colour Integer, between 1 and 9, which specifies the colour to be used
#'    from the `Set1` palette in the RColorBrewer package. Only applies to
#'    non-grouped plots. Default = 1 (red).
#' @param suppress_x_txt Logical, suppress x axis labels? Default = FALSE.
#' @param suppress_y_title Logical, suppress y axis title? Default = FALSE.
#' @param label Character, optional label to add to the plotting area.
#' @param highlight Character (optional) codons in this vector will be
#'    highlighted.
#' @param save Logical, should the plot be saved to file? Default
#'    = FALSE.
#'
#' @return A \code{ggplot} object.
#'
#' @export
setMethod(
    "codonPlot",
    "codonFreq",
    function(
        object, fname, units, width, height, dpi, groups, ptype, order,
        colour, suppress_x_txt, suppress_y_title, label, highlight, save
    ) {
        sc <- length(object@seqID)
        cdf <- as.data.frame(object@freq)
        cdf$Taxon <- object@seqID
        ## median usage per codon
        codMeds <- apply(cdf[, 1:(ncol(cdf)-1)], 2, FUN = median)
        brewer_pallette1 <- brewer.pal(9,"Set1")
        cc1 <- 12
        if (isTRUE(suppress_x_txt)) {
            xtxt = element_blank()
        } else {
            xtxt = element_text(
                colour = "black", size=cc1*1.2,
                angle = 90, hjust = 0.1, vjust = 0.5,
                margin = margin(10,0,0,0)
            )
        }
        if (isTRUE(suppress_y_title)) {
            ytit = element_blank()
        } else {
            ytit = element_text(
                margin = margin(0,20,0,0), size=cc1
            )
        }
        if (is.null(groups)) {
            cols <- rep(brewer_pallette1[colour], 64)
            codon_melt <- melt(cdf, variable.id = c("Taxon"))
            codon_melt$variable <- gsub("T", "U", codon_melt$variable)
            if (order[1] == "alphabetical") {
                codon_melt$variable <- factor(
                    codon_melt$variable,
                    levels = rev(unique(codon_melt$variable))
                )
            } else if (length(order) == 64) {
                codon_melt$variable <- factor(
                    codon_melt$variable,
                    levels = order
                )
            } else {
                codon_melt$variable <- factor(
                    codon_melt$variable,
                    levels = gsub(
                        "T", "U", names(sort(codMeds, decreasing = TRUE))
                    )
                )
            }
            p1 <- switch(ptype,
                "boxplot" = ggplot(
                    codon_melt, aes(x = variable, y = value, fill = variable)
                    ) +
                    geom_boxplot(
                        outlier.size = 2, lwd = 1#, fatten = 1
                    ) +
                    scale_fill_manual("", values = cols),
                "jitter" = ggplot(
                    codon_melt, aes(x = variable, y = value, colour = variable)
                    ) +
                    geom_jitter(position = position_jitter(0.3), size = 1) +
                    scale_colour_manual("", values = cols) +
                    guides(colour = guide_legend(override.aes = list(size=3)))
                )
            p1 <- p1 +
                labs(y = "Proportion of codons", x = "Codon") +
                theme_bw() +
                # coord_flip() +
                theme(
                    axis.title.x = element_blank(),
                    axis.title.y = ytit,
                    axis.text.y = element_text(colour = "black", size=cc1*1.2),
                    legend.position = "none"
                )
            if (!is.null(label)) {
                p1 <- p1 + annotate(
                    "text", label = label, x = 64,
                    y = max(codon_melt$value)-0.001, size = 6,
                    hjust = 1
                )
            }
            if (!is.null(highlight)) {
                highlight_idx <- which(
                    levels(codon_melt$variable) %in% highlight
                )
                fontFace <- rep("plain", nlevels(codon_melt$variable))
                fontFace[highlight_idx] <- "bold"
                fontCols <- rep("black", nlevels(codon_melt$variable))
                fontCols[highlight_idx] <- "red"
                print(highlight_idx)
                print(fontFace)
                print(levels((codon_melt$variable)))
                p1 <- p1 + theme(
                    axis.text.x = element_text(
                        colour = fontCols, size=cc1*1.2,
                        face = fontFace,
                        angle = 90, hjust = 0.1, vjust = 0.5,
                        margin = margin(10,0,0,0)
                    )
                )
            } else {
                p1 <- p1 + theme(axis.text.x = xtxt)
            }
        } else if (length(groups) == sc) {
            cdf <- cbind(cdf, groups)
            cdf$groups <- as.factor(cdf$groups)
            cols <- brewer_pallette1[1:nlevels(cdf$groups)]
            codon_melt <- melt(cdf, variable.id = c("Taxon", "groups"))
            codon_melt$variable <- gsub("T", "U", codon_melt$variable)
            if (order[1] == "alphabetical") {
                codon_melt$variable <- factor(
                    codon_melt$variable,
                    levels = rev(unique(codon_melt$variable))
                )
            } else if (length(order) == 64) {
                codon_melt$variable <- factor(
                    codon_melt$variable,
                    levels = order
                )
            } else {
                codon_melt$variable <- factor(
                    codon_melt$variable,
                    levels = gsub(
                        "T", "U", names(sort(codMeds, decreasing = TRUE))
                    )
                )
            }
            p1 <- switch(ptype,
                "boxplot" = ggplot(
                    codon_melt, aes(x = variable, y = value, fill = groups)
                    ) +
                    geom_boxplot(
                        outlier.size = 0.5, lwd = 0.5, fatten = 1,
                        position = position_dodge(0.7), alpha = 0.7
                    ) + scale_fill_manual("", values = cols)
                    ,
                "jitter" = ggplot(
                    codon_melt, aes(x = variable, y = value, colour = groups)
                    ) +
                    geom_jitter(position = position_jitter(0.3), size = 1) +
                    scale_colour_manual("", values = cols) +
                    guides(colour = guide_legend(override.aes = list(size=3)))
            )
            p1 <- p1 +
                labs(y = "Proportion of codons", x = "Codon") +
                theme_bw() +
                # coord_flip() +
                theme(
                    text = element_text(size=cc1),
                    axis.text.x = element_text(
                        colour = "black", size=cc1,
                        angle = 90, hjust = 0.1, vjust = 0.5,
                        margin = margin(10,0,0,0)
                    ),
                    # axis.title.x = element_text(
                    #     colour = "black", size=cc1,
                    #     hjust = 0.5, #vjust = 0.1
                    #     margin = margin(20,0,0,0)
                    # ),
                    axis.title.x = element_blank(),
                    axis.title.y = element_text(
                        margin = margin(0,20,0,0), size=cc1
                    ),
                    axis.text.y = element_text(colour = "black", size=cc1),
                    legend.position = "right",
                    legend.text = element_text(size = cc1)
                )
        } else {
            stop("`groups` vector must be same length as codonFreq object")
        }
        if (isTRUE(save)) {
            ggsave(
                p1,
                file = paste0(fname, ".png"),
                device = "png",
                units = units,
                width = width,
                height = height,
                dpi = dpi
            )
        }
        p1
})


#' Boxplots of normalised codon usage by amino acid (i.e codon bias) in a
#'    \code{codonFreq} object.
#'
#' @rdname plot-codonFreq
#'
#' @param object A \code{codonFreq} object.
#' @param fname Character, name of figure generated.
#' @param units Numeric, units to be used for defining the plot size.
#'    Options are "in" (default), "cm", and "mm".
#' @param width Numeric, width of the figure (in \code{units}).
#' @param height Numeric, height of the figure (in \code{units}).
#' @param dpi Numeric, resolution of the figure (default = 600).
#' @param groups Character, optional vector giving groups for each of the
#'    sequences in the \code{codonFreq} object. If not supplied, the sequences
#'    will be treated as a single group.
#' @param aa Character, amino acid for which codon bias plot will be made
#'    (must be specified using single letter IUPAC code).
#' @param norm Logical, should the data be normalised? If TRUE, a normalisation
#'    step will be performed. Default = FALSE.
#' @param label Character, optional label to add to the plotting area.
#' @param colours Optional vector of same length as `groups`, containing integers
#'    between 1 and 9, which specifies the colours to be used from the `Set1`
#'    palette in the RColorBrewer package. By default, the Set1 order is used.
#' @param suppress_y_txt Logical, suppress y axis labels? Default = FALSE.
#' @param suppress_y_title Logical, suppress y axis title? Default = FALSE.
#' @param legend Logical, plot a legend? Default = TRUE.
#' @param ylim Numeric vector, gives the y-axis limits (if not supplied, they
#'    will be chosen based on the data).
#' @param save Logical, plot the plot be saved to file? Default = FALSE.
#'
#' @return A \code{ggplot} object.
#'
#' @export
setMethod(
    "biasPlot",
    "codonFreq",
    function(
        object, fname, units, width, height, dpi, groups, aa, norm, label,
        colours, suppress_y_txt, suppress_y_title, legend, ylim, save
    ) {
        if (isTRUE(norm)) object <- normalise(object)
        sc <- length(object@seqID)
        cdf <- as.data.frame(object@freq)
        if (is.null(aa)) {
            stop("Amino acid must be specified, using the `aa` parameter")
        } else {
            keepCod <- rownames(stdgc)[grepl(aa, stdgc$AA)]
            cdf <- subset(cdf, select = keepCod)
        }
        cdf$Taxon <- object@seqID
        brewer_pallette1 <- brewer.pal(9,"Set1")
        cc1 <- 10
        if (isTRUE(suppress_y_txt)) {
            ytxt = element_blank()
        } else {
            ytxt = element_text(
                colour = "black", size=cc1,
                margin = margin(0,5,0,0)
            )
        }
        if (isTRUE(suppress_y_title)) {
            ytit = element_blank()
        } else {
            ytit = element_text(
                margin = margin(0,5,0,0), size=cc1
            )
        }
        if (length(groups) == sc) {
            cdf <- cbind(cdf, groups)
            cdf$groups <- factor(cdf$groups, levels = unique(groups))
        } else {
            stop("`groups` vector must be same length as codonFreq object")
        }
        if (is.null(colours)) {
            cols <- brewer_pallette1[1:nlevels(cdf$groups)]
        } else if (length(colours) == nlevels(cdf$groups)) {
            cols <- brewer_pallette1[colours]
        } else {
            stop("`colours` vector must be same length as number of groups")
        }
        codon_melt <- melt(cdf, variable.id = c("Taxon", "groups"))
        codon_melt$variable <- gsub("T", "U", codon_melt$variable)
        if (is.null(ylim)) {
            y1 <- c(
                min(codon_melt$value)-0.05, max(codon_melt$value)+0.05
            )
        } else {
            y1 <- ylim
        }
        lgd <- ifelse(isTRUE(legend), "right", "none")
        p1 <- ggplot(
                codon_melt, aes(x = variable, y = value, fill = groups)
            ) +
            geom_boxplot(
                outlier.size = 2, lwd = 1, #fatten = 1,
                position = position_dodge(1), alpha = 1
            ) +
            scale_fill_manual("", values = cols) +
            labs(y = "Proportion of codons", x = "Codon") +
            theme_bw() +
            ylim(y1) +
            theme(
                legend.text = element_text(size = cc1),
                text = element_text(size=cc1),
                axis.text.x = element_text(
                    colour = "black", size=cc1,
                    # angle = 90, hjust = 0.1, vjust = 0.5
                    margin = margin(5,0,0,0)
                ),
                axis.title.x = element_blank(),
                axis.title.y = ytit,
                axis.text.y = ytxt,
                legend.position = lgd
            )
        if (!is.null(label)) {
            p1 <- p1 + annotate(
                "text", label = label, x = length(keepCod)+0.5,
                y = y1[2]-0.01, size = 4.5,
                hjust = 1
            )
        }
        if (isTRUE(save)) {
            ggsave(
                p1,
                file = paste0(fname, ".png"),
                device = "png",
                units = units,
                width = width,
                height = height,
                dpi = dpi
            )
        }
        p1
})


#' GC3 plots: plot the 3rd-codon GC content for sequences in a \code{codonFreq}
#'   object.
#'
#' @rdname plot-codonFreq
#'
#' @param object A \code{codonFreq} object.
#' @param fname Character, name of figure generated.
#' @param units Numeric, units to be used for defining the plot size.
#'    Options are "in" (default), "cm", and "mm".
#' @param width Numeric, width of the figure (in \code{units}).
#' @param height Numeric, height of the figure (in \code{units}).
#' @param dpi Numeric, resolution of the figure (default = 600).
#' @param groups Character, optional vector giving groups for each of the
#'    sequences in the \code{codonFreq} object. If not supplied, the sequences
#'    will be treated as a single group.
#' @param aa Character, amino acid for which codon bias plot will be made
#'    (must be specified using single letter IUPAC code).
#' @param norm Logical, should the data be normalised? If TRUE, a normalisation
#'    step will be performed. Default = FALSE.
#' @param label Character, optional label to add to the plotting area.
#' @param colours Optional vector of same length as `groups`, containing integers
#'    between 1 and 9, which specifies the colours to be used from the `Set1`
#'    palette in the RColorBrewer package. By default, the Set1 order is used.
#' @param suppress_y_txt Logical, suppress y axis labels? Default = FALSE.
#' @param suppress_y_title Logical, suppress y axis title? Default = FALSE.
#' @param legend Logical, plot a legend? Default = TRUE.
#' @param xlim Numeric vector, gives the x-axis limits (if not supplied, they
#'    will be chosen based on the data).
#' @param outtab Character (optional); name of file to which GC3 data will be
#'     saved. Format is tab-delimited.
#' @param save Logical, save the plot to file? Default = FALSE.
#'
#' @return A \code{ggplot} object.
#'
#' @export
setMethod(
    "gcPlot",
    "codonFreq",
    function(
        object, fname, units, width, height, dpi, groups, aa, norm, label,
        colours, suppress_y_txt, suppress_y_title, legend, xlim, outtab, save
    ) {
        sc <- length(object@seqID)
        cdf <- as.data.frame(object@freq)
        cdf$Taxon <- object@seqID
        brewer_pallette1 <- brewer.pal(9,"Set1")
        cc1 <- 12
        if (isTRUE(suppress_y_txt)) {
            ytxt = element_blank()
        } else {
            ytxt = element_text(
                colour = "black", size=cc1*1.2,
                margin = margin(0,10,0,0)
            )
        }
        if (isTRUE(suppress_y_title)) {
            ytit = element_blank()
        } else {
            ytit = element_text(
                margin = margin(0,10,0,0), size=cc1*1.2
            )
        }
        if (is.null(groups)) {
            stop(
                "groups` vector must be either length 1 or same length as
                codonFreq object"
            )
        } else {
            cdf <- cbind(cdf, groups)
            cdf$groups <- factor(cdf$groups, levels = sort(unique(groups)))
        }
        if (is.null(colours)) {
            cols <- brewer_pallette1[1:nlevels(cdf$groups)]
        } else if (length(colours) == nlevels(cdf$groups)) {
            cols <- brewer_pallette1[colours]
        } else {
            stop("`colours` vector must be same length as number of groups")
        }
        ## get GC3% per taxon
        third_gc_idx <- which(
            str_sub(colnames(cdf)[1:64], -1, -1) %in% c("G", "C")
        )
        cdfGC <- cbind(
            (cdf)[,65:ncol(cdf)],
            rowSums(cdf[,third_gc_idx])
        )
        colnames(cdfGC) <- c("Taxon", "groups", "GC3_proportion")
        codon_melt <- melt(cdfGC, variable.id = c("Taxon", "groups"))
        if (is.null(xlim)) {
            x1 <- c(
                min(codon_melt$value)-0.05, max(codon_melt$value)+0.05
            )
        } else {
            x1 <- xlim
        }
        lgd <- ifelse(isTRUE(legend), "right", "none")
        if (!is.null(outtab)) {
            write.table(
                cdfGC, file = paste0(outtab, ".tsv"),
                quote = FALSE, sep = "\t", row.names = FALSE
             )
        }
        group_means <- cdfGC %>%
            group_by(groups) %>%
            dplyr::summarize(mean = mean(GC3_proportion, na.rm=TRUE))
        group_means <- group_means$mean
        p1 <- ggplot(
            codon_melt, aes(x = value, colour = groups)
            ) +
            geom_density(size = 1) +
            theme_bw() +
            xlim(x1) +
            geom_vline(
                xintercept = c(group_means),
                linetype = "dashed",
                color = c(cols),
                size = 1
            ) +
            # scale_y_continuous(
            #     labels = function(x) format(
            #         x, scientific = FALSE, digits = 1, nsmall = 4
            #     )
            # ) +
            labs(
                x = "GC3 proportion in codons", y = "Density"
            ) +
            scale_colour_manual("Group", values = cols) +
            theme(
                text = element_text(size=cc1),
                axis.text.x = element_text(colour = "black", size=cc1),
                axis.title.x = element_text(
                    colour = "black", size=cc1, margin = margin(5,0,0,0)
                ),
                axis.title.y = element_text(
                    colour = "black", size=cc1, margin = margin(0,5,0,0)
                ),
                axis.text.y = element_text(colour = "black", size=cc1),
                legend.position = "right"
            )
        if (isTRUE(save)) {
            ggsave(
                p1,
                file = paste0(fname, ".png"),
                device = "png",
                units = units,
                width = width,
                height = height,
                dpi = dpi
            )
        }
        p1
})



##------------------------------------------------------------------------------
## Accessor methods
##------------------------------------------------------------------------------
#' @describeIn codonFreq Returns the relative frequencies of codons per sequence
#'    in a \code{codonFreq} object.
#'
#' @param object An object of class \code{codonFreq}.
#'
#' @return Matrix, frequencies of codons in the \code{codonFreq} object.
#'    Codons are in columns and input sequences are in rows.
setMethod("getFreqs", "codonFreq", function(object) return(object@freq))

#' @describeIn codonFreq Returns the number of sequences in a \code{codonFreq}
#'    object.
#'
#' @inheritParams getFreqs
#'
#' @return Numeric, the number of sequences in the \code{codonFreq} object.
setMethod("nseq", "codonFreq", function(object) nrow(getFreqs(object)))

#' @describeIn codonFreq Returns the sequence identifiers in a \code{codonFreq}
#'    object.
#'
#' @inheritParams getFreqs
#'
#' @return Character, the identifiers of sequences in the \code{codonFreq}
#'    object.
setMethod("seqID", "codonFreq", function(object) return(object@seqID))

#' @describeIn codonFreq Returns the lengths of sequences (in codons) in a
#'    \code{codonFreq} object.
#'
#' @inheritParams getFreqs
#'
#' @return Numeric, lengths of sequences (in codons).
setMethod("seqlen", "codonFreq", function(object) return(object@ncod))


##------------------------------------------------------------------------------
## Subset methods
##------------------------------------------------------------------------------
#' Methods for subsetting a codonFreq object
#
#' @rdname extract-codonFreq
#'
#' @inheritParams base::Extract
#'
#' @seealso [base::Extract]
#'
#' @return A subset of the sequences in the \code{codonFreq} object.
#'
#' @export
setMethod("[", "codonFreq", function(x, i, j) {
    new(
        "codonFreq",
        seqID = x@seqID[i],
        freq = rbind(x@freq[i, j]),
        ncod = x@ncod[i]
    )
})

#' @rdname extract-codonFreq
#'
#' @export
setMethod("[[", "codonFreq", function(x, i, j) {
    new(
        "codonFreq",
        seqID = x@seqID[i],
        freq = rbind(x@freq[i, j]),
        ncod = x@ncod[i]
    )
})



##------------------------------------------------------------------------------
## Normalisation methods
##------------------------------------------------------------------------------
#' Normalise the codon frequencies in a codonFreq object by amino acid
#
#' @name normalise
#' @rdname normalise
#'
#' @param object An object of class \code{codonFreq}.
#'
#' @return An object of class \code{codonFreq}, containing normalised
#'    data in which the sum of the relative frequencies per amino acid for a
#'    given sequence is 1. The Standard Genetic code is currently the only
#'    code used for normalisation.
#'
#' @export
setMethod("normalise", "codonFreq", function(object) {
    cfdf <- cbind(stdgc, data.frame(t(object@freq)))
    cfdf$AA <- as.factor(cfdf$AA)
    cflist <- split(cfdf, cfdf$AA)
    cflist <- lapply(cflist, function(x) {
        sweep(
            x[2:(ncol((cflist)[[1]]))],
            2,
            colSums(x[2:ncol((cflist)[[1]])]),
            FUN="/"
        )
    })
    cfnorm <- do.call("rbind", cflist)
    ## for M and W amino acids, there is only a single codon; therefore, must
    ## re-insert the codons to rownames
    for (uniqaa in c("M", "W")) {
        rownames(cfnorm)[grepl(uniqaa, rownames(cfnorm))] <- paste0(
            rownames(cfnorm)[grepl(uniqaa, rownames(cfnorm))],
            ".",
            rownames(stdgc)[grepl(uniqaa, stdgc$AA)]
        )
    }
    rownames(cfnorm) <- str_split_fixed(rownames(cfnorm), ".", 3)[,3]
    cfnorm <- t(cfnorm[order(row.names(cfnorm)), ])
    cfnorm[is.nan(cfnorm)] <- NA ## NaN values occur when the AA is not used
    new(
        "codonFreq",
        seqID = object@seqID,
        freq = cfnorm,
        ncod = object@ncod
    )
})
