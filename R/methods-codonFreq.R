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
#' @importFrom Biostrings trinucleotideFrequency
#' @importFrom Biostrings readDNAStringSet
#' @importFrom Biostrings width
NULL

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
## Plot method
##------------------------------------------------------------------------------
#' Make a boxplot of the data in a \code{codonFreq} object.
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
#'
#' @export
setMethod(
    "plot",
    "codonFreq",
    function(
        object, fname, units, width, height, dpi) {
        sc <- length(object@seqID)
        cdf <- as.data.frame(object@freq)
        cdf$Taxon <- object@seqID
        codMeds <- apply(cdf[, 1:(ncol(cdf)-1)], 2, FUN = median)
        codon_melt <- melt(cdf, variable.id = c("Taxon"))
        codon_melt$variable <- gsub("T", "U", codon_melt$variable)
        codon_melt$variable <- factor(
            codon_melt$variable,
            levels = gsub("T", "U", names(sort(codMeds, decreasing = TRUE)))
        )
        brewer_pallette1 <- brewer.pal(9,"Set1")
        cols <- rep(brewer_pallette1[1], 64)
        cc1 <- 12
        p1 <- ggplot(
            codon_melt, aes(x = variable, y = value, fill = variable)
        ) +
            geom_boxplot() +
            scale_colour_manual("Codon", values = cols) +
            scale_fill_manual("Codon", values = cols) +
            labs(y = "Proportion of codons", x = "Codon") +
            theme_bw() +
            theme(
                legend.position = "none",
                text = element_text(size=cc1),
                axis.text.x = element_text(
                    colour = "black", size=cc1*0.75, angle=45,
                    margin = margin(15,0,0,0)
                ),
                axis.title.x = element_blank(),
                axis.title.y = element_text(margin = margin(0,15,0,0)),
                axis.text.y = element_text(colour = "black", size=cc1)
            )
        ggsave(
            p1,
            file = paste0(fname, ".png"),
            device = "png",
            units = units,
            width = width,
            height = height,
            dpi = dpi
        )
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
