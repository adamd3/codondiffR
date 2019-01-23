#' @import Biostrings
#' @import dplyr
#' @import stringr

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
    freqmat <- .codonFreq(object)
    freqmat <- freqmat[,order(colnames(freqmat))]
    new(
        "codonFreq",
        seqID = names(object),
        freq = freqmat,
        seqlen = width(object)
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
setMethod("seqlen", "codonFreq", function(object) return(object@seqlen))


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
        seqlen = x@seqlen[i]
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
        seqlen = x@seqlen[i]
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
        seqlen = object@seqlen
    )
})
