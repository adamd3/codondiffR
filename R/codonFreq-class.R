#' @import methods
#' @importFrom Biostrings trinucleotideFrequency
NULL

#' An S4 class to store relative codon frequencies for user-supplied sequences.
#'
#' @slot seqID Character vector of sequence identifiers, derived from Fasta
#'    identifier lines.
#' @slot codon_perc Matrix containing the relative frequencies of codons; each
#'    sequence is a row and each codon is a column.
#' @slot ntlen Numeric vector of sequence lengths (in nucleotides).
setClass(
    "codonFreq",
    slots = c(
        seqID = "character",
        freq = "matrix",
        ntlen = "numeric"
    ),
    prototype = list(
        name = NA_character_,
        freq = matrix(nrow = 0, ncol = 0),
        ntlen = NA_real_
    )
)

setValidity(
    "codonFreq",
    function(x) {
        errors <- character()
        if (!is.numeric(x@freq)) {
            error <- "Codon frequencies must be numeric values. \n"
            errors <- c(errors, error)
        }
        if (!all(width(tmp)%%3 == 0)) {
            error <- paste0(
                "All sequence lengths must be a multiple of three. ",
                "Partial codons and non-coding regions are not allowed. \n"
            )
            errors <- c(errors, error)
        }
        if (length(errors) == 0) TRUE else stop(errors)
    }
)

.codonFreq <- function(x) {
  trinucleotideFrequency(x, step = 3, as.prob = TRUE)
}

##------------------------------------------------------------------------------
## codonFreq constructor
##------------------------------------------------------------------------------
#' @rdname codonFreq-class
#' @export
setGeneric("codonFreq", function(x) standardGeneric("codonFreq"))

#' Create new objects of class \code{codonFreq}.
#'
#' @describeIn codonFreq
#'
#' @param x An object of class \code{DNAStringSet}.
#'
#' @return A \code{codonFreq} object.
setMethod("codonFreq", "DNAStringSet", function(x) {
    emptyidx <- which(nchar(x) == 0)
    if (length(emptyidx) > 0) {
      warning(
          paste0(
              "\nThe following sequences are empty and will be ignored: \n",
              names(x)[emptyidx],
              "\n"
          )
      )
      x <- x[-emptyidx]
    }
    freqmat <- .codonFreq(x)
    freqmat <- freqmat[,order(colnames(freqmat))]
    new(
        "codonFreq",
        seqID = names(x),
        freq = freqmat,
        ntlen = width(x)
    )
})



##------------------------------------------------------------------------------
## Show method
##------------------------------------------------------------------------------
#' Display the data in a \code{codonFreq} object.
#'
#' @describeIn codonFreq
#'
#' @param object A \code{codonFreq} object.
#'
#' @export
setMethod("show", "codonFreq", function(x) {
    sc <- length(c1@seqID)
    cat(
        "An object of class codonFreq, with data for",
        sc,
        "sequences.\n"
    )
    cdf <- as.data.frame(x@freq)
    rownames(cdf) <- c1@seqID
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
#' @rdname codonFreq-class
#'
#' @export
setGeneric("getFreqs", function(x) standardGeneric("getFreqs"))

#' Get relative frequencies from a \code{codonFreq} object.
#'
#' Returns the relative frequencies of codons per sequence in a \code{codonFreq}
#'    object.
#'
#' @describeIn codonFreq
#'
#' @param x An object of class \code{codonFreq}.
#'
#' @return Matrix, frequencies of codons in the \code{codonFreq} object.
#'    Codons are in columns and input sequences are in rows.
#'
#' @export
setMethod("getFreqs", "codonFreq", function(x) return(x@freq))


#' @rdname codonFreq-class
#'
#' @export
setGeneric("nseq", function(x) standardGeneric("nseq"))

#' Number of sequences in a \code{codonFreq} object.
#'
#' Returns the number of sequences in a \code{codonFreq} object.
#'
#' @describeIn codonFreq
#'
#' @export
#'
#' @inheritParams getFreqs
#'
#' @return Numeric, the number of sequences in the \code{codonFreq} object.
setMethod("nseq", "codonFreq", function(x) nrow(getFreqs(x)))


#' @rdname codonFreq-class
#'
#' @export
setGeneric("seqID", function(x) standardGeneric("seqID"))

#' Get the sequence identifiers in a \code{codonFreq} object.
#'
#' Returns the sequence identifiers in a \code{codonFreq} object.
#'
#' @describeIn codonFreq
#'
#' @export
#'
#' @inheritParams getFreqs
#'
#' @return Character, the identifiers of sequences in the \code{codonFreq}
#'    object.
setMethod("seqID", "codonFreq", function(x) return(x@seqID))

#' @rdname codonFreq-class
#'
#' @export
setGeneric("seqlen", function(x) standardGeneric("seqlen"))

#' Get the lengths of sequences (nucleotides) in a \code{codonFreq} object.
#'
#' Returns the lengths of sequences (nucleotides) in a \code{codonFreq} object.
#'
#' @describeIn codonFreq
#'
#' @export
#'
#' @inheritParams getFreqs
#'
#' @return Numeric, lengths of sequences (in nucleotides).
setMethod("seqlen", "codonFreq", function(x) return(x@ntlen))
