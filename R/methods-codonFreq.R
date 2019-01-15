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
        ntlen = width(object)
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
#' @param object An object of class \code{codonFreq}.
#'
#' @return Numeric, the number of sequences in the \code{codonFreq} object.
setMethod("nseq", "codonFreq", function(object) nrow(getFreqs(object)))

#' @describeIn codonFreq Returns the sequence identifiers in a \code{codonFreq}
#'    object.
#'
#' @param object An object of class \code{codonFreq}.
#'
#' @return Character, the identifiers of sequences in the \code{codonFreq}
#'    object.
setMethod("seqID", "codonFreq", function(object) return(object@seqID))

#' @describeIn codonFreq Returns the lengths of sequences (nucleotides) in a
#'    \code{codonFreq} object.
#'
#' @param object An object of class \code{codonFreq}.
#'
#' @return Numeric, lengths of sequences (in nucleotides).
setMethod("seqlen", "codonFreq", function(object) return(object@ntlen))
