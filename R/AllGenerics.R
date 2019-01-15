#' Constructor for objects of \code{codonFreq} class
#'
#' @rdname codonFreq-class
#'
#' @export
setGeneric("codonFreq", function(object) standardGeneric("codonFreq"))

#' Get relative frequencies from a \code{codonFreq} object.
#'
#' @rdname codonFreq-class
#'
#' @export
setGeneric("getFreqs", function(object) standardGeneric("getFreqs"))

#' Number of sequences in a \code{codonFreq} object.
#'
#' @rdname codonFreq-class
#'
#' @export
setGeneric("nseq", function(object) standardGeneric("nseq"))


#' Get the sequence identifiers in a \code{codonFreq} object.
#'
#' @rdname codonFreq-class
#'
#' @export
setGeneric("seqID", function(object) standardGeneric("seqID"))

#' Get the lengths of sequences (nucleotides) in a \code{codonFreq} object.
#'
#' @rdname codonFreq-class
#'
#' @export
setGeneric("seqlen", function(object) standardGeneric("seqlen"))
