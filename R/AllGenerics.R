##------------------------------------------------------------------------------
## Constructor
##------------------------------------------------------------------------------
#' @rdname codonFreq-class
#'
#' @export
setGeneric("codonFreq", function(object) standardGeneric("codonFreq"))



##------------------------------------------------------------------------------
## Accessors
##------------------------------------------------------------------------------
#' @rdname codonFreq-class
#'
#' @export
setGeneric("getFreqs", function(object) standardGeneric("getFreqs"))

#' @rdname codonFreq-class
#'
#' @export
setGeneric("nseq", function(object) standardGeneric("nseq"))

#' @rdname codonFreq-class
#'
#' @export
setGeneric("seqID", function(object) standardGeneric("seqID"))

#' @rdname codonFreq-class
#'#'
#' @export
setGeneric("seqlen", function(object) standardGeneric("seqlen"))



##------------------------------------------------------------------------------
## Normalisation
##------------------------------------------------------------------------------
#' @rdname codonFreq-class
#'
#' @export
setGeneric("normalise", function(object) standardGeneric("normalise"))



##------------------------------------------------------------------------------
## Codon comparisons
##------------------------------------------------------------------------------
#' @rdname codonCompare
#'
#' @export
setGeneric(
    "MCUFD",
    function(
        cFobj, exclude = character(length = 0), minlen = 600, norm = TRUE
    ) standardGeneric("MCUFD")
)
