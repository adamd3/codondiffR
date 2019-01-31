#' @import methods
#' @import Biostrings
NULL

#' A class for storing relative codon frequencies in user-supplied sequences.
#'
#' @slot seqID Character vector of sequence identifiers, derived from Fasta
#'    identifier lines.
#' @slot codon_perc Matrix containing the relative frequencies of codons; each
#'    sequence is a row and each codon is a column.
#' @slot ncod Numeric vector of sequence lengths (in codons).
setClass(
    "codonFreq",
    slots = c(
        seqID = "character",
        freq = "matrix",
        ncod = "numeric"
    ),
    prototype = list(
        name = NA_character_,
        freq = matrix(nrow = 0, ncol = 0),
        ncod = NA_real_
    )
)

setValidity(
    "codonFreq",
    function(object) {
        if (is.numeric(object@freq)) TRUE else
            stop("Codon frequencies must be numeric values. \n")
})
