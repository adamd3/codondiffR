#' @import methods
#' @importFrom Biostrings trinucleotideFrequency
NULL

#' A class for storing relative codon frequencies in user-supplied sequences.
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
    function(object) {
        errors <- character()
        if (!is.numeric(object@freq)) {
            error <- "Codon frequencies must be numeric values. \n"
            errors <- c(errors, error)
        }
        if (!all(object@ntlen %%3 == 0)) {
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
