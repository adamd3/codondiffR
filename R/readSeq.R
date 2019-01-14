#' Read nucleotide sequences from a file.
#'
#' Reads a fasta-format \code{file} containing one or more protein-coding
#'     DNA sequences.
#'
#' @param file Character, path to a file containing one or more protein-coding
#'     DNA sequences in fasta format.
#' @param ... Other arguments passed to readDNAStringSet or readRNAStringSet.
#'
#' @return Returns a \code{DNAStringSet} object.
#'
#' @examples
#' readSeq(file = "../inst/extdata/example.fna")
#'
#' @importFrom Biostrings readDNAStringSet
#'
#' @export
readSeq <- function(file = character(), ...) {
    SSobj <- readDNAStringSet(file, format="fasta")
    SSobj
}
