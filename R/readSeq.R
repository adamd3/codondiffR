#' Read nucleotide sequences from a file.
#'
#' Reads a fasta-format \code{file} containing one or more DNA or RNA sequences.
#'
#' @param file Character, path to a file containing sequences in fasta format.
#' @param type Character, indicates the nucleic acid type (DNA or RNA). The
#'     default is DNA.
#' @param ... Other arguments passed to readDNAStringSet or readRNAStringSet.
#'
#' @return Returns a \code{DNAStringSet} object or an \code{RNAStringSet}
#'     object.
#'
#' @examples
#' readSeq(file = "example.fna", type = "RNA")
#'
#' @importFrom Biostrings readDNAStringSet
#' @importFrom Biostrings readRNAStringSet
#'
#' @export
readSeq <- function(file = character(), type = "DNA", ...) {
    if (type == "DNA") {
        SSobj <- readDNAStringSet(file, format="fasta")
    } else if (type == "RNA") {
        SSobj <- readRNAStringSet(file, format="fasta")
    }
    SSobj
}

# to get lengths of sequences:
# fasta.seqlengths("reverse.ORF.RdRps.representative.fna")
