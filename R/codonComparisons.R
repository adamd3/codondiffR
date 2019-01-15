#' @include AllClasses.R
#' @include AllGenerics.R
#' @include methods-codonFreq.R
NULL

#' Calculate measures of codon usage similarity.
#'
#' Compare codon usage between sequences in a \code{codonFreq} object with those
#'    in a reference database.
#' Mean codon frequency usage difference (\code{MCUFD}) is calculated as
#'    described in \href{https://goo.gl/7sGZeE}{Stedman et al. (2013)}.
#'
#' @param cFobj An object of class \code{codonFreq}.
#' @param exclude A character vector of codons to be excluded from comparisons.
#' @param minlen Numeric (optional), the minimum length of sequence (in
#'    nucleotides) to be included in the analysis.
#'
#' @inheritParams Biostrings::getGeneticCode
#'
#' @return A data.frame containing values which indicate the degree of codon
#'    usage similarity between the sequences used to construct the
#'    \code{codonFreq} object and those used to construct the reference
#'    database. The reference database entries are in rows and the
#'    \code{codonFreq} sequences are in columns.
#'
#' @examples
#' ## To be added
#'
#' @name codonUsage
NULL
