#' @include AllClasses.R
#' @include AllGenerics.R
#' @include methods-codonFreq.R
NULL



##------------------------------------------------------------------------------
## Documentation for codon comparison methods
##------------------------------------------------------------------------------
#' Calculate measures of codon usage similarity.
#'
#' Compare codon usage between sequences in a \code{codonFreq} object with those
#'    in a reference database.
#' Mean codon usage frequency difference (\code{MCUFD}) is calculated as
#'    described in \href{https://goo.gl/7sGZeE}{Stedman et al. (2013)}.
#'
#' @param cFobj An object of class \code{codonFreq}.
#' @param exclude A character vector of codons to be excluded from comparisons.
#' @param minlen Numeric, the minimum length of sequence (in codons) to be
#'    included in the analysis. Default = 600.
#' @param norm Logical, should the data be normalised? Default = TRUE.
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
#' @name codonCompare
NULL



##------------------------------------------------------------------------------
## Mean codon usage frequency difference (MCUFD)
##------------------------------------------------------------------------------
#' @rdname codonCompare
setMethod("MCUFD",
    signature(cFobj = "codonFreq"),
    function(cFobj, exclude, minlen, norm) {
        if (isTRUE(norm)) {
            GBdat <- gbnorm
            cFdat <- normalise(cFobj)
        } else {
            GBdat <- GenbankCUTD
            cFdat <- cFobj
        }
        seqidx <- 1:nrow(cFdat@freq)
        if (length(exclude) > 0) {
            ## Subset the codonFreq object
            keepcF <- which(!(colnames(cFdat@freq) %in% exclude))
            cFdat <- cFdat[seqidx, keepcF]
            ## Subset the Reference taxa data
            keepRef <- which(!(colnames(GBdat) %in% exclude))
            GBdat <- GBdat[, keepRef]
        }
        if (minlen > 0) {
            codidx <- 1:ncol(cFdat@freq)
            keepSeq <- which(cFdat@seqlen > minlen)
            cFdat <- cFdat[keepSeq, codidx]
            GBdat <- subset(GBdat, X..Codons >= (minlen/3))
        }
        res <- list()
        res <- c(
            res,
            apply(cFdat@freq, 1, function(cFrow) {
                GBdat[, 6:ncol(GBdat)] <- abs(
                    sweep(GBdat[, 6:ncol(GBdat)], 2, cFrow)
                )
                GBdat$mcufd <- rowMeans(GBdat[, 6:ncol(GBdat)], na.rm = TRUE)
                GBdat <- GBdat[order(GBdat$mcufd),]
                return(GBdat)
            })
        )
        names(res) <- cFdat@seqID
        return(res)
})
