#' @include AllClasses.R
#' @include AllGenerics.R
#' @include methods-codonFreq.R
#' @import reshape2
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
#'    included in the analysis. Default = 500.
#' @param norm Logical, is the data pre-normalised? If FALSE, a normalisation
#'    step will be performed. Default = TRUE.
#' @param n Numeric, the top n-ranked taxa (ordered by MCUFD value, from low to
#'    high) to be returned. Default = NA (return all taxa).
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
.enrichTest <- function(df, n, rank) {
    df[, rank] <- as.character(df[, rank])
    fullCounts <- data.frame(df %>% group_by_(rank) %>% tally())
    subCounts <- data.frame(df[1:n,] %>% group_by_(rank) %>% tally())
    fullTax <- fullCounts[, rank]
    ntax <- length(na.omit(fullTax))
    resDF <- data.frame(
        Taxon = character(ntax),
        foldEnrich = numeric(ntax),
        pvalue = numeric(ntax),
        stringsAsFactors = FALSE
    )
    sapply(1:length(na.omit(fullTax)), function(x) {
        taxon <- na.omit(fullTax)[x]
        taxCountFull <- fullCounts[,2][grepl(taxon, fullCounts[,1])]
        taxCountSub <- subCounts[,2][grepl(taxon, subCounts[,1])]
        if (length(taxCountSub) == 0) taxCountSub <- 0
        otherCountFull <- sum(fullCounts[,2]) - taxCountFull
        otherCountSub <- sum(subCounts[,2]) - taxCountSub
        ratioFull <- taxCountFull/otherCountFull
        ratioSub <- taxCountSub/otherCountSub
        foldEnrich <- ratioSub - ratioFull
        fisherp <- fisher.test(
            matrix(
                c(taxCountFull, taxCountSub, otherCountFull, otherCountSub),
                nrow = 2, ncol = 2
            )
        )$p.value
        resDF[x, 1] <<- taxon
        resDF[x, 2] <<- foldEnrich
        resDF[x, 3] <<- fisherp
    })
    resDF
}


#' @rdname codonCompare
setMethod("mcufd",
    signature(cFobj = "codonFreq"),
    function(cFobj, exclude, minlen, norm) {
        if (isTRUE(norm)) {
            cFdat <- cFobj
        } else {
            cFdat <- normalise(cFobj)
        }
        seqidx <- 1:nrow(cFdat@freq)
        if (length(exclude) > 0) {
            keepcF <- which(!(colnames(cFdat@freq) %in% exclude))
            cFdat <- cFdat[seqidx, keepcF]
            keepRef <- which(!(colnames(gbnorm) %in% exclude))
            gbnorm <- gbnorm[, keepRef]
        }
        if (minlen > 0) {
            codidx <- 1:ncol(cFdat@freq)
            keepSeq <- which(cFdat@ncod > minlen)
            cFdat <- cFdat[keepSeq, codidx]
            gbnorm <- subset(gbnorm, X..Codons >= minlen)
        }
        resList <- vector("list", length = nseq(cFdat))
        names(resList) <- cFdat@seqID
        sapply(1:nrow(cFdat@freq), function(i) {
            cFrow <- cFdat@freq[i,]
            gbnorm[, 7:ncol(gbnorm)] <- abs(
                sweep(gbnorm[, 7:ncol(gbnorm)], 2, cFrow)
            )
            gbnorm$mcufd <- rowMeans(gbnorm[, 7:ncol(gbnorm)], na.rm = TRUE)
            gbnorm$seqid <- cFdat@seqID[i]
            gbnorm <- gbnorm[order(gbnorm$mcufd),]
            rownames(gbnorm) <- NULL
            resList[[i]] <<- gbnorm
        })
        resList
})


#' @rdname codonCompare
setMethod("mcufd_enrich",
    signature(cFres = "list"),
    function(cFres, n, rank, plot, fname, units, width, height, dpi) {
        resList <- vector("list", length(cFres))
        sapply(seq_along(cFres), function(i) {
            resList[[i]] <<- .enrichTest(cFres[[i]], n, rank)
        })
        if (plot == TRUE) {
            resMerge <- dplyr::bind_rows(cFres)
            resMelt <- melt(
                resMerge, id.vars = c("Taxon"), measure.vars = "foldEnrich"
            )
            resMelt$Taxon <- factor(
                resMelt$Taxon, levels = rev(levels(factor(resMelt$Taxon)))
            )
            plt <- ggplot(
                    resMelt, aes(x = Taxon, y = log(value))
                ) +
                geom_point(
                    shape = 21, size = 4, fill = "red", alpha = 0.5,
                    # position = position_dodge(width = 0.4),
                    show.legend = FALSE
                ) +
                theme_classic() +
                coord_flip() +
                geom_hline(
                    yintercept = 0, linetype = "dashed",
                    color = "grey", size = 1
                ) +
                xlab("Taxon") + #ylim(range(resMelt$value) + c(-0.1, 0.1)) +
                ylab(paste0("log(Fold enrichment in top ", n, " taxa)"))
            cat(paste0("Saving plot as ", fname, ".png"), "\n")
            ggsave(
                plt,
                file = paste0(fname, ".png"),
                device = "png",
                units = units,
                width = width,
                height = height,
                dpi = dpi
            )
        }
        resList
    }
)
