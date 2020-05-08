#' @include AllClasses.R
#' @include AllGenerics.R
#' @include methods-codonFreq.R
NULL


##------------------------------------------------------------------------------
## Principal component analysis (PCA) methods
##------------------------------------------------------------------------------
#' PCA of codon usage by taxon
#'
#' Perform PCA on codon usage in a reference
#'    database and sequences of unknown taxonomic affinity.
#'
#' @param exclude A character vector of codons to be excluded from comparisons.
#' @param minlen Numeric, the minimum length of sequence (in codons) to be
#'    included in the analysis. Default = 500.
#' @param norm Logical, should the codon abundances be normalised? If TRUE,
#'    codon abundances will be converted to codon bias scores, such that the sum
#'    of scores for each amino acid sum to 1. Default = FALSE.
#' @param rank Character, taxonomic rank to be used for categorisation. Options
#'    are "Domain", "Kingdom", and "Phylum". Default = "Phylum".
#' @param corCut Numeric, correlation cutoff used for dropping codons that
#'    exhibit collinearity. Default = 0.9.
#' @param includeTax Character, optional vector of taxa (at the specified rank)
#'    to be included in the PCA. If not supplied, all taxa will be used.
#'
#' @return An object of class \code{prcomp}.
#'
#' @examples
#'     PCA_tmp <- PCA(exclude = exclCod, rank = "Phylum", minlen = 600)
#'     class(PCA_tmp)
#'
#' @name PCA
#' @rdname PCA
#'
#' @export
setMethod("PCA", signature = character(),
    function(exclude, minlen, norm, rank, corCut, includeTax) {
        if (length(exclude) > 0) {
            keepRef <- which(!(colnames(gbnorm) %in% exclude))
            gbnorm <- gbnorm[, keepRef]
        }
        if (minlen > 0) gbnorm <- subset(gbnorm, X..Codons >= minlen)
        rmcols <- dplyr::setdiff(colnames(gbnorm)[1:6], rank)
        gbnorm <- gbnorm[, !(names(gbnorm) %in% rmcols)]
        ## remove NAs from taxon column
        completeVec <- complete.cases(gbnorm[, 1])
        gbnorm <- gbnorm[completeVec, ]
        ## remove variables with zero variance
        colVars <- apply(gbnorm[2:ncol(gbnorm)], 2, function (x) {
            var(x, na.rm = TRUE)
        })
        rmvars <- names(colVars[colVars == 0])
        ## drop collinear variables
        if (corCut < 1) {
            cormat <- cor(
                gbnorm[, 2:ncol(gbnorm)],
                gbnorm[, 2:ncol(gbnorm)],
                use = "pairwise.complete.obs"
            )
            rmcor <- findCorrelation(
                cormat,
                cutoff = corCut,
                verbose = FALSE,
                names = TRUE
            )
            rmvars <- c(rmvars, rmcor)
        }
        gbnorm <- gbnorm[, -which(names(gbnorm) %in% rmvars)]
        ## remove taxa in includeTax, if supplied
        if (!is.null(includeTax)) {
            if (all(includeTax %in% gbnorm[[eval(rank)]])) {
                gbnorm <- gbnorm[gbnorm[[eval(rank)]] %in% includeTax ,]
            } else {
                stop(
                    "All values in `includeTax` must be in the specified rank"
                )
            }

        }
        gbnorm[,1] <- as.factor(gbnorm[,1])
        m1 <- prcomp(gbnorm[,2:ncol(gbnorm)], center = TRUE, scale. = FALSE)
        # ggbiplot(m1, groups = gbnorm[,1], ellipse = TRUE, circle = TRUE)
        m1
})




##------------------------------------------------------------------------------
## Predict method and plots
##------------------------------------------------------------------------------
#' PCA plots
#'
#' Predict taxonomic classifications for sequences in a \code{codonFreq})
#' object, using linear discriminants from an \code{PCA}) model. Plot
#' disciminants and predictions.
#'
#' @param cFobj An object of class \code{codonFreq}.
#' @param pcaObj Object of class \code{prcomp}.
#' @param rank Character, taxonomic rank to be used for categorisation. Options
#'    are "Domain", "Kingdom", and "Phylum". Default = "Phylum".
#' @param minlen Numeric, the minimum length of sequence (in codons) to be
#'    included in the analysis. Default = 500.
#' @param fname Character, name of figure generated.
#' @param units Numeric, units to be used for defining the plot size.
#'    Options are "in" (default), "cm", and "mm".
#' @param width Numeric, width of the figure (in \code{units}).
#' @param height Numeric, height of the figure (in \code{units}).
#' @param dpi Numeric, resolution of the figure (default = 600).
#' @param norm Logical, should the codon abundances be normalised? If TRUE,
#'    codon abundances will be converted to codon bias scores, such that the sum
#'    of scores for each amino acid sum to 1. Default = FALSE.
#' @param save Logical, should the enrichment results plot be saved to a file?
#'    Default = FALSE.
#' @param identifier Character, optional group label to be assigned to sequences
#'    in the \code{codonFreq} object. If not supplied, each sequence will be
#'    labelled individually on the plot. Maximum number of unique identifiers
#'    = 6.
#' @param includeTax Character, optional vector of taxa (at the specified rank)
#'    to be included in the PCA. If not supplied, all taxa will be used.
#' @param colours Integer vector, containing values between 1 and 9, which
#'    specifies the colours to be used from the `Set1` palette in the
#'    `RColorBrewer` package.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#'    predPCA <- predict_PCA(
#'        tmp2norm, PCA_tmp2, rank = "Phylum", save = TRUE,
#'        minlen = 600, fname = "PCA_tmp2", height = 5, width = 7
#'    )
#'
#' @name predict_PCA
#' @rdname predict_PCA
#'
#' @export
setMethod("predict_PCA",
    signature(cFobj = "codonFreq"),
    function(
        cFobj, pcaObj, rank, minlen, fname, units, width,
        height, dpi, norm, save, identifier, colours
    ) {
        if (isTRUE(norm)) cFobj <- normalise(cFobj)
        if (minlen > 0) {
            codidx <- 1:ncol(cFobj@freq)
            keepSeq <- which(cFobj@ncod > minlen)
            cFobj <- cFobj[keepSeq, codidx]
        }
        ## codons included in PCA
        codonsInc <- rownames(factoextra::get_pca_var(pcaObj)$coord)
        keepcF <- which(colnames(cFobj@freq) %in% codonsInc)
        cFobj <- cFobj[, keepcF]
        rmcols <- dplyr::setdiff(colnames(gbnorm)[1:6], rank)
        gbnorm <- gbnorm[,!(colnames(gbnorm) %in% rmcols)]
        gbnorm <- cbind(
            gbnorm[,1], gbnorm[colnames(gbnorm) %in% codonsInc]
        )
        cFdat <- cbind(cFobj@seqID, as.data.frame(cFobj@freq))
        colnames(cFdat)[1] <- colnames(gbnorm)[1] <- "Taxon"
        if (is.null(identifier)) {
            cFdat[,1] <- as.factor(str_sub(cFdat[,1], 1, 20))
        } else {
            if (length(unique(identifier)) > 6) {
                stop("Maximum number of unique values in `identifier` is 6.")
            } else {
                cFdat[,1] <- as.factor(identifier)
            }
        }
        ## remove taxa in includeTax, if supplied
        if (!is.null(includeTax)) {
            if (all(includeTax %in% gbnorm$Taxon)) {
                gbnorm <- gbnorm[gbnorm$Taxon %in% includeTax ,]
            } else {
                stop(
                    "All values in `includeTax` must be in the specified rank"
                )
            }

        }
        gbnorm[,1] <- as.factor(gbnorm[,1])
        pred <- predict(pcaObj, newdata = cFdat[,2:ncol(cFdat)])
        cc1 <- 12
        brewer_pallette <- c(brewer.pal(9, "Set1"), "black")
        brewer_pallette2 <- brewer_pallette[c(1:4,7:9)] ##exclude yellow
        brewer_greys <- brewer.pal(9, "Greys")[3:8]
        ngp <- length(includeTax)
        if (!is.null(colours)) {
            colvec <- c(colours)
        } else {
            colvec <- c(
                brewer_pallette2[1:ngp],
                brewer_greys[1:length(unique(identifier))+1]
            )
        }
        p1 <- ggbiplot(
                pcaObj, groups = gbnorm[,1], ellipse = TRUE, circle = TRUE,
                var.axes = FALSE, alpha = 0
                ) +
                scale_colour_manual("Group", values = colvec, guide = FALSE) +
                scale_fill_manual("Group", values = colvec) +
                scale_shape_manual(values = c(rep(16, ngp), 17)) +
                theme_bw() +
                theme(
                  text = element_text(size = cc1),
                  axis.text.x = element_text(colour = "black", size=cc1),
                  axis.text.y = element_text(colour = "black", size=cc1)
              ) +
              geom_point(
                  size = 3, shape = 21, aes(fill = groups)
              )
        ## add predictions
        p1$data <- rbind(p1$data,
            data.frame(
                xvar = pred[, "PC1"],
                yvar = pred[, "PC2"],
                groups = cFdat[,1]
            )
        )
        if (isTRUE(save)) {
            ggsave(
                p1,
                file = paste0(fname, ".png"),
                device = "png",
                units = units,
                width = width,
                height = height,
                dpi = 300
            )
        }
        p1
    }
)
