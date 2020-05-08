#' @include AllClasses.R
#' @include AllGenerics.R
#' @include methods-codonFreq.R
NULL


##------------------------------------------------------------------------------
## Functions
##------------------------------------------------------------------------------
## Box-Cox transformation
.boxcox <- function(in_df, exclude_cols){
    transform <- preProcess(
        in_df[, -exclude_cols], c("BoxCox", "center", "scale")
    )
    out_df <- data.frame(trans = predict(transform, in_df[, -exclude_cols]))
    out_df <- cbind(in_df[, exclude_cols], out_df)
    colnames(out_df) <- colnames(in_df)
    out_df
}



##------------------------------------------------------------------------------
## Discriminant analysis methods
##------------------------------------------------------------------------------
#' Discriminant analysis of codon usage by taxon
#'
#' Perform linear discriminant analysis (LDA) on codon usage in a reference
#'    database and use it to classify sequences of unknown taxonomic affinity.
#'    Data can be optionally scaled using the
#'    \href{https://goo.gl/qzN63P}{Box-Cox power transformation}.
#'
#' @param exclude A character vector of codons to be excluded from comparisons.
#' @param minlen Numeric, the minimum length of sequence (in codons) to be
#'    included in the analysis. Default = 500.
#' @param trans Logical; if true, a Box-Cox transformation will be applied to
#'    the data. Default = FALSE.
#' @param propTrain Numeric, proportion of the reference database to use for
#'    the LDA training set (must be in the range [0,1]). Default = 1.
#' @param rank Character, taxonomic rank to be used for categorisation. Options
#'    are "Domain", "Kingdom", and "Phylum". Default = "Phylum".
#' @param corCut Numeric, correlation cutoff used for dropping codons that
#'    exhibit collinearity. Default = 0.9.
#'
#' @return An object of class \code{lda}, with the following components:
#' @return "prior": the prior probabilities used (determined from training set).
#' @return "means": the group means.
#' @return "scaling": a matrix which transforms observations to discriminant
#'     functions, normalized so that within groups covariance matrix is
#'     spherical.
#' @return "svd": the singular values, which give the ratio of the between- and
#'     within-group standard deviations on the linear discriminant variables.
#'     Their squares are the canonical F-statistics.
#' @return "N": The number of observations used.
#' @return "call": The (matched) function call to MASS:lda.
#'
#' @examples
#'     LDA_tmp <- LDA(
#'         exclude = exclCod, rank = "Phylum", trans = FALSE,
#'         propTrain = 1, corCut = 0.95, minlen = 600
#'     )
#'     names(LDA_tmp)
#'
#' @name LDA
#' @rdname LDA
#'
#' @export
setMethod("LDA", "data.frame",
    function(exclude, minlen, trans, propTrain, rank, corCut) {
        if (length(exclude) > 0) {
            keepRef <- which(!(colnames(gbnorm) %in% exclude))
            gbnorm <- gbnorm[, keepRef]
        }
        if (minlen > 0) {
            gbnorm <- subset(gbnorm, X..Codons >= minlen)
        }
        if (trans == TRUE) gbnorm <- .boxcox(gbnorm, 1:6)
        if ((0 > propTrain) | (1 < propTrain)) {
            stop("propTrain must be in the range [0,1]")
        }
        rmcols <- setdiff(colnames(gbnorm)[1:6], rank)
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
        gbnorm[,1] <- as.factor(gbnorm[,1])
        ntrain <- round((nrow(gbnorm) * propTrain), digits = 0)
        ntest <- nrow(gbnorm) - ntrain
        train <- sample(1:nrow(gbnorm), ntrain)
        f <- as.formula(substitute(x ~ ., list(x = as.name(rank))))
        m1 <- lda(f, gbnorm, na.action = "na.omit", subset = train)
        m1
})


#' Bootstrap method for linear discriminant analysis
#'
#' Perform random sampling using different subsets of the reference database
#'    to assess the impact on model accuracy.
#'
#' @param rep Numeric, number of bootstrap replicates to perform (default =
#'    100).
#' @param exclude A character vector of codons to be excluded from comparisons.
#' @param minlen Numeric, the minimum length of sequence (in codons) to be
#'    included in the analysis. Default = 500.
#' @param norm Logical, should the codon abundances be normalised? If TRUE,
#'    codon abundances will be converted to codon bias scores, such that the sum
#'    of scores for each amino acid sum to 1. Default = FALSE.
#' @param trans Logical; if true, a Box-Cox transformation will be applied to
#'    the data. Default = FALSE.
#' @param propTrain Numeric, proportion of the reference database to use for
#'    the LDA training set (must be in the range [0,1]). Default = 0.8.
#' @param rank Character, taxonomic rank to be used for categorisation. Options
#'    are "Domain", "Kingdom", and "Phylum". Default = "Phylum".
#' @param corCut Numeric, correlation cutoff used for dropping codons that
#'    exhibit collinearity. Default = 0.9.
#'
#' @return A named list with the following components:
#' @return "codons": the codons used for model construction, aftering filtering
#'    criteria have been applied
#' @return "taxa": taxa used for model construction/testing
#' @return "accuracy": the proportion of accurate classifications using the test
#'    subset of the data
#'
#' @examples
#'    exclCod <- c("ATT", "TGT")
#'    boot <- bootstrap_LDA(
#'        rep = 100, propTrain = 0.8, trans = FALSE, rank = "Phylum",
#'        exclude = exclCod, minlen = 600, corCut = 0.95
#'    )
#'    names(boot)
#'    mean(boot$Accuracy)
#'
#' @name LDA
#' @rdname LDA
#'
#' @export
setMethod("bootstrap_LDA",
    signature = c(rep = "numeric"),
    function(rep, trans, propTrain, rank, exclude, minlen, corCut, norm) {
        resList <- list(
            codons = character(),
            taxa = character(),
            accuracy = numeric(rep)
        )
        if (length(exclude) > 0) {
            keepRef <- which(!(colnames(gbnorm) %in% exclude))
            gbnorm <- gbnorm[, keepRef]
        }
        if (minlen > 0) gbnorm <- subset(gbnorm, X..Codons >= minlen)
        if (trans == TRUE) gbnorm <- .boxcox(gbnorm, 1:6)
        if ((0 > propTrain) | (1 < propTrain)) {
            stop("propTrain must be in the range [0,1]")
        }
        rmcols <- setdiff(colnames(gbnorm)[1:6], rank)
        gbnorm <- gbnorm[ , !(names(gbnorm) %in% rmcols)]
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
        resList$codons <- colnames(gbnorm)[2:ncol(gbnorm)]
        resList$taxa <- unique(gbnorm[,1])
        gbnorm[,1] <- as.factor(gbnorm[,1])
        ntrain <- round((nrow(gbnorm) * propTrain), digits = 0)
        ntest <- nrow(gbnorm) - ntrain
        sapply(1:rep, function(k) {
            assign("gbnorm", gbnorm, envir = parent.frame(n=2))
            train <- sample(1:nrow(gbnorm), ntrain)
            f <- as.formula(substitute(x ~ ., list(x = as.name(rank))))
            m1 <- lda(f, gbnorm, na.action = "na.omit", subset = train)
            predAll <- predict(m1, gbnorm[-train, ])
            tablin <- table(gbnorm[-train, 1], predAll$class)
            acc <- sum(diag(tablin))/sum(tablin)
            resList$accuracy[k] <<- acc
        })
        resList
})



##------------------------------------------------------------------------------
## Predict method and plots
##------------------------------------------------------------------------------
#' LDA plots
#'
#' Predict taxonomic classifications for sequences in a \code{codonFreq})
#' object, using linear discriminants from an \code{lda}) model. Plot
#' disciminants and predictions.
#'
#' @param cFobj An object of class \code{codonFreq}.
#' @param ldaObj Object of class \code{lda} - produced using the LDA() function.
#' @param rank Character, taxonomic rank to be used for categorisation of the
#'    CUFD hit sequences. Options are "Domain", "Kingdom", and "Phylum".
#'    Default = "Phylum".
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
#' @param plot Logical, should the enrichment results be plotted?
#'    Default = FALSE.
#' @param identifier Character, optional group label to be assigned to sequences
#'    in the \code{codonFreq} object. If not supplied, each sequence will be
#'    labelled individually on the plot.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#'    predLDA <- predict_LDA(
#'        tmp2norm, LDA_tmp2, rank = "Phylum", plot = TRUE,
#'        minlen = 600, fname = "lda_tmp2", height = 5, width = 7
#'    )
#'    head(sort(table(predLDA$class), decreasing = TRUE))
#'
#' @rdname plots
#'
#' @export
setMethod("predict_LDA",
    signature(cFobj = "codonFreq"),
    function(
        cFobj, ldaObj, rank, minlen,
        fname, units, width, height, dpi, norm, plot, identifier
    ) {
        if (isTRUE(norm)) cFobj <- normalise(cFobj)
        if (minlen > 0) {
            codidx <- 1:ncol(cFobj@freq)
            keepSeq <- which(cFobj@ncod > minlen)
            cFobj <- cFobj[keepSeq, codidx]
        }
        codonsInc <- colnames(ldaObj$means) ## codons included in LDA model
        keepcF <- which(colnames(cFobj@freq) %in% codonsInc)
        cFobj <- cFobj[, keepcF]
        rmcols <- setdiff(colnames(gbnorm)[1:6], rank)
        gbnorm <- gbnorm[,!(colnames(gbnorm) %in% rmcols)]
        gbnorm <- cbind(
            gbnorm[,1], gbnorm[colnames(gbnorm) %in% codonsInc]
        )
        cFdat <- cbind(cFobj@seqID, as.data.frame(cFobj@freq))
        colnames(cFdat)[1] <- colnames(gbnorm)[1] <- "Taxon"
        if (is.na(identifier)) {
            cFdat[,1] <-  as.factor(str_sub(cFdat[,1], 1, 20))
        } else {
            cFdat[,1] <- as.factor(identifier)
        }
        allDat <- rbind(gbnorm, cFdat)
        allDat[,1] <- as.factor(allDat[,1])
        ## remove NAs from taxon column
        completeVec <- complete.cases(allDat[, 1])
        allDat <- allDat[completeVec, ]
        ## get proportion of inter-group variance explained by LDs
        prop_lda <- ldaObj$svd^2/sum(ldaObj$svd^2)
        predcF <- predict(object = ldaObj, newdata = cFdat[, 2:ncol(cFdat)])
        if (isTRUE(plot)) {
            predAll <- predict(
                object = ldaObj, newdata = allDat[, 2:ncol(allDat)]
            )
            plot_df <- data.frame(
                Taxon = as.factor(allDat[,1]), lda = predAll$x
            )
            ## plot the n most abundant predicted categories
            if (rank == "Kingdom") nplot <- 3 else nplot <- 9
            abClass <- names(
                sort(table(predcF$class), decreasing = TRUE)
            )[1:nplot]
            abClass <- c(abClass, levels(cFdat[,1])[1])
            plot_df <- plot_df[plot_df$Taxon %in% abClass,]
            cc1 <- 12
            brewer_pallette <- c(brewer.pal(nplot, "Set1"), "black")
            p1 <- ggplot(plot_df) +
                geom_point(
                    aes_string(
                        lda.LD1, lda.LD2,
                        colour = Taxon, size = Taxon,
                        shape = Taxon, alpha = Taxon
                    )
                ) +
                scale_colour_manual(values = brewer_pallette) +
                scale_size_manual(values = c(rep(2, nplot), 3)) +
                scale_shape_manual(values = c(rep(16, nplot), 17)) +
                scale_alpha_manual(values = c(rep(0.7, nplot), 1)) +
                labs(
                    x = paste("LD1 (", percent(prop_lda[1]), ")", sep=""),
                    y = paste("LD2 (", percent(prop_lda[2]), ")", sep="")
                ) +
                theme_bw() +
                theme(
                  text = element_text(size = cc1),
                  axis.text.x = element_text(colour = "black", size=cc1),
                  axis.text.y = element_text(colour = "black", size=cc1)
                )
            ggsave(
                p1,
                file = paste0(fname, ".png"),
                device = "png",
                units = units,
                width = width,
                height = height,
                dpi = dpi
            )
        }
        predcF
    }
)
