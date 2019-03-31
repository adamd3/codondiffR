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
#'
#' @export
setGeneric("seqlen", function(object) standardGeneric("seqlen"))



##------------------------------------------------------------------------------
## Plots
##------------------------------------------------------------------------------
#' @rdname codonFreq-class
#'
#' @export
setGeneric(
    "codonPlot",
    function(
        object, fname = NA_character_, units = "in", width = 10,
        height = 7, dpi = 600, groups = NULL, ptype = "boxplot",
        order = "median", colour = 1, suppress_x_txt = FALSE,
        suppress_y_title = FALSE, label = NULL, highlight = NULL,
        save = FALSE
    ) standardGeneric("codonPlot")
)

setGeneric(
    "biasPlot",
    function(
        object, fname = NA_character_, units = "in", width = 10,
        height = 7, dpi = 600, groups = NULL, aa = NULL, norm = FALSE,
        label = NULL, colours = NULL, suppress_y_txt = FALSE,
        suppress_y_title = FALSE, legend = TRUE, ylim = NULL, save = FALSE
    ) standardGeneric("biasPlot")
)

setGeneric(
    "gcPlot",
    function(
        object, fname = NA_character_, units = "in", width = 10,
        height = 7, dpi = 600, groups = NULL, aa = NULL, norm = FALSE,
        label = NULL, colours = NULL, suppress_y_txt = FALSE,
        suppress_y_title = FALSE, legend = TRUE, xlim = NULL,
        outtab = NULL, save = FALSE
    ) standardGeneric("gcPlot")
)


##------------------------------------------------------------------------------
## Normalisation
##------------------------------------------------------------------------------
#' @rdname codonFreq-class
#'
#' @export
setGeneric("normalise", function(object) standardGeneric("normalise"))



##------------------------------------------------------------------------------
## MCUFD
##------------------------------------------------------------------------------
#' @rdname MCUFD
#'
#' @export
setGeneric(
    "MCUFD",
    function(
        cFobj, exclude = character(length = 0),
        minlen = 600, norm = FALSE
    ) standardGeneric("MCUFD")
)


#' @rdname MCUFD
#'
#' @export
setGeneric(
    "MCUFD_enrich",
    function(
        cFres, n = NA_real_, rank = "Phylum",
        save = FALSE, pthresh = NA_real_, fname = NA_character_,
        units = "in", width = 10, height = 7, dpi = 600,
        ptype = "heatmap", outtab = NULL
    ) standardGeneric("MCUFD_enrich")
)

#' @rdname MCUFD
#'
#' @export
setGeneric(
    "MCUFD_plot",
    function(
        cFres, type = "bar", n = NA_real_, rank = "Phylum",
        fname = NA_character_, units = "in", width = 10, height = 7, dpi = 600,
        save = FALSE
    ) standardGeneric("MCUFD_plot")
)



##------------------------------------------------------------------------------
## PCA
##------------------------------------------------------------------------------
#' @rdname PCA
#'
#' @export
setGeneric(
    "PCA",
    function(
        exclude = character(length = 0), minlen = 600,
        norm = FALSE, rank = "Phylum", corCut = 0.9, includeTax = NULL
    ) standardGeneric("PCA")
)


#' @rdname predict_PCA
#'
#' @export
setGeneric(
    "predict_PCA",
    function(
        cFobj, pcaObj, rank = "Phylum",
        minlen = 600, fname = NA_character_, units = "in",
        width = 10, height = 7, dpi = 600, norm = FALSE, save = FALSE,
        identifier = NULL, includeTax = NULL
    ) standardGeneric("predict_PCA")
)



##------------------------------------------------------------------------------
## LDA
##------------------------------------------------------------------------------
#' @rdname LDA
#'
#' @export
setGeneric(
    "LDA",
    function(
        exclude = character(length = 0), minlen = 600,
        trans = FALSE, propTrain = 1,
        rank = "Phylum", corCut = 0.9
    ) standardGeneric("LDA")
)

#' @rdname LDA
#'
#' @export
setGeneric(
    "bootstrap_LDA",
    function(
        rep = 100, trans = FALSE, propTrain = 0.8, rank = "Phylum",
        exclude = character(length = 0), minlen = 600,
        corCut = 0.9, norm = FALSE
    )
    standardGeneric("bootstrap_LDA")
)


#' @rdname LDA
#'
#' @export
setGeneric(
    "predict_LDA",
    function(
        cFobj, ldaObj, rank = "Phylum",
        minlen = 600, fname = NA_character_, units = "in",
        width = 10, height = 7, dpi = 600, norm = FALSE, plot = FALSE,
        identifier = NA_character_
    ) standardGeneric("predict_LDA")
)
