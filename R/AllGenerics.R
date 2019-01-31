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
#'#'
#' @export
setGeneric("seqlen", function(object) standardGeneric("seqlen"))



##------------------------------------------------------------------------------
## Normalisation
##------------------------------------------------------------------------------
#' @rdname codonFreq-class
#'
#' @export
setGeneric("normalise", function(object) standardGeneric("normalise"))



##------------------------------------------------------------------------------
## Codon comparisons
##------------------------------------------------------------------------------
#' @rdname codonCompare
#'
#' @export
setGeneric(
    "mcufd",
    function(
        cFobj, exclude = character(length = 0),
        minlen = 600, norm = TRUE
    ) standardGeneric("mcufd")
)



##------------------------------------------------------------------------------
## Enrichment testing
##------------------------------------------------------------------------------
#' MCUFD enrichment testing
#'
#' Compare taxon proportions in the top n hits with those in the full set of
#'    taxa.
#'
#' @param cFres List of data frames containing MCUFD results.
#' @param n Numeric, the number of top-ranked reference taxa to be plotted
#'    per input sequence. Default = 100.
#' @param rank Character, taxonomic rank to be used for categorisation of the
#'    CUFD hit sequences. Options are "Domain", "Kingdom", and "Phylum".
#'    Default = "Phylum".
#' @param plot Logical, should the enrichment results be plotted?
#'    Default = FALSE.
#' @param fname Character, name of the output file (used if plot == TRUE).
#' @param units Numeric, units to be used for defining the plot size.
#'    Options are "in" (default), "cm", and "mm".
#' @param width Numeric, width of the figure (in \code{units}).
#' @param height Numeric, height of the figure (in \code{units}).
#' @param dpi Numeric, resolution of the figure (default = 600).
#'
#' @return A list of data frames containing enrichment results. Each data frame
#'
#' @rdname plots
#'
#' @export
setGeneric(
    "mcufd_enrich",
    function(
        cFres, n = NA_real_, rank = "Phylum",
        plot = FALSE, fname = NA_character_,
        units = "in", width = 10, height = 7, dpi = 600
    ) standardGeneric("mcufd_enrich")
)



##------------------------------------------------------------------------------
## Visualisation of results
##------------------------------------------------------------------------------
#' MCUFD plot
#'
#' Plot the distribution of MCUFD values per taxon
#'
#' @param cFres List of data frames containing MCUFD results.
#' @param type Character, the type of plot to make. Options are "bar" (default)
#'    and "line".
#' @param n Numeric, the number of top-ranked reference taxa to be plotted
#'    per input sequence. Default = 100.
#' @param rank Character, taxonomic rank to be used for categorisation of the
#'    CUFD hit sequences. Options are "Domain", "Kingdom", and "Phylum".
#'    Default = "Phylum".
#' @param fname Character, name of figure generated.
#' @param units Numeric, units to be used for defining the plot size.
#'    Options are "in" (default), "cm", and "mm".
#' @param width Numeric, width of the figure (in \code{units}).
#' @param height Numeric, height of the figure (in \code{units}).
#' @param dpi Numeric, resolution of the figure (default = 600).
#'
#' @return A \code{ggplot} object.
#'
#' @rdname plots
#'
#' @export
setGeneric(
    "mcufd_plot",
    function(
        cFres, type = "bar", n = NA_real_, rank = "Phylum",
        fname = NA_character_, units = "in", width = 10, height = 7, dpi = 600
    ) standardGeneric("mcufd_plot")
)
