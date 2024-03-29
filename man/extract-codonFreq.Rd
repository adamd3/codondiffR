% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods-codonFreq.R
\name{[,codonFreq,ANY,ANY,ANY-method}
\alias{[,codonFreq,ANY,ANY,ANY-method}
\alias{[[,codonFreq-method}
\title{Methods for subsetting a codonFreq object}
\usage{
\S4method{[}{codonFreq,ANY,ANY,ANY}(x, i, j)

\S4method{[[}{codonFreq}(x, i, j)
}
\arguments{
\item{x}{
    object from which to extract element(s) or in which to replace element(s).
  }

\item{i}{
    indices specifying elements to extract or replace.  Indices are
    \code{numeric} or \code{character} vectors or empty (missing) or
    \code{NULL}.  Numeric values are coerced to integer as by
    \code{\link[base]{as.integer}} (and hence truncated towards zero).
    Character vectors will be matched to the \code{\link[base]{names}} of the
    object (or for matrices/arrays, the \code{\link[base]{dimnames}}):
    see \sQuote{Character indices} below for further details.

    For \code{[}-indexing only: \code{i}, \code{j}, \code{\dots} can be
    logical vectors, indicating elements/slices to select.  Such vectors
    are recycled if necessary to match the corresponding extent.
    \code{i}, \code{j}, \code{\dots} can also be negative integers,
    indicating elements/slices to leave out of the selection.

    When indexing arrays by \code{[} a single argument \code{i} can be a
    matrix with as many columns as there are dimensions of \code{x}; the
    result is then a vector with elements corresponding to the sets of
    indices in each row of \code{i}.

    An index value of \code{NULL} is treated as if it were \code{integer(0)}.
  }

\item{j}{
    indices specifying elements to extract or replace.  Indices are
    \code{numeric} or \code{character} vectors or empty (missing) or
    \code{NULL}.  Numeric values are coerced to integer as by
    \code{\link[base]{as.integer}} (and hence truncated towards zero).
    Character vectors will be matched to the \code{\link[base]{names}} of the
    object (or for matrices/arrays, the \code{\link[base]{dimnames}}):
    see \sQuote{Character indices} below for further details.

    For \code{[}-indexing only: \code{i}, \code{j}, \code{\dots} can be
    logical vectors, indicating elements/slices to select.  Such vectors
    are recycled if necessary to match the corresponding extent.
    \code{i}, \code{j}, \code{\dots} can also be negative integers,
    indicating elements/slices to leave out of the selection.

    When indexing arrays by \code{[} a single argument \code{i} can be a
    matrix with as many columns as there are dimensions of \code{x}; the
    result is then a vector with elements corresponding to the sets of
    indices in each row of \code{i}.

    An index value of \code{NULL} is treated as if it were \code{integer(0)}.
  }
}
\value{
A subset of the sequences in the \code{codonFreq} object.
}
\description{
Methods for subsetting a codonFreq object
}
\seealso{
[base::Extract]
}
