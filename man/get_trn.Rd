% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools.R
\name{get_trn}
\alias{get_trn}
\title{Get a TRN at a given FDR.}
\usage{
get_trn(x, fdr = 0.01, method = "bonferroni", cutoff = NULL)
}
\arguments{
\item{x}{a TRN model.}

\item{fdr}{FDR level.}

\item{method}{method for p.adjust().}

\item{cutoff}{magnitude cutoff to filter coefficients.}
}
\description{
Get a TRN at a given FDR.
}
