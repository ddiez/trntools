% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot_trn}
\alias{plot_trn}
\alias{plot_trn.data.frame}
\alias{plot_trn.tbl_graph}
\title{Plot transcriptional regulatory networks.}
\usage{
plot_trn(x, ...)

\method{plot_trn}{data.frame}(
  x,
  node.size = 10,
  label.size = 3,
  regulator.color = "palevioletred1",
  target.color = "turquoise3",
  color.by = NULL,
  color.low = "blue",
  color.mid = "white",
  color.high = "red"
)

\method{plot_trn}{tbl_graph}(x, layout = "nicely", ...)
}
\arguments{
\item{x}{a network object with attributes.}
}
\description{
Plot transcriptional regulatory networks.
}
