% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/libraryQC.R
\name{coverageSubmap}
\alias{coverageSubmap}
\title{Helper function to draw a subsection of the coverage map}
\usage{
coverageSubmap(
  startPos,
  endPos,
  aaMarginal,
  seps = NULL,
  thresholds = c(1e-06, 5e-05)
)
}
\arguments{
\item{startPos}{start position for subsection of coverage map}

\item{endPos}{end position for subsection of coverage map}

\item{aaMarginal}{marginal frequency table at amino acid-level}

\item{seps}{tile separator lines}

\item{thresholds}{heatmap colour transition thresholds}
}
\description{
Helper function to draw a subsection of the coverage map
}
