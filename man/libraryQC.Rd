% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/libraryQC.R
\name{libraryQC}
\alias{libraryQC}
\title{run library quality control (QC)}
\usage{
libraryQC(dataDir, paramFile = paste0(dataDir, "parameters.json"),
  logger = NULL)
}
\arguments{
\item{dataDir}{working data directory}

\item{paramFile}{input parameter file. defaults to <dataDir>/parameters.json}
}
\value{
NULL. Results are written to file.
}
\description{
run library quality control (QC)
}