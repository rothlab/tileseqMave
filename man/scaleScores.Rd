% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scaling.R
\name{scaleScores}
\alias{scaleScores}
\title{perform scaling of logPhi into final scores}
\usage{
scaleScores(
  dataDir,
  scoreDir = NA,
  outDir = NA,
  paramFile = paste0(dataDir, "parameters.json"),
  srOverride = FALSE,
  bnOverride = FALSE,
  autoPivot = FALSE,
  useSimpleAverages = FALSE,
  codonQuorum = c("off", "simple", "harsh")
)
}
\arguments{
\item{dataDir}{working data directory}

\item{scoreDir}{input directory for scores, defaults to subdirectory with latest timestamp ending in _scores.}

\item{outDir}{output directory, defaults to name of input directory with _QC tag attached.}

\item{paramFile}{input parameter file. defaults to <dataDir>/parameters.json}

\item{srOverride}{single replicate override}

\item{bnOverride}{bottleneck filter override}
}
\value{
NULL. Results are written to file.
}
\description{
perform scaling of logPhi into final scores
}
