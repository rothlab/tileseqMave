#!/usr/bin/env Rscript

# Copyright (C) 2018  Jochen Weile, Roth Lab
#
# This file is part of tileseqMave.
#
# tileseqMave is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# tileseqMave is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with tileseqMave.  If not, see <https://www.gnu.org/licenses/>.

#####################################################
# This is a command line wrapper for scaling.R
#####################################################


options(
  stringsAsFactors=FALSE,
  ignore.interactive=TRUE,
  future.cmdargs=c()
)

#load libraries
library(tileseqMave)
library(argparser)
library(yogilog)

#process command line arguments
p <- arg_parser(
  "Performs scaling and scoring on the output of calcEnrichment.R.",
  name="scaleScores.R"
)
p <- add_argument(p, "--workspace", help="workspace data directory. Defaults to current working directory")
p <- add_argument(p, "--scores", help="input directory containing the score data. Defaults to subdirectory with latest timestamp ending in _scores")
p <- add_argument(p, "--output", help="output directory. Defaults to name of input directory with _scores tag")
p <- add_argument(p, "--parameters", help="parameter file. Defaults to parameters.json in the data directory.")
p <- add_argument(p, "--logfile", help="log file. Defaults to selectionQC.log in the same directory")
p <- add_argument(p, "--srOverride", help="Manual override to allow singleton replicates. USE WITH EXTREME CAUTION!",flag=TRUE)
p <- add_argument(p, "--bnOverride", help="Manual override to disable bottleneck filtering.",flag=TRUE)
p <- add_argument(p, "--autoPivot", help="Automatically determine scale pivots with requiring definitions in parameter sheet",flag=TRUE)
args <- parse_args(p)

#Workaround for bug in future package, that re-uses command line arguments:
#Override commandArgs function with dummy that returns nothing
commandArgs <- function(trailingOnly=FALSE) {
  character()
}

#ensure datadir ends in "/" and exists
if (is.na(args$workspace)) {
  dataDir <- getwd()
} else {
  dataDir <- args$workspace
}
if (!grepl("/$",dataDir)) {
  dataDir <- paste0(dataDir,"/")
}
if (!dir.exists(dataDir)) {
  #logger cannot initialize without dataDirectory, so just a simple exception here.
  stop("Data folder does not exist!")
}
paramFile <- if (is.na(args$parameters)) paste0(dataDir,"parameters.json") else args$parameters
logfile <- if (is.na(args$logfile)) paste0(dataDir,"scaling.log") else args$logfile

#set up logger and shunt it into the error handler
logger <- new.logger(logfile)
registerLogger(logger)
registerLogErrorHandler(logger)
logVersion()

#run the actual function
invisible(
  scaleScores(
    dataDir, scoreDir=args$scores, outDir=args$output, 
    paramFile=paramFile, srOverride=args$srOverride,
    bnOverride=args$bnOverride,autoPivot=args$autoPivot
  )
)

  
  

