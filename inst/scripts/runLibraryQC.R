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
# This is a command line wrapper for libraryQC
#####################################################

options(
  stringsAsFactors=FALSE,
  ignore.interactive=TRUE
)

#load libraries
library(tileseqMave)
library(argparser)
library(yogilog)

#process command line arguments
p <- arg_parser(
  "Performs a library QC analysis on the output of joinCounts.R, yielding informative plots.",
  name="runLibraryQC.R"
)
p <- add_argument(p, "--workspace", help="workspace data directory. Defaults to current working directory")
p <- add_argument(p, "--input", help="input directory containing the count data. Defaults to subdirectory with latest timestamp ending in _mut_count")
p <- add_argument(p, "--output", help="output directory. Defaults to name of input directory with _QC tag")
p <- add_argument(p, "--parameters", help="parameter file. Defaults to parameters.json in the data directory.")
p <- add_argument(p, "--logfile", help="log file. Defaults to libraryQC.log in the same directory")
p <- add_argument(p, "--cores", default=6L, help="number of CPU cores to use in parallel for multi-threading")
p <- add_argument(p, "--srOverride", help="Manual override to allow singleton replicates. USE WITH EXTREME CAUTION!",flag=TRUE)
p <- add_argument(p, "--allConditions", help="Run on all conditions instead of just nonselect.",flag=TRUE)
p <- add_argument(p, "--depthFilterOverride", help="Disable minimum effective depth filter.",flag=TRUE)
#Workaround: Setting the default value to String type, to avoid bug. Validate manually later
p <- add_argument(p, "--wmThreshold", default="5e-5", help="Define the marginal frequency threshold for well-measuredness.")
p <- add_argument(p, "--silent", help="Turn off message printing to stdout",flag=TRUE)
args <- parse_args(p)

#Manully validate wmThreshold as part of workaround
args$wmThreshold <- as.numeric(args$wmThreshold)
if (is.na(args$wmThreshold)) {
  stop("argument --wmThreshold must be numeric!")
}

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
paramfile <- if (is.na(args$parameters)) paste0(dataDir,"parameters.json") else args$parameters
logfile <- if (is.na(args$logfile)) paste0(dataDir,"libraryQC.log") else args$logfile
mc.cores <- if (is.na(args$cores)) 6 else args$cores

#set up logger and shunt it into the error handler
logger <- new.logger(logfile,stdout=!args$silent)
registerLogger(logger)
registerLogErrorHandler(logger)
logVersion()

#run the actual function
invisible(
  libraryQC(dataDir,inDir=args$input,outDir=args$output,paramfile,mc.cores,
            srOverride=args$srOverride,wmThreshold=args$wmThreshold,
            allCondOverride=args$allConditions,depthFilterOverride=args$depthFilterOverride
  )
)

