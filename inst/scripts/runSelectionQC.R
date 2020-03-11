#!/usr/bin/Rscript

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
# This is a command line wrapper for selectionQC
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
	"Performs a selection QC analysis on the output of runScoring.R, yielding informative plots.",
	name="runSelectionQC.R"
)
p <- add_argument(p, "dataDir", help="workspace data directory")
p <- add_argument(p, "--parameters", help="parameter file. Defaults to parameters.json in the data directory.")
p <- add_argument(p, "--sdThreshold", default=0.3, help="Stdev threshold for determination of syn/stop medians. Default 0.3")
p <- add_argument(p, "--logfile", help="log file. Defaults to selectionQC.log in the same directory")
p <- add_argument(p, "--srOverride", help="Manual override to allow singleton replicates. USE WITH EXTREME CAUTION!",flag=TRUE)
args <- parse_args(p)

#ensure datadir ends in "/" and exists
dataDir <- args$dataDir
if (!grepl("/$",dataDir)) {
	dataDir <- paste0(dataDir,"/")
}
if (!dir.exists(dataDir)) {
	#logger cannot initialize without dataDirectory, so just a simple exception here.
	stop("Data folder does not exist!")
}
paramfile <- if (is.na(args$parameters)) paste0(dataDir,"parameters.json") else args$parameters
logfile <- if (is.na(args$logfile)) paste0(dataDir,"selectionQC.log") else args$logfile

#set up logger and shunt it into the error handler
logger <- new.logger(logfile)
registerLogErrorHandler(logger)

#run the actual function
invisible(
	selectionQC(
		dataDir,paramfile,logger,args$sdThreshold,srOverride=args$srOverride
	)
)


