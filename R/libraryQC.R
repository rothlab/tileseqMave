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

#' run library quality control (QC)
#' 
#' @param dataDir working data directory
#' @param paramFile input parameter file. defaults to <dataDir>/parameters.json
#' @return NULL. Results are written to file.
#' @export
libraryQC <- function(dataDir,paramFile=paste0(dataDir,"parameters.json"),logger=NULL) {

	op <- options(stringsAsFactors=FALSE)

	library(yogitools)
	library(hgvsParseR)

	if (!is.null(logger)) {
		stopifnot(inherits(logger,"yogilogger"))
	}
	logInfo <- function(...) {
		if (!is.null(logger)) {
			logger$info(...)
		} else {
			do.call(cat,c(list(...),"\n"))
		}
	}
	logWarn <- function(...) {
		if (!is.null(logger)) {
			logger$warning(...)
		} else {
			do.call(cat,c("Warning:",list(...),"\n"))
		}
	}

	#make sure data and out dir exist and ends with a "/"
	if (!grepl("/$",dataDir)) {
		dataDir <- paste0(dataDir,"/")
	}
	if (!dir.exists(dataDir)) {
		#we don't use the logger here, assuming that whichever script wraps our function
		#catches the exception and writes to the logger (or uses the log error handler)
		stop("Data folder does not exist!")
	}
	if (!canRead(paramFile)) {
		stop("Unable to read parameter file!")
	}


	logInfo("Reading parameters")
	params <- parseParameters(paramFile)

	#identify nonselect conditions
	nsConditions <- with(as.data.frame(params$conditions$definitions),{
		`Condition 2`[which(Relationship == "is_selection_for")]
	})

	bigTableFile <- paste0(dataDir,"counts/allCounts.csv")
	allCounts <- read.csv(bitTableFile)

	#TODO: 
	# * Filter down to nonselect and their WT control
	# * Average over replicates
	# * subtract WT control from nonselect and floor
	# * Iterate over (possibly multiple) nonselects and analyze:
	#   * Coverage map and census for each tile
	#   * Extrapolated census
	#   * Complexity analysis
	#   * Nucleotide bias






	options(op)
	return(NULL)
}
