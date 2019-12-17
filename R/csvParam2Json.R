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

#' convert CSV input parameter file to JSON format
#'
#' @param infile the input CSV file
#' @return nothing, but an output file will be written in the same directory as the input file
#' @export
csvParam2Json <- function(infile,outfile=NULL) {

	#for writing JSON output
	library(RJSONIO)
	#for helper functions
	library(yogitools)

	#check that the file is indeed a csv file and can be read
	stopifnot(grepl("\\.csv$",infile), canRead(infile))

	if (is.null(outfile)) {
		outfile <- sub("csv$","json",infile)
	}

	#read the file into a list of lists and extract the first column
	csv <- strsplit(scan(infile,what="character",sep="\n",quiet=TRUE),",")
	col1 <- sapply(csv,`[[`,1)

	#helper function to locate a named row
	getRow <- function(rowname) {
		i <- which (col1==rowname)
		if (is.null(i))  stop("Missing field: ",rowname)
		return(i)
	}

	extractTable <- function(firstField,nextSection) {
		iHead <- getRow(firstField)
		#table header
		headers <- csv[[iHead]]
		headers <- headers[!(headers == "" | is.na(headers))]
		#find the end of the table
		iEnd <- iHead
		while (iEnd < length(col1) && !(col1[[iEnd+1]] %in% c("",nextSection))) iEnd <- iEnd+1
		#extract the table data and apply formatting
		rawTable <- do.call(rbind,csv[(iHead+1):iEnd])
		# formTable <- apply(rawTable[,1:length(headers)],c(1,2),as.integer)
		formTable <- rawTable[,1:length(headers),drop=FALSE]
		colnames(formTable) <- headers
		return(formTable)
	}

	#prepare output data structure
	output <- list()

	#extract project name
	output$project <- csv[[getRow("Project name:")]][[2]]

	#extract template sequence
	output$template <- list()
	output$template$seq <- csv[[getRow("Sequence:")]][[2]]
	output$template$cds_start <- as.integer(csv[[getRow("CDS start:")]][[2]])
	output$template$cds_end <- as.integer(csv[[getRow("CDS end:")]][[2]])

	if (is.na(output$template$cds_start) || is.na(output$template$cds_end)) {
		stop("CDS start and end must be integer numbers!")
	}

	#extract conditions, replicates and timepoints
	conditions <- csv[[getRow("List of conditions:")]][-1]
	conditions <- conditions[!(conditions=="" | is.na(conditions))]
	output$conditions <- list()
	output$conditions$names <- conditions
	output$numReplicates <-  csv[[getRow("Number of Replicates:")]][[2]]
	output$numTimepoints <-  csv[[getRow("Number of time points:")]][[2]]

	if (is.na(output$template$cds_start) || is.na(output$template$cds_end)) {
		stop("CDS start and end must be integer numbers!")
	}

	#extract regions table
	regionTable <- extractTable(firstField="Region Number",nextSection="Sequencing Tiles")
	regionTable <- apply(regionTable,2,as.integer)
	output$regions <- regionTable

	#extract tile table
	tileTable <- extractTable(firstField="Tile Number",nextSection="Condition definitions")
	tileTable <- apply(tileTable,2,as.integer)
	output$tiles <- tileTable

	#extract condition definitions
	conditionTable <- extractTable(firstField="Condition 1",nextSection="Time point definitions")
	#check for validity
	if (!all(conditionTable[,"Condition 1"] %in% conditions) || 
		!all(conditionTable[,"Condition 2"] %in% conditions)) {
		stop("Undeclared condition name in definitions!")
	}
	relationships <- c("is_selection_for","is_wt_control_for")
	if (!all(conditionTable[,"Relationship"] %in% relationships)) {
		stop("Invalid relationship in condition definitions!")
	}
	output$conditions$definitions <- conditionTable

	#Time point definitions
	timeTable <- as.data.frame(extractTable(firstField="Time point name",nextSection="Sequencing samples"))
	timeTable$Time <- as.numeric(timeTable$Time)
	output$timepoints <- timeTable


	#Extract sample sheet
	sampleTable <- as.data.frame(extractTable(firstField="Sample ID",nextSection=""))
	sampleTable$`Tile ID` <- as.integer(sampleTable$`Tile ID`)
	sampleTable$Replicate <- as.integer(sampleTable$Replicate)

	#validate the sample sheet
	if (!all(sampleTable[,"Tile ID"] %in% tileTable[,"Tile Number"])) {
		stop("Undeclared tiles found in sample sheet!")
	}
	if (!all(sampleTable[,"Condition"] %in% conditions)) {
		stop("Undeclared conditions found in sample sheet!")
	}
	if (!all(sampleTable[,"Time point"] %in% timeTable[,"Time point name"])) {
		stop("Undeclared time points found in sample sheet!")
	}
	if (!all(sampleTable[,"Replicate"] <= output$numReplicates)) {
		stop("Undeclared time points found in sample sheet!")
	}
	output$samples <- sampleTable



	#convert output to JSON and write to file
	con <- file(outfile,open="w")
	writeLines(toJSON(output),con)
	close(con)

	cat("Conversion successful!\n")
}