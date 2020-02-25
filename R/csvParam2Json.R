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

#' run validation checks on parameters object
#'
#' @param params parameter object to validate
#' @return TRUE if everything checks out, otherwise it throws errors
#' @export
validateParameters <- function(params) {

	#Validate template sequence
	if (is.na(params$template$cds_start) || is.na(params$template$cds_end)) {
		stop("CDS start and end must be integer numbers!")
	}
	if (params$template$cds_start < 1 || params$template$cds_start > nchar(params$template$seq)) {
		stop("CDS start is out of bounds!")
	}
	if (params$template$cds_end < 1 || params$template$cds_end > nchar(params$template$seq)) {
		stop("CDS end is out of bounds!")
	}
	if (!grepl("^[ACGT]{3,}$",params$template$seq)) {
		stop("Template sequence is not a valid DNA sequence!")
	}
	cdsLength <- params$template$cds_end - params$template$cds_start + 1
	if (cdsLength %% 3 != 0) {
		stop("CDS end is not in-frame with start position!")
	}
	startCodon <- substr(params$template$seq,params$template$cds_start,params$template$cds_start+2)
	if (startCodon != "ATG") {
		stop("CDS start position is not a start codon (ATG)")
	}

	#Validate Uniprot Accession via Regex
	uniprotRX <- "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"
	if (!grepl(uniprotRX,params$template$uniprot)) {
		stop("Invalid Uniprot Accession!")
	}

	#Check assay parameters
	if (!params$assay[["selection"]] %in% c("Positive","Negative")) {
		stop("Assay Selection must be either 'Positive' or 'Negative'!")
	}

	#check that at least one condition is defined
	if (length(params$conditions$names) < 1) {
		stop("Must define at least one condition!")
	}

	#check condition table for validity
	#for QC runs, no relationships are declared, but in any other case, the following rules must apply:
	if (nrow(params$conditions$definitions) > 0) {
		#check that definitions use valid condition names
		if (!all(params$conditions$definitions[,"Condition 1"] %in% params$conditions$names) || 
			!all(params$conditions$definitions[,"Condition 2"] %in% params$conditions$names)) {
			stop("Undeclared condition name in definitions!")
		}
		#and that they use valid relationship names
		relationships <- c("is_selection_for","is_wt_control_for")
		if (!all(params$conditions$definitions[,"Relationship"] %in% relationships)) {
			stop("Invalid relationship in condition definitions!")
		}
		#check that each conditon has at least one relationship
		used <- sapply(params$conditions$names, function(cname) {
			cname %in% params$conditions$definitions[,"Condition 1"] || 
			cname %in% params$conditions$definitions[,"Condition 2"]
		})
		if (!all(used)) {
			stop("Conditions not in relationships: ",paste(params$conditions$names[!used],collapse=", "))
		}
		#at least one selection condition must be declared:
		if (!("is_selection_for" %in% params$conditions$definitions[,"Relationship"])) {
			stop("No select-nonselect relationship was defined!")
		}
		#Check that each sel/nonsel condition has a WT control
		mainConds <- c(getSelects(params),getNonselects(params))
		hasWT <- sapply(mainConds, function(cond) {
			length(getWTControlFor(cond,params)) > 0
		})
		if (!all(hasWT)) {
			stop("No WT control defined for: ",paste(mainConds[!hasWT],collapse=", "))
		}
	} else {
		warning("No condition definitions detected! Is this a QC run?")
	}
	

	if (is.na(params$numReplicates)) {
		stop("Number of replicates must be an integer number!")
	}
	if (params$numReplicates < 2) {
		stop("Number of replicates must be at least 2!")
	}

	#validate the sample sheet
	if (!all(params$samples[,"Tile ID"] %in% params$tiles[,"Tile Number"])) {
		stop("Undeclared tiles found in sample sheet!")
	}
	if (!all(grepl("^[A-Za-z0-9]+$",params$samples[,"Sample ID"]))) {
		stop("Sample IDs must not contain special characters!")
	}
	if (!all(params$samples[,"Condition"] %in% params$conditions$names)) {
		stop("Undeclared conditions found in sample sheet!")
	}
	if (!all(params$samples[,"Time point"] %in% params$timepoints[,"Time point name"])) {
		stop("Undeclared time points found in sample sheet!")
	}
	if (!all(params$samples[,"Replicate"] <= params$numReplicates)) {
		stop("Undeclared time points found in sample sheet!")
	}

	#TODO: Validate that the entire CDS in covered with tiles and regions

	return(TRUE)
}


#' convert CSV input parameter file to JSON format
#'
#' @param infile the input CSV file
#' @param outfile the output JSON file. Defaults to parameters.json in the same directory
#' @param logger yogilog logger. Defaults to NULL and writes to stdout.
#' @return NULL. Results are written to file.
#' @export
csvParam2Json <- function(infile,outfile=sub("[^/]+$","parameters.json",infile),logger=NULL) {

	op <- options(stringsAsFactors=FALSE)

	#for writing JSON output
	library(RJSONIO)
	#for helper functions
	library(yogitools)

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
	logErr <- function(...) {
		if (!is.null(logger)) {
			logger$error(...)
		} else {
			do.call(cat,c("ERROR:",list(...),"\n"))
		}
	}


	#check that the file is indeed a csv file and can be read
	if(!grepl("\\.csv$",infile)) {
		stop("Input file must be CSV file!")
	}
	if (!canRead(infile)) {
		stop("Input file cannot be read!")
	}

	#read the file into a list of lists and extract the first column
	csv <- strsplit(scan(infile,what="character",sep="\n",quiet=TRUE),",")
	#remove any potential quotes introduced by Excel
	csv <- lapply(csv,function(row) {
		sapply(row,function(cell) {
			gsub("\"","",cell)
		})
	})
	#Extract the first column
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
	output$template$geneName <- csv[[getRow("Gene name:")]][[2]]
	output$template$seq <- csv[[getRow("Sequence:")]][[2]]
	output$template$cds_start <- as.integer(csv[[getRow("CDS start:")]][[2]])
	output$template$cds_end <- as.integer(csv[[getRow("CDS end:")]][[2]])
	output$template$uniprot <- csv[[getRow("Uniprot Accession:")]][[2]]

	#extract assay parameters
	output$assay <- list()
	output$assay$type <- csv[[getRow("Assay Type:")]][[2]]
	output$assay$selection <- csv[[getRow("Selection:")]][[2]]

	#extract conditions, replicates and timepoints
	conditions <- csv[[getRow("List of conditions:")]][-1]
	conditions <- conditions[!(conditions=="" | is.na(conditions))]
	output$conditions <- list()
	output$conditions$names <- as.vector(conditions)
	output$numReplicates <-  as.integer(csv[[getRow("Number of Replicates:")]][[2]])
	output$numTimepoints <-  as.integer(csv[[getRow("Number of time points:")]][[2]])

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
	output$conditions$definitions <- conditionTable

	#Time point definitions
	timeTable <- as.data.frame(extractTable(firstField="Time point name",nextSection="Sequencing samples"))
	timeTable$Time <- as.numeric(timeTable$Time)
	output$timepoints <- timeTable


	#Extract sample sheet
	sampleTable <- as.data.frame(extractTable(firstField="Sample ID",nextSection=""))
	sampleTable$`Tile ID` <- as.integer(sampleTable$`Tile ID`)
	sampleTable$Replicate <- as.integer(sampleTable$Replicate)
	output$samples <- sampleTable

	#Run validation on all parameters
	validateParameters(output)

	#convert output to JSON and write to file
	logInfo("Writing output to",outfile,"\n")
	
	con <- file(outfile,open="w")
	writeLines(toJSON(output),con)
	close(con)

	options(op)

	logInfo("Conversion successful!\n")
	invisible(return(NULL))
}

#' parse JSON parameter file
#'
#' @param filename the input CSV file
#' @return the parameter object, as a list of lists
#' @export
parseParameters <- function(filename) { 

	op <- options(stringsAsFactors=FALSE)

	#for writing JSON output
	library(RJSONIO)
	#for helper functions
	library(yogitools)

	#check that the file is indeed a json file and can be read
	stopifnot(grepl("\\.json$",filename), canRead(filename))

	#parse JSON to list of lists
	params <- fromJSON(filename)

	#rebuild tables and dataframes from lists
	params$conditions$definitions <- do.call(rbind,params$conditions$definitions)
	if (!inherits(params$regions,"list")) {
		params$regions <- list(params$regions)
	}
	params$regions <- do.call(rbind,params$regions)
	if (!inherits(params$tiles,"list")) {
		params$tiles <- list(params$tiles)
	}
	params$tiles <- do.call(rbind,params$tiles)

	timeCols <- names(params$timepoints)
	params$timepoints <- as.data.frame(params$timepoints)
	colnames(params$timepoints) <- timeCols

	sampleCols <- names(params$samples)
	params$samples <- as.data.frame(params$samples)
	colnames(params$samples) <- sampleCols

	#run full validation of the parameter object
	validateParameters(params)

	#calculate CDS and protein sequence
	cdsLength <- params$template$cds_end - params$template$cds_start + 1
	cdsSeq <- substr(params$template$seq,params$template$cds_start,params$template$cds_end)

	data(trtable)
	proteinSeq <- sapply(seq(1,cdsLength,3),function(pos) trtable[[substr(cdsSeq,pos,pos+2)]])
	proteinLength <- length(proteinSeq)
	if (proteinSeq[[proteinLength]] == "*") {
		proteinSeq <- proteinSeq[-proteinLength]	
		proteinLength <- proteinLength-1
	}

	params$template$cdsSeq <- cdsSeq
	params$template$cdsLength <- cdsLength
	params$template$proteinSeq <- paste(proteinSeq,collapse="")
	params$template$proteinLength <- proteinLength


	#convert region and tile positions to nucleotide positions
	params$regions <- cbind(params$regions, 
		`Start NC in CDS` = params$regions[,"Start AA"]*3-2,
		`End NC in CDS` = params$regions[,"End AA"]*3
	)
	params$regions <- cbind(params$regions, 
		`Start NC in Template` = params$regions[,"Start NC in CDS"]+ params$template$cds_start - 1,
		`End NC in Template` = params$regions[,"End NC in CDS"]+ params$template$cds_start - 1
	)

	params$tiles <- cbind(params$tiles, 
		`Start NC in CDS` = params$tiles[,"Start AA"]*3-2,
		`End NC in CDS` = params$tiles[,"End AA"]*3
	)
	params$tiles <- cbind(params$tiles, 
		`Start NC in Template` = params$tiles[,"Start NC in CDS"]+ params$template$cds_start - 1,
		`End NC in Template` = params$tiles[,"End NC in CDS"]+ params$template$cds_start - 1
	)


	#Process the condition definitions to annotate the sample table
	sampleTable <- params$samples
	#Assign correct nonselect sample for each selection sample
	sampleTable$nonselectRef <- sapply(1:nrow(sampleTable), function(i)with(sampleTable[i,], {
		nsCond <- with(as.data.frame(params$conditions$definitions),{
			`Condition 2`[which(Relationship == "is_selection_for" & `Condition 1` == Condition)]
		})
		if (length(nsCond) == 0) return(NA)
		rows <- which(
			sampleTable$`Tile ID` == `Tile ID` & 
			sampleTable$Replicate == Replicate & 
			sampleTable$`Time point` == `Time point` & 
			sampleTable$Condition == nsCond
		)
		if (length(rows) == 0) {
			return(NA)
		} else if (length(rows) > 1) {
			logWarn("More than one nonselect match found for sample",`Sample ID`)
			paste(sampleTable$`Sample ID`[rows],collapse=",")
		} else {
			sampleTable$`Sample ID`[rows]
		}
	}))
	#Assign correct WT control sample for each main sample
	sampleTable$wtCtrlRef <- sapply(1:nrow(sampleTable), function(i)with(sampleTable[i,], {
		wtCond <- with(as.data.frame(params$conditions$definitions),{
			`Condition 1`[which(Relationship == "is_wt_control_for" & `Condition 2` == Condition)]
		})
		if (length(wtCond) == 0) return(NA)
		rows <- which(
			sampleTable$`Tile ID` == `Tile ID` & 
			sampleTable$Replicate == Replicate & 
			sampleTable$`Time point` == `Time point` & 
			sampleTable$Condition == wtCond
		)
		if (length(rows) == 0) {
			return(NA)
		} else if (length(rows) > 1) {
			logWarn("More than one WT control match found for sample",`Sample ID`)
			paste(sampleTable$`Sample ID`[rows],collapse=",")
		} else {
			sampleTable$`Sample ID`[rows]
		}
	}))
	params$samples <- sampleTable

	options(op)

	return(params)
}

#' Convenience function to calculate the list of SNV-reachable AA changes
#'
#' Given a parameter object containing a coding sequence, this function will calculate a table
#' of all SNV-reachable amino acid changes, detailing which codon changes correspond to each.
#'
#' @param param the parameter object
#' @return a data.frame containing wt and mutant codons as well as wt and mutant AAs for all possible SNVs
#' @export
reachableChanges <- function(params) {
	data(trtable)
	library(hgvsParseR)
	hgvsp <- new.hgvs.builder.p(aacode=3)
	codons <- sapply(
		seq(1,params$template$cdsLength,3),
		function(s) substr(params$template$cdsSeq,s,s+2)
	)
	changes <- expand.grid(i=1:3,base=toChars("ACGT"),stringsAsFactors=FALSE)
	reachable <- do.call(rbind,lapply(1:length(codons), function(pos) {
		wtcodon <- codons[[pos]]
		wtaa <- trtable[[wtcodon]]
		muts <- as.df(lapply(1:nrow(changes),function(k) with(changes[k,],{
			mutcodon <- wtcodon
			substr(mutcodon,changes[k,"i"],changes[k,"i"]) <- changes[k,"base"]
			mutaa <- trtable[[mutcodon]]
			return(list(
				wtcodon=wtcodon,pos=pos,mutcodon=mutcodon,
				wtaa=wtaa,mutaa=mutaa
			))
		})))
		muts <- muts[with(muts,wtaa!=mutaa),]
		as.df(with(muts,tapply(1:nrow(muts),mutaa,function(is) {
			list(
				wtcodon=unique(wtcodon[is]),pos=unique(pos[is]),
				mutcodons=paste0(mutcodon[is],collapse="|"),
				wtaa=unique(wtaa[is]),mutaa=unique(mutaa[is])
			)
		})))
	}))
	reachable$aachange <- with(reachable,paste0(wtaa,pos,mutaa))
	reachable$hgvsp <- sapply(1:nrow(reachable),function(i) with(reachable[i,],hgvsp$substitution(pos,wtaa,mutaa)))

	return(reachable)
}

#' Convenience function get all selective conditions
#'
#' @param param the parameter object
#' @return the list of selective conditions
#' @export
getSelects <- function(params) {
	unique(with(as.data.frame(params$conditions$definitions),{
		`Condition 1`[which(Relationship == "is_selection_for")]
	}))
}

#' Convenience function get all non-selective conditions
#'
#' @param param the parameter object
#' @return the list of non-selective conditions
#' @export
getNonselects <- function(params) {
	unique(with(as.data.frame(params$conditions$definitions),{
		`Condition 2`[which(Relationship == "is_selection_for")]
	}))
}

#' Convenience function get the matching non-selective condition for the given input condition
#'
#' @param param the parameter object
#' @return the list of matching non-selective conditions
#' @export
getNonselectFor <- function(cond,params) {
	unique(with(as.data.frame(params$conditions$definitions),{
		`Condition 2`[which(Relationship == "is_selection_for" & `Condition 1` == cond)]
	}))
}

#' Convenience function get the matching WT control condition for the given input condition
#'
#' @param param the parameter object
#' @return the list of matching WT control conditions
#' @export
getWTControlFor <- function(cond,params) {
	unique(with(as.data.frame(params$conditions$definitions),{
		`Condition 1`[which(Relationship == "is_wt_control_for" & `Condition 2` == cond)]
	}))
}



