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
#' @param srOverride single replicate override - allows single replicates
#' @return TRUE if everything checks out, otherwise it throws errors
#' @export
validateParameters <- function(params,srOverride=FALSE) {

	if (srOverride) {
		logWarn("SINGLE REPLICATE OVERRIDE HAS BEEN ENABLED. This means:
 * No quality filtering can be performed!
 * The scoring function will be unable to generate error estimates!
 * The selectionQC function will be unable to create correlation plots!
 * Imputation and refinement will not be possible!
 * Downstream scripts and functions may be unable to process the data!
")
	}

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

	if (!all(grepl("^[A-Za-z][A-Za-z0-9]*$",params$conditions$names))) {
		stop("Condition names must be strictly alpha-numeric and start with a letter.")
	}

	#check condition table for validity
	#for QC runs, no relationships are declared, but in any other case, the following rules must apply:
	if (nrow(params$conditions$definitions) > 0) {
		#check that definitions use valid condition names
		defconds <- unique(c(
			params$conditions$definitions[,"Condition 1"],
			params$conditions$definitions[,"Condition 2"]
		))
		if (!all(defconds %in% params$conditions$names)) {
			culprits <- defconds[which(!(defconds %in% params$conditions$names))]
			stop("The following condition(s) listed in the relationship definitions was(were) not declared in the list of conditions: ",
				paste(culprits,collapse=", ")
			)
		}
		#and that they use valid relationship names
		relationships <- c("is_selection_for","is_wt_control_for")
		if (!all(params$conditions$definitions[,"Relationship"] %in% relationships)) {
			rels <- params$conditions$definitions[,"Relationship"]
			culprits <- unique(rels[which(!(rels %in% relationships))])
			stop("Invalid relationship(s) in condition definitions: ",paste(culprits,collapse=", "))
		}
		#check that each condition has at least one relationship
		used <- sapply(params$conditions$names, function(cname) {
			cname %in% params$conditions$definitions[,"Condition 1"] || 
			cname %in% params$conditions$definitions[,"Condition 2"]
		})
		if (!all(used)) {
			stop(
				"All conditions must have at least one defined relationship! Missing: ",
				paste(params$conditions$names[!used],collapse=", ")
			)
		}
		#at least one selection condition must be declared:
		if (!("is_selection_for" %in% params$conditions$definitions[,"Relationship"])) {
			logWarn("No select-nonselect relationship was defined! Is this a QC run?")
		} else {
  		#make sure that no selection condition has more than one nonselect partner
  		numNS <- sapply(getSelects(params), function(sCond) length(getNonselectFor(sCond,params)))
  		if (any(numNS > 1)) {
  			culprits <- getSelects(params)[which(numNS > 1)]
  			stop("Each selection conditon may only have one corresponding nonselect condition! Violations: ",culprits)
  		}
		}
		#Check that each sel/nonsel condition has a WT control
		mainConds <- c(getSelects(params),getNonselects(params))
		hasWT <- sapply(mainConds, function(cond) {
			length(getWTControlFor(cond,params)) > 0
		})
		if (!all(hasWT)) {
			logWarn("No WT control defined for: ",paste(mainConds[!hasWT],collapse=", "))
		}
	} else {
		logWarn("No condition definitions detected! Is this a QC run?")
	}
	

	if (any(is.na(params$numReplicates))) {
		stop("Number of replicates must be integer numbers for each condition!")
	}
	minRep <- if (!srOverride && nrow(params$conditions$definitions) > 0 && length(getSelects(params)) > 0) 2 else 1
	if (any(params$numReplicates < minRep)) {
		stop("Number of replicates must be at least ",minRep,"!")
	}

	if (any(is.na(params$numTimepoints))) {
		stop("Number of time points must be integer numbers for each condition!")
	}
	if (any(params$numTimepoints < 1)) {
		stop("Each condition must have at least 1 timepoint!")
	}
	if (!all(params$numTimepoints == params$numTimepoints[[1]])) {
		stop("Differing numbers of timepoints per condition are not yet supported!")
	}



	#validate the sample sheet
	if (!all(params$samples[,"Tile ID"] %in% params$tiles[,"Tile Number"])) {
		tid <- params$samples[,"Tile ID"]
		culprits <- unique(tid[which(!(tid %in% params$tiles[,"Tile Number"]))])
		stop("Undeclared tiles found in sample sheet:",
			paste(culprits,collapse=", "),
			"\nPlease declare them in the Tile table!"
		)
	}
	if (!all(grepl("^[A-Za-z0-9-]+$",params$samples[,"Sample ID"]))) {
		sid <- params$samples[,"Sample ID"]
		culprits <- sid[which(!grepl("^[A-Za-z0-9-]+$",sid))]
		stop("Sample IDs must not contain special characters except minus signs!\n",
			"Violations: ",paste(culprits,collapse=", ")
		)
	}
	if (!all(params$samples[,"Condition"] %in% params$conditions$names)) {
		cid <- params$samples[,"Condition"]
		culprits <- unique(cid[which(!(cid %in% params$conditions$names))])
		stop("Undeclared conditions found in sample sheet:",
			paste(culprits,collapse=", "),
			"\nPlease declare them in the list of conditions!"
		)
	}
	if (!all(params$samples[,"Time point"] %in% params$timepoints[,"Time point name"])) {
		tid <- params$samples[,"Time point"]
		culprits <- unique(tid[which(!(tid %in% params$timepoints[,"Time point name"]))])
		stop("Undeclared time points found in sample sheet:",
			paste(culprits,collapse=", "),
			"\nPlease declare them in the time point definition section!"
		)
	}
	repsPerCon <- tapply(params$samples[,"Replicate"],params$samples[,"Condition"],table)
	invisible(lapply(params$conditions$names, function(cond) {
		if (!(length(repsPerCon[[cond]]) == params$numReplicates[[cond]])) {
			stop("Number of replicate samples for condition ",cond," does not match declared number of replicates!")
		}
		if (!all(repsPerCon[[cond]] == repsPerCon[[cond]][[1]])) {
			stop("Inconsistent number of samples for condition ",cond,"! (Some samples have more replicates than others!)")
		}
	}))

	#validate metaparameters
	if (is.na(params$varcaller$posteriorThreshold)) {
		stop("Posterior threshold must be numeric!")
	}
	if (params$varcaller$posteriorThreshold < 0.5 || params$varcaller$posteriorThreshold >= 1) {
		stop("Posterior threshold must be within greater or equal to 0.5 and less than 1!")
	}
	if (is.na(params$varcaller$minCover)) {
		stop("Minimum cover parameter must be numeric!")
	}
	if (params$varcaller$minCover < 0 || params$varcaller$minCover > 1) {
		stop("Minimum cover must be between 0 and 1!")
	}
	if (is.na(params$varcaller$mutRate)) {
		stop("Mutation rate must be numeric!")
	}
	if (params$varcaller$mutRate < 0 || params$varcaller$mutRate > 1) {
		stop("Mutation rate must be between 0 and 1!")
	}
	if (params$varcaller$mutRate > 0.1) {
		logWarn("Mutation rate parameter is set unusually high! Are you sure about this?")
	}
	if (is.na(params$scoring$countThreshold)) {
		stop("Minimum read count must be an integer!")
	}
	if (params$scoring$countThreshold < 1) {
		stop("Minimum read count must be at least 1!")
	}
	if (is.na(params$scoring$pseudo.n)) {
		stop("Pseudo-replicates must be an integer!")
	}
	if (params$scoring$pseudo.n < 1) {
		stop("Pseudo-replicates must be at least 1!")
	}
	if (is.na(params$scoring$sdThreshold)) {
		stop("SD threshold must be numeric!")
	}
	if (params$scoring$sdThreshold <= 0) {
		stop("SD threshold must be greater than 0!")
	}
	if (is.na(params$scoring$wtQuantile)) {
	  stop("WT filter quantile must be numeric!")
	}
	if (params$scoring$wtQuantile < 0.5 || params$scoring$wtQuantile >= 1) {
	  stop("WT filter quantile only allows values from interval [0.5;1[ ")
	}
	if (is.na(params$scoring$cvDeviation)) {
	  stop("Replicate disagreement factor (cvDeviation) must be numeric!")
	}
	if (params$scoring$cvDeviation < 1) {
	  stop("Replicate disagreement factor (cvDevation) must be greater than 1!")
	}

	#validate normalization table
	if (length(params$normalization) > 0 && nrow(params$normalization) > 0) {
		if (!all(params$normalization$Condition %in% getSelects(params))) {
			stop("Normalization overrides can only be defined for selection conditions.")
		}
		if (!all(params$normalization$`Time point` %in% params$timepoints$`Time point name`)) {
			stop("Invalid time-point in normalization override table!")
		}
		if (!all(params$normalization$Region %in% params$regions$`Region Number`)) {
			stop("Invalid region(s) in normalization override table!")
		}
		if (!all(params$normalization$Type %in% c("synonymous","nonsense"))) {
			stop("Invalid entries in normalization override table: Type must be 'synonymous' or 'nonsense'!")
		}
		if (any(is.na(params$normalization$Value))) {
			stop("Invalid entries in normalization override table: Value must be numeric!")
		}
	}

	#TODO: Validate that the entire CDS in covered with tiles and regions

	return(TRUE)
}


#' convert CSV input parameter file to JSON format
#'
#' @param infile the input CSV file
#' @param outfile the output JSON file. Defaults to parameters.json in the same directory
#' @return NULL. Results are written to file.
#' @export
csvParam2Json <- function(infile,outfile=sub("[^/]+$","parameters.json",infile),srOverride=FALSE) {

	op <- options(stringsAsFactors=FALSE)

	#for writing JSON output
	library(RJSONIO)
	#for helper functions
	library(yogitools)


	#check that the file is indeed a csv file and can be read
	if(!grepl("\\.csv$",infile)) {
		stop("Input file must be CSV file!")
	}
	if (!canRead(infile)) {
		stop("Input file cannot be read!")
	}

	#read the file into a list of lists
	csv <- strsplit(scan(infile,what="character",sep="\n",quiet=TRUE),",")
	#remove any potential quotes introduced by Excel
	csv <- lapply(csv,function(row) {
		sapply(row,function(cell) {
			trimws(gsub("\"","",cell))
		})
	})
	#Extract the first column
	col1 <- sapply(csv,`[[`,1)

	#Helper function to check if row exists
	hasRow <- function(rowname) {
		any(col1==rowname)
	}

	#helper function to locate a named row
	getRow <- function(rowname) {
		if (!hasRow(rowname)) {
			stop("Missing section: ",rowname)
		}
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
		#if the table is empty...
		if (iEnd == iHead) {
			return(data.frame())
			# return(matrix(ncol=length(headers),nrow=0,dimnames=list(NULL,headers)))
		} else {
			#extract the table data and apply formatting
			rawTable <- do.call(rbind,csv[(iHead+1):iEnd])
			# formTable <- apply(rawTable[,1:length(headers)],c(1,2),as.integer)
			formTable <- rawTable[,1:length(headers),drop=FALSE]
			colnames(formTable) <- headers
			return(formTable)
		}
	}

	#prepare output data structure
	output <- list()

	#extract project name
	output$project <- csv[[getRow("Project name:")]][[2]]

	#extract template sequence
	output$template <- list()
	output$template$geneName <- csv[[getRow("Gene name:")]][[2]]
	output$template$seq <- toupper(csv[[getRow("Sequence:")]][[2]])
	output$template$cds_start <- as.integer(csv[[getRow("CDS start:")]][[2]])
	output$template$cds_end <- as.integer(csv[[getRow("CDS end:")]][[2]])
	output$template$uniprot <- csv[[getRow("Uniprot Accession:")]][[2]]

	#extract assay parameters
	output$assay <- list()
	output$assay$name <- csv[[getRow("Assay Name:")]][[2]]
	# output$assay$description <- csv[[getRow("Assay Description:")]][[2]]
	output$assay$selection <- switch(csv[[getRow("Negative selection?")]][[2]],
		Yes="Negative", No="Positive", stop("Negative selection must be 'Yes' or 'No'!")
	)
	

	#extract conditions, replicates and timepoints
	conditions <- csv[[getRow("Condition IDs:")]][-1]
	conditions <- conditions[!(conditions=="" | is.na(conditions))]
	output$conditions <- list()
	output$conditions$names <- as.vector(conditions)
	reps <- as.integer(csv[[getRow("Number of Replicates:")]][-1])
	if (length(reps) < length(conditions)) {
		stop("Replicate numbers missing for some conditions!")
	}
	reps <- setNames(reps[1:length(conditions)],conditions)
	output$numReplicates <- reps
	tps <- as.integer(csv[[getRow("Number of time points:")]][-1])
	if (length(tps) < length(conditions)) {
		stop("Replicate numbers missing for some conditions!")
	}
	tps <- setNames(tps[1:length(conditions)],conditions)
	output$numTimepoints <- tps

	#extract regions table
	regionTable <- as.data.frame(extractTable(firstField="Region Number",nextSection="Sequencing Tiles"))
	# regionTable <- apply(regionTable,2,as.integer)
	regionTable$`Region Number` <- as.integer(regionTable$`Region Number`)
	regionTable$`Start AA` <- as.integer(regionTable$`Start AA`)
	regionTable$`End AA` <- as.integer(regionTable$`End AA`)
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

	#Extract meta-parameters
	output$varcaller <- list()
	if (hasRow("Posterior threshold:")) {
		output$varcaller$posteriorThreshold <- as.numeric(csv[[getRow("Posterior threshold:")]][[2]])
	} else {
		logWarn("No posterior threshold specified. Defaulting to 0.5")
		output$varcaller$posteriorThreshold <- 0.5
	}
	if (hasRow("Minimum coverage:")) {
		output$varcaller$minCover <- as.numeric(csv[[getRow("Minimum coverage:")]][[2]])
	} else {
		logWarn("No minimum coverage specified. Defaulting to 0.6")
		output$varcaller$minCover <- 0.6
	}
	if (hasRow("Per-base mutation rate:")) {
		output$varcaller$mutRate <- as.numeric(csv[[getRow("Per-base mutation rate:")]][[2]])
	} else {
		logWarn("No mutation rate specified. Defaulting to 0.0025")
		output$varcaller$mutRate <- 0.0025
	}
	output$scoring <- list()
	if (hasRow("Minimum read count:")) {
		output$scoring$countThreshold <- as.integer(csv[[getRow("Minimum read count:")]][[2]])
	} else {
		logWarn("No minimum read count specified. Defaulting to 1")
		output$scoring$countThreshold <- 1L
	}
	if (hasRow("Pseudo-replicates:")) {
		output$scoring$pseudo.n <- as.integer(csv[[getRow("Pseudo-replicates:")]][[2]])
	} else {
		logWarn("No pseudo-replicates specified. Defaulting to 8")
		output$scoring$pseudo.n <- 8L
	}
	if (hasRow("SD threshold:")) {
		output$scoring$sdThreshold <- as.numeric(csv[[getRow("SD threshold:")]][[2]])
	} else {
		logWarn("No SD threshold specified. Defaulting to 0.3")
		output$scoring$sdThreshold <- 0.3
	}
	if (hasRow("WT filter quantile:")) {
	  output$scoring$wtQuantile <- as.numeric(csv[[getRow("WT filter quantile:")]][[2]])
	} else {
	  logWarn("No WT filter quantile specified. Defaulting to 0.95")
	  output$scoring$wtQuantile <- 0.95
	}
	if (hasRow("Replicate disagreement factor:")) {
	  output$scoring$cvDeviation <- as.numeric(csv[[getRow("Replicate disagreement factor:")]][[2]])
	} else {
	  logWarn("No replicate disagreement factor specified. Defaulting to 10.")
	  output$scoring$cvDeviation <- 10
	}
	
	#Extract normalization overrides
	output$normalization <- list()
	if (hasRow("Score normalization overrides") && hasRow("Condition")) {
		overrideTable <- as.data.frame(extractTable(firstField="Condition",nextSection=""))
		overrideTable$`Region` <- as.numeric(overrideTable$`Region`)
		overrideTable$`Value` <- as.numeric(overrideTable$`Value`)
		output$normalization <- overrideTable
	}


	#Run validation on all parameters
	withCallingHandlers(
		validateParameters(output,srOverride=srOverride),
		warning=function(w)logWarn(conditionMessage(w))
	)

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
parseParameters <- function(filename,srOverride=FALSE) { 

	op <- options(stringsAsFactors=FALSE)

	#for writing JSON output
	library(RJSONIO)
	#for helper functions
	library(yogitools)

	#check that the file is indeed a json file and can be read
	stopifnot(grepl("\\.json$",filename), canRead(filename))

	#parse JSON to list of lists
	params <- fromJSON(filename,nullValue=NA)

	#rebuild tables and dataframes from lists
	params$conditions$definitions <- do.call(rbind,params$conditions$definitions)
	if (is.null(params$conditions$definitions)) {
		params$conditions$definitions <- matrix(nrow=0,ncol=3,dimnames=list(NULL,c("Condition 1","Relationship","Condition 2")))
	}
	# if (!inherits(params$regions,"list")) {
	# 	params$regions <- list(params$regions)
	# }
	# params$regions <- do.call(rbind,params$regions)
	# if (!inherits(params$tiles,"list")) {
	# 	params$tiles <- list(params$tiles)
	# }
	regCols <- names(params$regions)
	params$regions <- as.data.frame(as.list(params$regions))
	colnames(params$regions) <- regCols

	params$tiles <- do.call(rbind,params$tiles)

	timeCols <- names(params$timepoints)
	params$timepoints <- as.data.frame(params$timepoints)
	colnames(params$timepoints) <- timeCols

	sampleCols <- names(params$samples)
	params$samples <- as.data.frame(params$samples)
	colnames(params$samples) <- sampleCols

	if (length(params$normalization) > 0 && length(params$normalization[[1]]) > 0) {
		normCols <- names(params$normalization)
		params$normalization <- as.data.frame(params$normalization)
		colnames(params$normalization) <- normCols
	} else {
		params$normalization <- NULL
	}

	#make sure lists remain lists
	params$varcaller <- as.list(params$varcaller)
	params$scoring <- as.list(params$scoring)

	#run full validation of the parameter object
	validateParameters(params,srOverride=srOverride)

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

	#convenience functions to retrieve tile and region information by name
	params$regi <- function(region) {
		params$regions[which(params$regions[,"Region Number"] %in% region),]
	}
	params$tili <- function(tile) {
		params$tiles[which(params$tiles[,"Tile Number"] %in% tile),,drop=FALSE]
	}
	params$pos2tile <- function(pos) {
		rows <- sapply(pos,function(pos) {
			i <- which(params$tiles[,"Start AA"] <= pos & params$tiles[,"End AA"] >= pos)
			if (length(i)==0) NA else i
		}) 
		params$tiles[rows,"Tile Number"]
	}
	params$pos2reg <- function(pos) {
		rows <- sapply(pos,function(pos) {
			i <- which(params$regions[,"Start AA"] <= pos & params$regions[,"End AA"] >= pos)
			if (length(i)==0) NA else i
		})
		params$regions[rows,"Region Number"]
	}

	options(op)

	return(params)
}

#' Convenience function to find the tiles that belong to each region.
#'
#' Given the parameter object returns a list of vectors containing the tile IDs for each region ID
#'
#' @param param the parameter object
#' @return a  list of vectors containing the tile IDs for each region ID
#' @export
tilesInRegions <- function(params) {
	setNames(lapply(params$regions$`Region Number`, function(ri) {
		# rs <- params$regions[ri,"Start AA"]
		rs <- params$regi(ri)[["Start AA"]]
		# re <- params$regions[ri,"End AA"]
		re <- params$regi(ri)[["End AA"]]
		tileRows <- which(sapply(params$tiles[,"Start AA"], function(ts){
			ts >=rs && ts < re
		}))
		params$tiles[tileRows,"Tile Number"]
	}),params$regions$`Region Number`)
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

#' Convenience function to find a subdirectory matching the given pattern.
#' 
#' The subfolders should follow the scheme [<label>_]<timestamp>_<outputType>/
#' 
#' @param parentDir the parent folder in which to look for matching subfolders
#' @param pattern a regex pattern used to identify the correct subfolder
#' @return a vector containing the path to the subfolder, its timestamp and any possible label. 
#' @export
latestSubDir <- function(parentDir,pattern="_mut_call$|mut_count$") {
	subDirs <- list.dirs(parentDir,recursive=FALSE)
	subDirs <- subDirs[grepl(pattern,subDirs)]
	if (length(subDirs) == 0) {
		stop("No applicable input folder found!")
	}
	#extract time stamp
	labelsAndTimes <- extract.groups(subDirs,"([^/]+_)?(\\d{4}-\\d{2}-\\d{2}-\\d{2}-\\d{2}-\\d{2})")
	#select the newest dataset
	latestIdx <- which.max(order(labelsAndTimes[,2]))
	return(c(
		dir=subDirs[[latestIdx]],
		timeStamp=labelsAndTimes[latestIdx,2],
		label=labelsAndTimes[latestIdx,1]
	))
}
