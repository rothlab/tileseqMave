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
buildJointTable <- function(dataDir,paramFile=paste0(dataDir,"parameters.json"),logger=NULL,mc.cores=6) {

	op <- options(stringsAsFactors=FALSE)

	library(yogitools)
	library(hgvsParseR)
	library(parallel)

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

	#Read parameters
	logInfo("Reading parameters")
	params <- parseParameters(paramFile)

	######################
	# BUILD SAMPLE TABLE #
	######################

	#prepare a table of all samples
	logInfo("Accounting for all input count files")
	sampleTable <- params$samples

	#find frequency table files and assign each to its respective sample
	freqfiles <- list.files(paste0(dataDir,"counts"),pattern="freqs_sample\\d+\\.csv",full.names=TRUE)
	filesamples <- as.integer(extract.groups(freqfiles,"freqs_sample(\\d+)\\.csv")[,1])
	sampleTable$countfile <- sapply(sampleTable$`Sample ID`, function(sid) {
		i <- which(filesamples == sid)
		if (is.null(i)) {
			NA
		} else {
			freqfiles[[i]]
		}
	})

	#if any samples/files are unaccounted for, throw an error!
	if (any(is.na(sampleTable$countfile))) {
		missingFiles <- with(sampleTable,`Sample ID`[which(is.na(countfile))])
		stop("Frequency files for samples",paste(missingFiles,collapse=", "),"are missing!")
	}


	#########################
	# Read all input tables #
	#########################

	logInfo("Reading count data")
	allCounts <- lapply(sampleTable$countfile,parseCountFile)


	################
	# Merge tables #
	################

	logInfo("Merging tables")
	#find the union of all variants
	allVars <- Reduce(union,lapply(allCounts,function(x)x$HGVS))
	#and build index of the matching translations
	trIdx <- hash()
	for (i in 1:length(allCounts)) {
		values(trIdx,keys=allCounts[[i]][,1]) <- allCounts[[i]][,2]
	}

	#make lists of unique conditions, timepoints and replicates
	conditions <- unique(sampleTable$Condition)
	timepoints <- unique(sampleTable$`Time point`)
	replicates <- unique(sampleTable$Replicate)

	#helper function to combine a set of count/frequency tables, ordered by allVars
	joinTables <- function(tables, allVars) {
		#extract the union of variants
		cnts <- hash(allVars,0)
		freqs <- hash(allVars,0)
		for (i in 1:length(tables)) {
			for (j in 1:nrow(tables[[i]])) {
				mut <- tables[[i]][j,"HGVS"]
				#TODO: think about whether summing is appropriate or if max() is better
				cnts[[mut]] <- cnts[[mut]] + tables[[i]][j,"count"]
				freqs[[mut]] <- freqs[[mut]] + tables[[i]][j,"frequency"]
			}
		}
		data.frame(
			count=values(cnts,keys=allVars),
			frequency=values(freqs,keys=allVars)
		)
	}

	#join all samples together into one comprehensive table
	condTables <- mclapply(conditions, function(cond) {
		tpTables <- lapply(timepoints, function(tp) {
			replTables <- lapply(replicates, function(repl) {
				#find the set of samples that represent the tiles for this specific
				# combination of condition / timepoint / replicate
				is <- with(sampleTable,{
					which(Condition == cond & `Time point` == tp & Replicate == repl)
				})
				#And join into a single table
				joinTables(allCounts[is],allVars)
			})
			#bind columns together
			names(replTables) <- paste0("rep",replicates)
			do.call(cbind,replTables)
		})
		#bind columns together
		names(tpTables) <- paste0("t",timepoints)
		do.call(cbind,tpTables)
	},mc.cores=min(length(conditions),mc.cores))
	names(condTables) <- conditions
	#add back HGVS strings and translations in the final column binding step
	jointTable <- cbind(
		HGVS=allVars,
		HGVS_pro=values(trIdx,keys=allVars),
		do.call(cbind,condTables)
	)

	########################
	# WRITE OUTPUT TO FILE #
	########################
	logInfo("Writing results to file.")
	outfile <- paste0(dataDir,"counts/allCounts.csv")
	write.csv(jointTable,outfile,row.names=FALSE)
	
	logInfo("Done.")

	options(op)
	return(NULL)
}
