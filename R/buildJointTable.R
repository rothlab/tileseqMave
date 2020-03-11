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
buildJointTable <- function(dataDir,paramFile=paste0(dataDir,"parameters.json"),logger=NULL,mc.cores=6,srOverride=FALSE) {

	op <- options(stringsAsFactors=FALSE)

	library(yogitools)
	library(hash)
	library(hgvsParseR)
	library(parallel)
	library(pbmcapply)

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
	params <- parseParameters(paramFile,srOverride=srOverride)

	######################
	# BUILD SAMPLE TABLE #
	######################

	#prepare a table of all samples
	logInfo("Accounting for all input count files")
	sampleTable <- params$samples

	#find counts folder
	subDirs <- list.dirs(dataDir,recursive=FALSE)
	countDirs <- subDirs[grepl("_mut_call$",subDirs)]
	if (length(countDirs) == 0) {
		stop("No mutation call output found!")
	}
	latestCountDir <- sort(countDirs,decreasing=TRUE)[[1]]
	#extract time stamp
	timeStamp <- extract.groups(latestCountDir,"/(\\d{4}-\\d{2}-\\d{2}-\\d{2}-\\d{2}-\\d{2})")[1,1]

	#find count table files and assign each to its respective sample
	countfiles <- list.files(latestCountDir,pattern="counts_sample\\d+\\.csv$",full.names=TRUE)
	filesamples <- as.integer(extract.groups(countfiles,"counts_sample(\\d+)\\.csv")[,1])
	sampleTable$countfile <- sapply(sampleTable$`Sample ID`, function(sid) {
		i <- which(filesamples == sid)
		if (length(i) == 0 || is.null(i)) {
			NA
		} else {
			countfiles[[i]]
		}
	})

	#if any samples/files are unaccounted for, throw an error!
	if (any(is.na(sampleTable$countfile))) {
		missingFiles <- with(sampleTable,`Sample ID`[which(is.na(countfile))])
		stop("Missing count files for the following samples: ",paste(missingFiles,collapse=", "))
	}


	#####################################
	# Read and prepare all input tables #
	#####################################

	logInfo("Reading count data")
	allCounts <- lapply(sampleTable$countfile,parseCountFile)

	#extract sequencing depths for each sample
	sampleTable$depth <- sapply(allCounts,function(counts) as.integer(attr(counts,"depth"))) 
	#and save the sample table for future reference
	logInfo("Exporting sequencing depth information.")
	outfile <- paste0(latestCountDir,"/sampleDepths.csv")
	write.csv(sampleTable,outfile,row.names=FALSE)
	

	#Calculate frequencies for all counts
	allCounts <- lapply(allCounts,function(counts) {
		counts$frequency <- counts$count/as.integer(attr(counts,"depth"))
		counts
	})


	#find the union of all variants
	allVars <- Reduce(union,lapply(allCounts,function(x)x$HGVS))
	logInfo(sprintf(
		"Translating %d unique variant calls to protein level. This may take some time...", 
		length(allVars)
	))
	#translate them to amino acid level
	builder <- new.hgvs.builder.p(aacode=3)
	cbuilder <- new.hgvs.builder.c()
	# transTable <- as.df(pbmclapply(allVars, translateHGVS, params, builder, cbuilder, mc.cores=mc.cores))
	transTable <- as.df(pbmclapply(allVars, function(mut) {
		tryCatch({
			translateHGVS(mut, params, builder, cbuilder)
		},error=function(e){
			c(
				hgvsp=as.character(e),codonChanges=NA,codonHGVS=NA,
				aaChanges=NA,aaChangeHGVS=NA
			)
		})
	},mc.cores=mc.cores))

	if (any(is.na(transTable[,2]))) {
		for (i in which(is.na(transTable[,2]))) {
			logWarn("Translation for variant ",allVars[[i]],"failed: ",transTable[i,1])
		}
		stop("Translations for ",sum(is.na(transTable[,2]))," variants failed! See log for details.")
	}


	#TODO: detect failed translations


	#and build index of the translations
	# trIdx <- hash(allVars,1:nrow(transTable))
	# for (i in 1:length(allCounts)) {
	# 	values(trIdx,keys=allCounts[[i]][,1]) <- allCounts[[i]][,2]
	# }

	################
	# Merge tables #
	################

	logInfo("Merging tables...")

	#make lists of unique conditions, timepoints and replicates
	conditions <- params$conditions$names
	# conditions <- unique(sampleTable$Condition)
	# timepoints <- unique(sampleTable$`Time point`)
	# replicates <- unique(sampleTable$Replicate)

	#helper function to combine a set of count/frequency tables, ordered by allVars
	joinTables <- function(tables, allVars) {
		out <- data.frame(count=rep(0,length(allVars)),frequency=rep(0,length(allVars)),row.names=allVars)
		for (i in 1:length(tables)) {
			out[tables[[i]]$HGVS,] <- out[tables[[i]]$HGVS,] + tables[[i]][,c("count","frequency")]
		}
		return(out)
	}

	#join all samples together into one comprehensive table
	condTables <- lapply(conditions, function(cond) {
		#get the timepoints that are valid for the current condition
		timepoints <- params$timepoints$`Time point name`[1:params$numTimepoints[[cond]]]
		tpTables <- lapply(timepoints, function(tp) {
			logInfo("Processing condition",cond,"timepoint",tp)
			#get the replicates for this condition
			replicates <- 1:params$numReplicates[[cond]]
			replTables <- lapply(replicates, function(repl) {
				#find the set of samples that represent the tiles for this specific
				# combination of condition / timepoint / replicate
				is <- with(sampleTable,{
					which(Condition == cond & `Time point` == tp & Replicate == repl)
				})
				#And join into a single table
				joinTables(allCounts[is],allVars)
			# },mc.cores=min(length(replicates),mc.cores))
			})
			#bind columns together
			names(replTables) <- paste0("rep",replicates)
			do.call(cbind,replTables)
		})
		#bind columns together
		names(tpTables) <- paste0("t",timepoints)
		do.call(cbind,tpTables)
	})
	names(condTables) <- conditions
	#add back HGVS strings and translations in the final column binding step
	jointTable <- cbind(
		HGVS=allVars,
		# HGVS_pro=values(trIdx,keys=allVars),
		transTable,
		do.call(cbind,condTables)
	)

	########################
	# WRITE OUTPUT TO FILE #
	########################
	logInfo("Writing results to file.")
	outfile <- paste0(latestCountDir,"/allCounts.csv")
	write.csv(jointTable,outfile,row.names=FALSE)
	
	##################################
	# CALCULATE MARGINAL FREQUENCIES #
	##################################

	logInfo("Calculating marginal frequencies")

	codonChangeHGVSs <- strsplit(gsub("c\\.|c\\.\\[|\\]","",transTable$codonHGVS),";")
	codonChangeHGVSs <- lapply(codonChangeHGVSs, function(x) paste0("c.",x))
	aaChangeHGVSs <- strsplit(gsub("p\\.|p\\.\\[|\\]","",transTable$aaChangeHGVS),";")
	aaChangeHGVSs <- lapply(aaChangeHGVSs, function(x) paste0("p.",x))
	codonChangeStrs <- strsplit(transTable$codonChanges,"\\|")
	aaChangeStrs <- strsplit(transTable$aaChanges,"\\|")
	
	#build codon change index
	ccIdx <- hash()
	ccStrIdx <- hash()
	aaIdx <- hash()
	aaStrIdx <- hash()
	for (i in 1:length(codonChangeHGVSs)) {
		js <- min(length(codonChangeHGVSs[[i]]), length(aaChangeHGVSs[[i]]))
		for (j in 1:js) {
			cc <- codonChangeHGVSs[[i]][[j]]
			if (has.key(cc,ccIdx)) {
				ccIdx[[cc]] <- c(ccIdx[[cc]],i)
			} else {
				ccIdx[[cc]] <- i
				aaIdx[[cc]] <- aaChangeHGVSs[[i]][[j]]
				ccStrIdx[[cc]] <- codonChangeStrs[[i]][[j]]
				aaStrIdx[[cc]] <- aaChangeStrs[[i]][[j]]
			}
		}
	}

	#build marginals using index
	marginalCCs <- keys(ccIdx)
	marginalCounts <- as.df(pbmclapply(marginalCCs, function(cc) {
		c(
			list(hgvsc=cc,hgvsp=aaIdx[[cc]]),
			codonChange=ccStrIdx[[cc]],aaChange=aaStrIdx[[cc]],
			colSums(jointTable[ccIdx[[cc]],-(1:6)],na.rm=TRUE)
		)
	},mc.cores=mc.cores))

	logInfo("Writing results to file.")
	outfile <- paste0(latestCountDir,"/marginalCounts.csv")
	write.csv(marginalCounts,outfile,row.names=FALSE)


	logInfo("Done.")

	options(op)
	return(NULL)
}
