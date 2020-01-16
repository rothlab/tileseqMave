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

#' parse a count table
#' 
parseCountFile <- function(filename) {
	op <- options(stringsAsFactors=FALSE)

	lines <- scan(filename,what="character",sep="\n",quiet=TRUE)
	header <- do.call(c,lapply(strsplit(lines[grep("^#",lines)],":\\s+"),function(xs)setNames(xs[[2]],xs[[1]])))
	metadata <- list(
		sample=header[["#Sample"]],
		tile=as.integer(header[["#Tile"]]),
		condition=header[["#Condition"]],
		replicate=as.integer(header[["#Replicate"]]),
		timepoint=as.integer(header[["#Timepoint"]]),
		depth=as.integer(header[["#Read-depth"]])
	)

	countTable <- read.csv(textConnection(lines),comment.char="#")

	attributes(countTable) <- c(attributes(countTable),metadata)

	options(op)

	return(countTable)
}

# translateCountFile <- function(filename,params) {
# 	counts <- parseCountFile(filename)
# 	builder <- new.hgvs.builder.p(aacode=3)
# 	hgvsp <- sapply(foo$HGVS,translateHGVS,params,builder)
# 	counts$HGVS_pro <- hgvsp
# 	counts$frequency <- counts$count/attr(counts,"depth")
# 	return(counts[,c("HGVS","HGVS_pro","count","frequency")])
# }

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
	# #identify corresponding wt controls
	# nsControls <- with(as.data.frame(params$conditions$definitions),{
	# 	sapply(nsConditions,function(nsCondition) {
	# 		`Condition 1`[which(Relationship == "is_wt_control_for" & `Condition 2` == nsCondition)]
	# 	})
	# })
	

	#prepare a table of all samples with and corresponding metadata
	logInfo("Accounting for all input count files")
	sampleTable <- params$samples
	countfiles <- list.files(paste0(dataDir,"counts"),pattern="counts_sample\\d+\\.csv",full.names=TRUE)
	count.samples <- as.integer(extract.groups(countfiles,"counts_sample(\\d+)\\.csv")[,1])

	#confirm that countfiles exist and assign them to each sample
	sampleTable$countfile <- sapply(sampleTable$`Sample ID`, function(sid) {
		i <- which(count.samples == sid)
		if (is.null(i)) {
			NA
		} else {
			countfiles[[i]]
		}
	})

	#filter down to nonselect and their wt controls
	relevantSamples <- with(sampleTable,c(
		`Sample ID`[Condition %in% nsConditions],
		wtCtrlRef[Condition %in% nsConditions]
	))
	relevantSamples <- sampleTable[sampleTable$`Sample ID` %in% relevantSamples,]

	logInfo("Reading count data")
	allCounts <- lapply(relevantSamples$countfile,function(cfile) {
		counts <- parseCountFile(cfile)
		depth <- attr(counts,"depth")
		return(cbind(counts,frequency=counts$count/depth))
	}) 

	joinTables <- function(tables,fun) {
		#extract the union of variants
		freqs <- hash()
		for (i in 1:length(tables)) {
			for (j in 1:nrow(tables[[i]])) {
				mut <- tables[[i]][j,"HGVS"]
				f <- tables[[i]][j,"frequency"]
				if (has.key(mut,freqs)) {
					freqs[[mut]] <- c(freqs[[mut]],f)
				} else {
					freqs[[mut]] <- f
				}
			}
		}
		data.frame(HGVS=keys(freqs),frequency=sapply(keys(freqs),function(k)fun(freqs[[k]])))
	}

	#collapse count tables w.r.t tiles and replicates
	conditions <- unique(relevantSamples$Condition)
	tablesByCondition <- lapply(conditions, function(condID) {
		#join tiles together
		logInfo("joining tiles for",condID)
		freqTables <- lapply(unique(relevantSamples$Replicate), function(replID) {
			is <- with(relevantSamples,{
				which(Condition == condID & Replicate == replID)
			})
			#TODO: think about whether summing is appropriate or if max() is better
			joinTables(allCounts[is],sum)
		})
		logInfo("joining replicates for",condID)
		#average over replicates
		freqTable <- joinTables(freqTables, mean)
	})
	names(tablesByCondition) <- conditions


	#subtract wt controls from nonselect
	normalizedFrequencies <- lapply(nsConditions,function(condID) {
		logInfo("correcting for wt controls in",condID)
		#find the corresponding WT control condition
		wtCondID <- with(as.data.frame(params$conditions$definitions),{
			`Condition 1`[which(Relationship == "is_wt_control_for" & `Condition 2` == condID)]
		})
		#load wtControl into hash
		wtFreqs <- with(tablesByCondition[[wtCondID]],hash(HGVS,frequency))
		currTable <- tablesByCondition[[condID]]
		#subtract from nonselect and floor
		normFreqs <- sapply(1:nrow(currTable), function(i) with(currTable[i,],{
			if (has.key(HGVS,wtFreqs)) {
				max(0,frequency - wtFreqs[[HGVS]])
			} else {
				frequency
			}
		}))
		return(cbind(currTable,normalized.frequency=normFreqs))
	})
	names(normalizedFrequencies) <- nsConditions

	
	#iterate over (possible multiple) nonselective conditions to analyze
	for (condID in nsConditions) {
		normFreqs <- normalizedFrequencies[[condID]]
		#interpret variant descriptors
		logInfo("interpreting variant descriptors for",condID)
		# mutTable <- parseHGVS(normFreqs$HGVS)
		#translate
	}	

	options(op)
	return(NULL)
}
