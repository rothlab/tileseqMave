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
#' @param inDir input directory, defaults to subdirectory with latest timestamp ending in _mut_count.
#' @param outDir output directory, defaults to input directory.
#' @param paramFile input parameter file. defaults to <dataDir>/parameters.json
#' @param mc.cores the number of CPU cores to use in parallel
#' @param srOverride the single-replicate override flag
#' @return NULL. Results are written to file.
#' @export
buildJointTable <- function(dataDir,inDir=NA,outDir=NA,paramFile=paste0(dataDir,"parameters.json"),
                            mc.cores=6,srOverride=FALSE) {

	op <- options(stringsAsFactors=FALSE)

	library(yogitools)
	library(hash)
	library(hgvsParseR)
	library(parallel)
	library(pbmcapply)

	#make sure data and out dir exist and ends with a "/"
	if (!grepl("/$",dataDir)) {
		dataDir <- paste0(dataDir,"/")
	}
	if (!dir.exists(dataDir)) {
		#we don't use the logger here, assuming that whichever script wraps our function
		#catches the exception and writes to the logger (or uses the log error handler)
		stop("Workspace folder ",dataDir," does not exist!")
	}
	
	
	#Read parameters
	if (!canRead(paramFile)) {
	  stop("Unable to read parameter file!")
	}
	logInfo("Reading parameters from",normalizePath(paramFile))
	params <- withCallingHandlers(
		parseParameters(paramFile,srOverride=srOverride),
		warning=function(w)logWarn(conditionMessage(w))
	)
	
	#configure input and output folders
	if (is.na(inDir)) {
  	latest <- latestSubDir(parentDir=dataDir,pattern="_mut_call$|mut_count$")
  	inDir <- latest[["dir"]]
	} else {
	  if (!dir.exists(inDir)) {
	    stop("Input folder ",inDir," does not exist!")
	  }
	}
	if (!grepl("/$",inDir)) {
	  inDir <- paste0(inDir,"/")
	}
	if (is.na(outDir)) {
	  outDir <- inDir
	} else if (!dir.exists(outDir)) {
	  dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
	}

	if (!grepl("/$",outDir)) {
	  outDir <- paste0(outDir,"/")
	}
	
	logInfo("Using input directory",inDir,"and output directory",outDir)
	
	
	######################
	# BUILD SAMPLE TABLE #
	######################
	
	#prepare a table of all samples
	logInfo("Accounting for all input count files")
	sampleTable <- params$samples
	
	#find count table files and assign each to its respective sample
	countfiles <- list.files(inDir,pattern="counts_sample_.+\\.csv$",full.names=TRUE)
	if (length(countfiles)==0) {
		error("No count files found in ",inDir,"! (Are they named correctly?)")
	}
	filesamples <- extract.groups(countfiles,"counts_sample_(.+)\\.csv$")[,1]
	sampleTable$countfile <- sapply(as.character(sampleTable$`Sample ID`), function(sid) {
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
	allCounts <- lapply(sampleTable$countfile,function(cfile) {
		tryCatch(
			parseCountFile(cfile),
			error=function(e) stop("Error reading file ",cfile," : ",conditionMessage(e))
		)
	})

	#extract sequencing depths for each sample
	sampleTable$depth <- sapply(allCounts,function(counts) as.integer(attr(counts,"depth"))) 
	#and save the sample table for future reference
	logInfo("Exporting sequencing depth information.")
	outfile <- paste0(outDir,"/sampleDepths.csv")
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

	#Check for errors
	if (any(is.na(transTable[,2]))) {
		for (i in which(is.na(transTable[,2]))) {
			logWarn("Translation for variant ",allVars[[i]],"failed: ",transTable[i,1])
		}
		stop("Translations for ",sum(is.na(transTable[,2]))," variants failed! See log for details.")
	}


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
	  timepoints <- getTimepointsFor(cond,params)
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
	#after joining we don't need the individual tables anymore, so we can remove them
	#to save RAM
	rm(condTables)

	########################
	# WRITE OUTPUT TO FILE #
	########################
	logInfo("Writing results to file.")
	outfile <- paste0(outDir,"/allCounts.csv")
	cat("# COMBINATORY VARIANT COUNTS #",
	    "\n# project name:", params$project,
	    "\n# gene name:",params$template$geneName,
	    "\n# tileseqMave version:",as.character(packageVersion("tileseqMave")),
	    "\n# parameter sheet:",normalizePath(paramFile),"\n",
	    file=outfile
	)
	suppressWarnings(
	  write.table(jointTable,outfile,sep=",",append=TRUE,row.names=FALSE,qmethod="double")
	)
	
	##################################
	# CALCULATE MARGINAL FREQUENCIES #
	##################################

	logInfo("Calculating marginal frequencies")

	#split codon change HGVSs by codon
	codonChangeHGVSs <- strsplit(gsub("c\\.|c\\.\\[|\\]","",transTable$codonHGVS),";")
	#and attach the proper prefix
	codonChangeHGVSs <- lapply(codonChangeHGVSs, function(x) paste0("c.",x))
	#do the same for AA change HGVSs
	aaChangeHGVSs <- strsplit(gsub("p\\.|p\\.\\[|\\]","",transTable$aaChangeHGVS),";")
	aaChangeHGVSs <- lapply(aaChangeHGVSs, function(x) paste0("p.",x))
	#now split the native codon changes and AA changes
	codonChangeStrs <- strsplit(transTable$codonChanges,"\\|")
	aaChangeStrs <- strsplit(transTable$aaChanges,"\\|")
	
	cisMerge <- function(hs) {
	  if (length(hs) < 2) return(hs)
	  paste0("c.[",paste(substr(hs,3,nchar(hs)),collapse=";"),"]")
	}
	
	#build codon change index
	ccIdx <- hash()
	ccStrIdx <- hash()
	aaIdx <- hash()
	aaStrIdx <- hash()
	for (i in 1:length(codonChangeHGVSs)) {
		# js <- min(length(codonChangeHGVSs[[i]]), length(aaChangeHGVSs[[i]]))
		js <- length(aaChangeHGVSs[[i]])
		#in case of frameshifts, all other mutations get ignored, so we need to figure out 
		#which one triggered the frameshift(s). In that case, we have fewer aa changes than codon changes listed
		if (length(codonChangeHGVSs[[i]]) != js) {
			#extract positions from aaChanges and  codon changes and find the match
			aapos <- as.integer(extract.groups(aaChangeStrs[[i]],"(\\d+)")[,1])
			ncpos <- as.integer(extract.groups(codonChangeHGVSs[[i]],"(\\d+)")[,1])
			ccpos <- floor((ncpos-1)/3+1)
			#ks are the codon change HGVSs that match the positions for each j
			ks <- lapply(aapos,function(ap) {
			  ds <- abs(ccpos-ap)
			  matches <- which(ds==min(ds))
			  matches[order(ncpos[matches])]
		  })
		} else {
			ks <- as.list(1:js)
		}
		#now iterate over aa changes and index them
		for (j in 1:js) {
			cc <- cisMerge(codonChangeHGVSs[[i]][ks[[j]]])
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
	outfile <- paste0(outDir,"/marginalCounts.csv")
	cat("# MARGINAL VARIANT COUNTS #",
	    "\n# project name:", params$project,
	    "\n# gene name:",params$template$geneName,
	    "\n# tileseqMave version:",as.character(packageVersion("tileseqMave")),
	    "\n# parameter sheet:",normalizePath(paramFile),"\n",
	    file=outfile
	)
	suppressWarnings(
	  write.table(marginalCounts,outfile,sep=",",append=TRUE,row.names=FALSE,qmethod="double")
	)

	logInfo("Done.")

	options(op)
	return(NULL)
}
