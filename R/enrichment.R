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

#' calculate logPhi enrichment ratios
#' 
#' performs error modelling and regularization, filtering and calculates logPhi
#' enrichment values
#' 
#' @param dataDir working data directory
#' @param inDir input directory, defaults to subdirectory with latest timestamp ending in _mut_count.
#' @param outDir output directory, defaults to name of input directory with _QC tag attached.
#' @param paramFile input parameter file. defaults to <dataDir>/parameters.json
#' @param mc.cores number of CPU cores to use in parallel
#' @param srOverride single replicate override
#' @param bnOverride bottleneck filter override
#' @return NULL. Results are written to file.
#' @export
calcEnrichment <- function(dataDir,inDir=NA,outDir=NA,paramFile=paste0(dataDir,"parameters.json"),
	mc.cores=6, srOverride=FALSE, bnOverride=FALSE, useWTfilter=FALSE, nbs=1e4, pessimistic=TRUE,
	useQuorum=FALSE) {

	op <- options(stringsAsFactors=FALSE)

	library(yogitools)
	library(hgvsParseR)
	library(pbmcapply)
	library(optimization)

	if (nbs != 0) {
  	if (nbs < 1e2) {
  	  stop("Must use at least 100 bootstrap samples!")
  	}
  	if (nbs > 1e5) {
  	  logWarn("Large number of bootrap samples requested!",
  	          "This may take a very long time and possibly exceed memory limitations!"
  	  )
  	}
	} else {
    logWarn("Bootstrapping disabled. Using heuristic error propagation instead.")
  }
	
	#make sure data and out dir exist and ends with a "/"
	if (!grepl("/$",dataDir)) {
		dataDir <- paste0(dataDir,"/")
	}
	if (!dir.exists(dataDir)) {
		#we don't use the logger here, assuming that whichever script wraps our function
		#catches the exception and writes to the logger (or uses the log error handler)
		stop("Workspace data folder ",dataDir," does not exist!")
	}
	if (!canRead(paramFile)) {
		stop("Unable to read parameter file!")
	}


	logInfo("Reading parameters from",normalizePath(paramFile))
	params <- withCallingHandlers(
		parseParameters(paramFile,srOverride=srOverride),
		warning=function(w)logWarn(conditionMessage(w))
	)
	#cache function parameter in params object
	params$scoring$useQuorum <- useQuorum
	
	
	regmode <- if (pessimistic) "pessimistic" else "optimistic"
	logInfo("Enrichment function uses the following parameters:")
	logInfo("countThreshold =",params$scoring$countThreshold)
	logInfo("WT filter quantile =",params$scoring$wtQuantile)
	logInfo("pseudoReplicates (pseudo.n) =",params$scoring$pseudo.n)
	logInfo("sdThreshold =",params$scoring$sdThreshold)
	logInfo("cvDeviation =",params$scoring$cvDeviation)
	logInfo("assay direction =",params$assay[["selection"]])
	logInfo("regularization mode: ",regmode)
	
	if (useWTfilter) {
	  logInfo("WT excess filter enabled!")
	}

	if (bnOverride) {
		logWarn(
			"WARNING: Bottleneck override has been enabled!\n",
			" --> The final scores will not be filtered for bottlenecked variants!"
		)
	}

	
	#find counts folder
	if (is.na(inDir)) {
	  latest <- latestSubDir(parentDir=dataDir,pattern="_mut_call$|mut_count$")
	  inDir <- latest[["dir"]]
	} else { #if custom input dir was provided
	  #make sure it exists
	  if (!dir.exists(inDir)) {
	    stop("Input folder ",inDir," does not exist!")
	  }
	}
	#make sure it ends in "/"
	if (!grepl("/$",inDir)) {
	  inDir <- paste0(inDir,"/")
	}
	
	#if no output directory was defined
	if (is.na(outDir)) {
	  #derive one from the input
	  if (grepl("_mut_count/$",inDir)) {
	    outDir <- sub("_mut_count/$","_scores/",inDir)
	  } else {
	    outDir <- sub("/$","_scores/",inDir)
	  }
	} 
	#make sure it ends in "/"
	if (!grepl("/$",outDir)) {
	  outDir <- paste0(outDir,"/")
	}
	#make sure outdir exists
	dir.create(outDir,recursive=TRUE,showWarnings=FALSE)
	
	logInfo("Using input directory",inDir,"and output directory",outDir)
	
	
	#identify selection conditions
	selectConds <- getSelects(params)

	#load the marginal counts	
	marginalCountFile <- paste0(inDir,"/marginalCounts.csv")
	if (!file.exists(marginalCountFile)) {
	  stop("Invalid input folder ",inDir," ! Must contain marginalCounts.csv !")
	}
	marginalCounts <- read.csv(marginalCountFile,comment.char="#")

	#filter out frameshifts and indels
	toAA <- extract.groups(marginalCounts$aaChange,"\\d+(.*)$")
	indelIdx <- which(toAA=="-" | nchar(toAA) > 1)
	silentIdx <- which(marginalCounts$aaChange=="silent")
	if (length(indelIdx) > 0 || length(silentIdx) > 0) {
  	marginalCounts <- marginalCounts[-union(indelIdx,silentIdx),]
	}

	#extract variant positions and assign to regions and tiles
	marginalCounts$position <- as.integer(extract.groups(marginalCounts$codonChange,"(\\d+)"))
	marginalCounts$region <- params$pos2reg(marginalCounts$position)
	marginalCounts$tile <- params$pos2tile(marginalCounts$position)
	
	#load raw depth table
	depthTableFile <- paste0(inDir,"/sampleDepths.csv")
	depthTable <- read.csv(depthTableFile)
	#calculate drops in effective depth
	depthDrops <- calcDepthDrops(marginalCounts,marginalCounts$tile,depthTable,mc.cores)
	rownames(depthDrops) <- marginalCounts$hgvsc
	#TODO: eliminate individual replicates with excessive drops
	#note variants for which all variants have been dropped as filtered
	#adjust degrees of freedem as appropriate!

	# iterate selection conditions -----------------------------------------------
	for (sCond in selectConds) {

	  null2na <- function(x) if (length(x) == 0) NA else x
	  
		#pull up matching nonselect and WT controls
		nCond <- getNonselectFor(sCond,params)
		condQuad <- c(
			select=sCond,
			nonselect=nCond,
			selWT=null2na(getWTControlFor(sCond,params)),
			nonWT=null2na(getWTControlFor(nCond,params))
		)
		condNames <- names(condQuad)

		#Workaround: Since R doesn't like variable names starting with numerals, 
		# we need to adjust any of those cases
		if (any(grepl("^\\d",condQuad))) {
			culprits <- which(grepl("^\\d",condQuad))
			#we add the prefix "X" to the name, just like the table parser does.
			condQuad[culprits] <- paste0("X",condQuad[culprits])
		}

		# iterate time points -----------------------------------------------------
		for (tp in params$timepoints$`Time point name`) {

			logInfo("Processing selection",sCond,"; time point",tp)

			#variable for collecting error model parameters as they become available
			allModelParams <- NULL

			# iterate regions to process separately. --------------
			# regions <- 1:nrow(params$regions)
			regions <- params$regions[,"Region Number"]
			scoreTable <- do.call(rbind,lapply(regions, function(region) {

				#complain if there's no data for this region
				if (!any(marginalCounts$region == region)) {
					logWarn("No data for region ",region)
					return(NULL)
				} else {
					logInfo("Processing region",region)
				}

				#extract relevant count data for this region
			  regIdx <- which(marginalCounts$region == region)
				regionalCounts <- marginalCounts[regIdx,]
				regionalDrops <- depthDrops[regIdx,]

				# means, stdevs and average count calculation for each condition----------
				msc <- do.call(cbind,lapply(
				  condQuad, mean.sd.count, regionalCounts, regionalDrops, tp, params
				))
				
				#determine mutation types
				aac <- regionalCounts$aaChange
				fromAA <- substr(aac,1,1)
				toAA <- substr(aac,nchar(aac),nchar(aac))
				msc$type <- as.vector(mapply(function(from,to) {
				  if (from==to) "synonymous" else if (to=="*") "nonsense" else "missense"
				},fromAA,toAA))

				#for multi-condition maps, filter out rows that don't occur at all in this condition
				unrelated <- which(msc$nonselect.count == 0)
				if (length(unrelated) > 0) {
  				regionalCounts <- regionalCounts[-unrelated,]
  				regionalDrops <- regionalDrops[-unrelated,]
  				msc <- msc[-unrelated,]
				}

				#error modeling only if more than one replicate exists
				# if (params$numReplicates[[sCond]] > 1) {
				if (!srOverride) {

					# error model fitting -------------------------------------
					logInfo("Fitting error models for each tile")
					models <- runModelFits(msc,condNames,regionalCounts$tile, mc.cores=mc.cores)

					#tabulate model parameters
					modelParams <- lapply(models,function(ms) do.call(rbind,lapply(ms,function(ls)do.call(c,ls[-1]))) )
					for (cond in condNames) {
						colnames(modelParams[[cond]]) <- paste0(cond,".",colnames(modelParams[[cond]]))
					}
					modelParams <- do.call(cbind,modelParams)
					allModelParams <<- rbind(allModelParams,modelParams)

					#extract model functions
					modelFunctions <- lapply(models,function(ms) lapply(ms,`[[`,1))

					# regularize stdev of raw frequencies --------------------------------
					logInfo("Performing error regularization")
					msc <- cbind(msc,
						regularizeRaw(msc, condNames, 
							tiles=regionalCounts$tile,
							n=params$numReplicates[condQuad], pseudo.n=params$scoring$pseudo.n, 
							modelFunctions=modelFunctions,pessimistic=pessimistic
						)
					)
				}

				# filtering (count and frequency thresholds met) -----------------------
				logInfo("Filtering...")
				msc$filter <- rawFilter(msc,
				  countThreshold=params$scoring$countThreshold,
          wtq=params$scoring$wtQuantile,
          cvm=params$scoring$cvDeviation,
          srOverride=srOverride,
          useWTfilter=useWTfilter
				)

				
				# log(phi) calculation and error propagation ---------------------------
				logInfo("Calculating logPhi values...")
				if (nbs < 10 || srOverride) {
				  logInfo("Bootstrapping disabled. Defaulting to heuristic error propagation.")
				  msc <- cbind(msc,calcPhi(msc))
				} else {
				  msc <- cbind(msc,calcPhiBootstrap(msc,N=nbs))
				}

				#if this is a negative assay, invert the scores
				if (params$assay[["selection"]] == "Negative") {
					msc$logPhi <- -msc$logPhi
				}

				# Bias correction -----------------------------
				logInfo("Performing bias correction...")
				msc <- cbind(msc,biasCorrection(msc))
				
				
				# # Scaling to synonymous and nonsense medians -----------------------------
				# logInfo("Normalizing...")
				# #check if overrides were provided for the distribution modes
				# normOvr <- getNormOverrides(params,sCond,tp,region)
				# #and apply
				# msc <- cbind(msc,normalizeScores(msc,regionalCounts$aaChange,params$scoring$sdThreshold,normOvr))
				# 
				# # flooring and stderr calculation only where more than one replicate exists ----
				# if (params$numReplicates[[sCond]] > 1) {
				# 	#Flooring and SD adjustment
				# 	msc <- cbind(msc,flooring(msc,params))
				# 
				# 	#add stderr
				# 	msc$se.floored <- msc$sd.floored/sqrt(params$numReplicates[[sCond]])
				# }
				
				#attach labels and return
				return(cbind(regionalCounts[,1:4],msc))
			}))

			# write output ----------------------------------------------
			#export model parameters
			outFile <- paste0(outDir,sCond,"_t",tp,"_errorModel.csv")
			write.csv(as.data.frame(allModelParams),outFile)

			logInfo("Apply formatting...")

			#apply some vanity formatting
			type.idx <- which(colnames(scoreTable)=="type")
			filter.idx <- which(colnames(scoreTable)=="filter")
			new.order <- c(1:4,type.idx,filter.idx,setdiff(5:ncol(scoreTable),c(type.idx,filter.idx)))
			scoreTable <- scoreTable[,new.order]

			logInfo("Writing full table to file.")
			timestamp <- if (exists("latest")) latest[["timeStamp"]] else "?"
			#export to file
			paramHeader <- paste0(
			  "# project name: ", params$project,"\n",
			  "# gene name: ",params$template$geneName,"\n",
			  "# tileseqMave version: ",packageVersion("tileseqMave"),"\n",
			  "# variant call date: ",timestamp,"\n",
			  "# enrichment processing date: ",Sys.time(),"\n",
			  "# condition: ",sCond,"\n",
			  "# time point: ",tp,"\n",
			  "# parameter sheet: ",normalizePath(paramFile),"\n",
			  "# count threshold: ",params$scoring$countThreshold,"\n",
			  "# wt filter quantile: ",params$scoring$wtQuantile,"\n",
			  "# pseudo-replicates: ",params$scoring$pseudo.n,"\n",
			  "# regularization mode: ",regmode,"\n",
			  "# syn/non sd threshold: ",params$scoring$sdThreshold,"\n",
			  "# cv deviation threshold: ",params$scoring$cvDeviation,"\n",
			  "# assay direction: ",params$assay[["selection"]],"\n",
			  "# bootstrap samples: ",nbs,"\n"
			)
			outFile <- paste0(outDir,sCond,"_t",tp,"_enrichment.csv")
			cat("# COMPLETE UNFILTERED ENRICHMENT RATIOS #\n",paramHeader,file=outFile,sep="")
			suppressWarnings(
			  write.table(scoreTable,outFile,sep=",",append=TRUE,row.names=FALSE,qmethod="double")
			)
# 
# 			#configure filter level for export
# 			exportFilter <- if (bnOverride) {
# 				sapply(scoreTable$filter,function(f) is.na(f) || grepl(f=="bottleneck"))
# 			} else {
# 				is.na(scoreTable$filter)
# 			}
# 
# 			#simplified (MaveDB-compatible) format
# 			simpleCols <- if (srOverride) {
# 				c("hgvsc","hgvsp","score","score.sd","score.sd")
# 			} else {
# 				c("hgvsc","hgvsp","score.floored","sd.floored","se.floored")
# 			}
# 			simpleTable <- scoreTable[
# 				exportFilter,
# 				simpleCols
# 			]
# 			colnames(simpleTable) <- c("hgvs_nt","hgvs_pro","score","sd","se")
# 
# 			#export to file
# 			logInfo("Writing simplified table to file.")
# 			outFile <- paste0(outDir,sCond,"_t",tp,"_simple.csv")
# 			cat("# NUCLEOTIDE-LEVEL SCORES #\n",paramHeader,file=outFile,sep="")
# 			write.table(simpleTable,outFile,sep=",",append=TRUE,row.names=FALSE,qmethod="double")
# 
# 			#collapse by amino acid change
# 			logInfo("Collapsing amino acid changes...")
# 			flooredAA <- collapseByAA(scoreTable,params,sCond,"score.floored","sd.floored",srOverride,bnOverride)
# 			unflooredAA <- collapseByAA(scoreTable,params,sCond,"score","score.sd",srOverride,bnOverride)
# 
# 			logInfo("Writing AA-centric table to file.")
# 			outFile <- paste0(outDir,sCond,"_t",tp,"_simple_aa_floored.csv")
# 			cat("# FLOORED AMINO ACID-LEVEL SCORES #\n",paramHeader,file=outFile,sep="")
# 			write.table(flooredAA,outFile,sep=",",append=TRUE,row.names=FALSE,qmethod="double")
# 
# 			outFile <- paste0(outDir,sCond,"_t",tp,"_simple_aa.csv")
# 			cat("# UNFLOORED AMINO ACID-LEVEL SCORES #\n",paramHeader,file=outFile,sep="")
# 			write.table(unflooredAA,outFile,sep=",",append=TRUE,row.names=FALSE,qmethod="double")

		}

	}

	options(op)
	logInfo("Scoring complete.")
	return(NULL)
}

# 
# #' Collapse score table by amino acid changes
# #'
# #' @param scoreTable the full score
# #' @param params the parameter object
# #' @param sCond the current selective condition
# #' @param scoreCol the name of the column containing the scores
# #' @param sdCol the name of the column containing the stdev
# #' @param sdOverride the sdOverride flag
# #' @return a \code{data.frame} containing the collapsed score table
# #' @export
# collapseByAA <- function(scoreTable,params,sCond,scoreCol="score",sdCol="score.sd",srOverride=FALSE,bnOverride=FALSE) {
# 
# 	#configure filter level for export
# 	exportFilter <- if (bnOverride) {
# 		sapply(scoreTable$filter,function(f) is.na(f) || f=="bottleneck")
# 	} else {
# 		is.na(scoreTable$filter)
# 	}
# 
# 	filteredTable <- scoreTable[exportFilter,]
# 
# 	if (!srOverride) {
# 		aaTable <- as.df(tapply(1:nrow(filteredTable),filteredTable$hgvsp, function(is) {
# 			joint <- join.datapoints(
# 				filteredTable[is,scoreCol],
# 				filteredTable[is,sdCol],
# 				rep(params$numReplicates[[sCond]],length(is))
# 			)
# 			list(
# 				hgvs_pro=unique(filteredTable[is,"hgvsp"]),
# 				score=joint[["mj"]],
# 				sd=joint[["sj"]],
# 				df=joint[["dfj"]],
# 				se=joint[["sj"]]/sqrt(joint[["dfj"]])
# 			)
# 		}))
# 	} else {
# 		aaTable <- as.df(tapply(1:nrow(filteredTable),filteredTable$hgvsp, function(is) {
# 			list(
# 				hgvs_pro=unique(filteredTable[is,"hgvsp"]),
# 				score=mean(fin(filteredTable[is,scoreCol])),
# 				sd=NA,
# 				df=length(fin(filteredTable[is,scoreCol]))
# 			)
# 		}))
# 	}
# 	return(aaTable)
# }

#' Error propagation formula for ratios of Random Variables
#'
#' @param m1 the mean of the numerator RV
#' @param m2 the mean of the denominator RV
#' @param sd1 the stdev of the numerator RV
#' @param sd2 the stdev of the denominator RV
#' @param cov12 the covariance of the two variables
ratio.sd <- function(m1,m2,sd1,sd2,cov12=0) {
	abs(m1/m2) * sqrt( (sd1/m1)^2 + (sd2/m2)^2 - 2*cov12/(m1*m2) )
}

#' Error propagation formula for the logarithm of a Random Variable
#'
#' @param m1 the mean of the random variable
#' @param sd1 the standard deviation of the RV
#' @param base the logarithm base (default 10)
log.sd <- function(m1,sd1,base=10) {
	abs(sd1/(m1*log(base)))
}

#' Helper function to bottom-out a vector of numbers, so that no element is smaller than 'bottom'
#'
#' @param xs the input numerical vector
#' @param bottom the minimum to which the vector will be forced.
floor0 <- function(xs,bottom=0) sapply(xs, function(x) if (x < bottom) bottom else x)

#' Baldi & Long's formula
#'
#' @param pseudo.n the number of pseudo-replicates (the weight of the prior)
#' @param n the number of empirical replicates
#' @param model.sd the prior expectation of the standard deviation
#' @param empiric.sd the empircal standard deviation
bnl <- function(pseudo.n,n,model.sd,empiric.sd) {
	sqrt((pseudo.n * model.sd^2 + (n - 1) * empiric.sd^2)/(pseudo.n + n - 2))
}

#' Function to form mean, stdev and average raw count for a given condition
#'
#' @param cond the names of the condition (e.g. select, nonselect)
#' @param regionalCounts the region-specific subset of the marginal counts dataframe
#' @param tp the time point ID
#' @param params the global parameters object
mean.sd.count <- function(cond,regionalCounts,regionalDrops,tp,params) {
  
  #detect dummy conditions (when no WT control exists)
  if (is.na(cond)) {
    #zeroes
    z <- rep(0,nrow(regionalCounts))
    return(data.frame(
      mean=z,sd=z,cv=z,count=z,df=1
    ))
  }
  
  #get the number of replicates for this condition
	nrep <- params$numReplicates[[cond]]
	#validate against applicable time points
	tps <- getTimepointsFor(cond,params)
	if (!(tp %in% tps)) {
	  if (cond %in% getSelects(params)) {
	    stop("Unrecognized time point ",tp," for selection condition ",cond,"! ",
	         "This was not declared in the parameter sheet!")
	  } else {
	    #if this is a nonselect or WT condition, only one time point may exist here, so we default to that one
	    logWarn("No matching timepoint",tp,"found for condition",cond,
	            "! Using timepoint", tps[[1]], "instead! This could lead to incorrect results!"
	    )
	    tp <- tps[[1]]
	  }
	}
	
	maxFactor <- params$scoring$cvDeviation
	maxDrop <- params$varcaller$maxDrop
	
	freqs <- regionalCounts[,sprintf("%s.t%s.rep%d.frequency",cond,tp,1:params$numReplicates[[cond]])]
	counts <- regionalCounts[,sprintf("%s.t%s.rep%d.count",cond,tp,1:params$numReplicates[[cond]])]
	drops <- regionalDrops[,sprintf("%s.t%s.rep%d.effectiveDepth",cond,tp,1:params$numReplicates[[cond]])]
	dropFilter <- drops < maxDrop
	
	# mCount <- if(nrep>1)rowMeans(counts,na.rm=TRUE) else counts
	# cvPois <- 1/sqrt(mCount+.1)
	
	out <- as.df(lapply(1:nrow(regionalCounts), function(i) {
	  cols <- which(dropFilter[i,])
	  dfInit <- length(cols)
	  mCount <- mean(unlist(counts[i,cols]),na.rm=TRUE)
	  fs <- unlist(freqs[i,cols])
	  if (dfInit==0) {
	    list(mean=NA,sd=NA,cv=NA,count=NA,df=0)
	  } else if (dfInit==1) {
	    list(mean=fs,sd=NA,cv=NA,count=mCount,df=1)
	  } else if (dfInit==2||!params$scoring$useQuorum) {
	    m <- mean(fs,na.rm=TRUE)
	    s <- sd(fs,na.rm=TRUE)
	    list(mean=m,sd=s,cv=s/m,count=mCount,df=dfInit)
	  } else {#i.e. dfInit > 2
	    med <- median(fs,na.rm=TRUE)
	    valid <- (fs/med <= maxFactor) & (fs/med >= 1/maxFactor)
	    if (sum(valid,na.rm=TRUE) > 1) {
	      m <- mean(fs[which(valid)])
	      s <- sd(fs[which(valid)])
	      list(mean=m,sd=s,cv=s/m,
	        count=mean(counts[i,cols][which(valid)],na.rm=TRUE),
	        df=sum(valid,na.rm=TRUE)
	      )
	    } else {
	      m <- mean(fs,na.rm=TRUE)
	      s <- sd(fs,na.rm=TRUE)
	      list(mean=m,sd=s,cv=s/m,count=mCount,df=dfInit)
	    }
	  }
	}))
	
	# if (nrep == 1) {
	#   means <- freqs
	#   sds <- rep(NA,length(freqs))
	# } else if (nrep == 2 ||!params$scoring$useQuorum) {
	#   means <- rowMeans(freqs,na.rm=TRUE)
	#   sds <- apply(freqs,1,sd,na.rm=TRUE)
	# } else if (nrep > 2) {
	#   #FIXME: Track degrees of freedom!
	#   majorities <- t(apply(freqs,1,function(fs){
	#     m <- median(fs,na.rm=TRUE)
	#     valid <- (fs/m <= maxFactor) & (fs/m >= 1/maxFactor)
	#     if (sum(valid,na.rm=TRUE) > 1) {
	#       c(
	#         mean=mean(fs[which(valid)]),
	#         sd=sd(fs[which(valid)]),
	#         outliers=sum(!valid,na.rm=TRUE)
	#       )
	#     } else {
	#       c(
	#         mean=mean(fs),
	#         sd=sd(fs),
	#         outliers=length(valid)
	#       )
	#     }
	#   }))
	#   means <- majorities[,"mean"]
	#   sds <- majorities[,"sd"]
	# }
	
	# sds <- if (nrep > 1) apply(freqs,1,sd,na.rm=TRUE) else rep(NA,length(freqs))
	
	# out <- data.frame(
	# 	mean=means,
	# 	sd=sds,
	# 	cv=sds/means,
	# 	# cv=mapply(function(s,m) if(s==0) 0 else s/m,sds,means),
	# 	# count=if(nrep>1)rowMeans(counts,na.rm=TRUE) else counts
	# 	count=mCount
	# 	# csd=apply(counts,1,sd,na.rm=TRUE)
	# )
	return(out)
}

#' Function to apply filtering based on raw counts / frequency
#'
#' @param msc dataframe containing the output of the 'mean.sd.count()' function.
#' @param countThreshold the threshold of raw counts that must be met
#' @param wtq wild-type quantile, the quantile of the WT above which the nonselect or select
#'  have to be to not be filtered (lest they be considered spurious), also the quantile of the 
#'  median deviation-informed distribution around zero, at which we consider "excessive" WT levels.
#' @param cvm coefficient of variation multiplier. Up to how much more than the 
#'  expected CV do we accept as normal?
#' @return a vector listing for each row in the table which (if any) filters apply, otherwise NA.
rawFilter <- function(msc,countThreshold,wtq=0.95,cvm=10,srOverride=FALSE,useWTfilter=FALSE) {
	#if no error estimates are present, pretend it's 0
	if (all(c("nonWT.sd.bayes", "selWT.sd.bayes") %in% colnames(msc))) {
		sd.sWT <- msc$selWT.sd.bayes
		sd.nWT <- msc$nonWT.sd.bayes
	} else if (!all(is.na(msc$nonWT.sd))) {
	  sd.sWT <- msc$selWT.sd
	  sd.nWT <- msc$nonWT.sd
	} else {
	  #last fallback for when there are no replicates
	  sd.sWT <- sd.nWT <- 0
	}
  sQuant <- qnorm(wtq,msc$selWT.mean,sd.sWT)
  nQuant <- qnorm(wtq,msc$nonWT.mean,sd.nWT)
	#calculate nonselect filter using WT control
	nsFilter <- with(msc, 
		nonselect.count < countThreshold | nonselect.mean <= nQuant
		# nonselect.count < countThreshold | nonselect.mean <= nonWT.mean + 3*sd.nWT
	)
	#and select filter using WT control
	sFilter <- with(msc, 
	  #using `<=` instead of `<` to catch zeroes!
		select.mean <= sQuant
		# select.mean <= selWT.mean + 3*sd.sWT
	)
	
	#depth filter
	dfCols <- grep("\\.df$",colnames(msc))
	dFilter <- apply(msc[,dfCols],1,function(dfs) any(dfs < 1))
	
	#aberrant WT count filter
	if (useWTfilter) {
  	medianDeviation <- function(xs) {
  	  m <- median(xs,na.rm=TRUE)
  	  mean(abs(xs-m),na.rm=TRUE)
  	}
  	nonwtCutoff <- qnorm(wtq,0,medianDeviation(msc$nonWT.mean))
  	selwtCutoff <- qnorm(wtq,0,medianDeviation(msc$selWT.mean))
  	wFilter <- with(msc,
  	  nonWT.mean > nonwtCutoff | selWT.mean > selwtCutoff
  	)
	} else {
	  wFilter <- rep(FALSE,nrow(msc))
	}
	
	#determine replicate bottlenecks
	if (!srOverride) {
  	non.cv.poisson <- 1/sqrt(msc[,"nonselect.count"])
  	sel.cv.poisson <- 1/sqrt(msc[,"select.count"]+.1)
  	non.cv <- msc$nonselect.cv
  	sel.cv <- msc$select.cv
  	#remove NAs caused by mean = 0 (their sd is also zero)
  	sel.cv[is.na(sel.cv)] <- 0
	  rFilter <- non.cv > cvm*non.cv.poisson | sel.cv > cvm*sel.cv.poisson
	} else {
	  #if single-replicate override is enabled, we can't use the replicate filter
	  rFilter <- rep(FALSE,length(wFilter))
	}
	
	
	mapply(function(ns,s,rf,w,d) {
	  if (w) "wt_excess"
	  else if (ns) "frequency" 
	  else if (s) "bottleneck:select" 
	  else if (rf) "bottleneck:rep"
	  else if (d) "depth"
	  else NA
	},nsFilter,sFilter,rFilter,wFilter,dFilter)
}

#' Run model fitting on all tiles in all given conditions
#' 
#' @param msc dataframe containing the output of the 'mean.sd.count()' function.
#' @param condNames the condition names
#' @param tiles vector of tile assignments for each row in msc
#' @param mc.cores the number of CPU cores to use in parallel
runModelFits <- function(msc,condNames, tiles, mc.cores) {
	modelSets <- pbmclapply(condNames, function(cond) {
		tapply(1:nrow(msc),tiles,function(i) {
			tryCatch({
				fit.cv.model(msc[i,],cond)
			},error=function(e) {
				list(cv.model=NA,static=NA,additive=NA,multiplicative=NA)
			})
		})
	},mc.cores=mc.cores)
	names(modelSets) <- condNames
	return(modelSets)
}

#' Create a model of the error in for a subset of marginal frequencies
#' 
#' @param  subcores a subset of the count table (for a given tile)
#' @return a list containing the model function (mapping counts to expected CV), and the model parameters
fit.cv.model <- function(subscores,cond,noTryCutoff=10) {
  
  #if there is less than "noTryCutoff" counts each, abort the modeling process
  if (all(subscores[,paste0(cond,".count")] < noTryCutoff)) {
    stop("Aborting modeling due to insufficient counts")
  }
  
	#poisson model of CV
	cv.poisson <- 1/sqrt(subscores[,paste0(cond,".count")])
	#filter out NAs and infinites
	filter <- which(is.na(cv.poisson) | is.infinite(cv.poisson) | 
		is.na(subscores[,paste0(cond,".cv")]) | is.infinite(subscores[,paste0(cond,".cv")]))
	if (length(filter)> 0) {
		subscores <- subscores[-filter,]
		cv.poisson <- cv.poisson[-filter]
	}
	#model function
	log.cv.model <- function(log.cv.poisson,static,additive,multiplicative) {
		sapply(log.cv.poisson, function(x) max(static, multiplicative*x + additive))
	}
	#objective function
	objective <- function(theta) {
		reference <- log10(subscores[,paste0(cond,".cv")])
		prediction <- log.cv.model(log10(cv.poisson),theta[[1]],theta[[2]],theta[[3]])
		sum((prediction-reference)^2,na.rm=TRUE)
	}
	#run optimization
	theta.start <- c(static=-2,additive=0,multiplicative=1)
	z <- optim_nm(objective,start=theta.start)
	theta.optim <- setNames(z$par,names(theta.start))
	#wrap cv model using optimal parameters
	cv.model <- function(count) {
		10^log.cv.model(
			log10(1/sqrt(count)),
			theta.optim[["static"]],
			theta.optim[["additive"]],
			theta.optim[["multiplicative"]]
		)
	}

	# plot(cv.poisson,subscores$nonselect.cv,log="xy")
	# lines(runningFunction(cv.poisson,subscores$nonselect.cv,nbins=20,logScale=TRUE),col=2)
	# abline(0,1,col="green",lty="dashed")
	# lines(seq(0,1,0.01), 10^log.cv.model(
	# 	log10(seq(0,1,0.01)),theta.optim[[1]],theta.optim[[2]],theta.optim[[3]]
	# ),col="blue")
	
	#and return output
	return(c(list(cv.model=cv.model),theta.optim))
}

#' Function to perform Baldi & Long's error regularization
#'
#' @param msc dataframe containing the output of the 'mean.sd.count()' function
#' @param condNames the names of the different count conditions (select, nonselect etc)
#' @param n the number of replicates present for each datapoint in each condition
#' @param pseudo.n the number of pseudo-replicates that detemine the weight of the prior
#'    in the regularization
#' @return a data.frame containing for each condition the possion prior, bayesian regularized, and
#'    pessimistic (maximum) standard deviation.
regularizeRaw <- function(msc,condNames,tiles,n,pseudo.n, modelFunctions,pessimistic=TRUE) {
	#index replicates so we can look them up
	names(n) <- condNames
	#handle missing WT conditions 
	if (any(is.na(n))) {
	  n[is.na(n)] <- 1
	}
	#iterate over conditions and join as columns
	regul <- do.call(cbind,lapply(condNames, function(cond) {
	  #determine pseudocount
	  smallestSelect <- unique(sort(na.omit(msc[,paste0(cond,".mean")])))[[2]]
	  pseudoCount <- 10^floor(log10(smallestSelect))
	  #determine prior expectation
		sd.prior <- sapply(1:nrow(msc), function(i) {
		  #add pseudocounts to get valid priors for zeroes
			count <- msc[i,paste0(cond,".count")]+0.1
			freq <- msc[i,paste0(cond,".mean")]+pseudoCount
			#extract correct model funciton and use it to calculate the prior
			tile <- as.character(tiles[[i]])
			modelfun <- modelFunctions[[cond]][[tile]]
			if (is.function(modelfun)) {
				cv.prior <- modelfun(count)
			} else {
				#if model fitting failed, use simple poisson model
				cv.prior <- 1/sqrt(count)
			}
			return(cv.prior * freq)
		})
		# sd.prior[which(msc[,paste0(cond,".mean")] == 0)] <- 0

		sd.empiric <- msc[,paste0(cond,".sd")]
		if (pessimistic) {
		  sd.bayes <- mapply(max,sd.prior,sd.empiric,na.rm=TRUE)
		} else {
		  sd.bayes <- bnl(pseudo.n,n[[cond]],sd.prior,sd.empiric)
		}
		out <- cbind(sd.prior=sd.prior,sd.bayes=sd.bayes)
		colnames(out) <- paste0(cond,c(".sd.prior",".sd.bayes"))
		out
	}))
	return(regul)
}

# sample from beta distribution using mean and sd
rbeta2 <- function(n,m,s) {
  maxVar <- m*(1-m)
  if (s^2 >= maxVar) {
    warning("Variance too large!")
    s <- sqrt(maxVar)
  }
  nu <- (m*(1-m))/(s^2)-1
  a <- m*nu
  b <- (1-m)*nu
  rbeta(n,a,b)
}

calcPhiBootstrap <- function(msc,N=10000) {
  
  #determine pseudocounts
  smallestSelect <- unique(sort(na.omit(msc$select.mean)))[[2]]
  pseudoCount <- 10^floor(log10(smallestSelect))
  
  #calculate phi
  select <- with(msc,floor0(select.mean-selWT.mean))
  nonselect <- with(msc,floor0(nonselect.mean-nonWT.mean))
  phi <- (select+pseudoCount)/nonselect
  logPhi <- log10(phi)
  
  logInfo("Bootstrapping for error propagation...")
  boot <- pbmclapply(1:N, function(i) {
    selectSam <- with(msc,rnorm(
      length(select.mean),select.mean,
      select.sd.bayes/sqrt(select.df)
    ))
    selWTSam <- with(msc,rnorm(
      length(selWT.mean),selWT.mean,
      selWT.sd.bayes/sqrt(selWT.df)
    ))
    nonselectSam <- with(msc,rnorm(
      length(nonselect.mean),nonselect.mean,
      nonselect.sd.bayes/sqrt(nonselect.df)
    ))
    nonWTSam <- with(msc,rnorm(
      length(nonWT.mean),nonWT.mean,
      nonWT.sd.bayes/sqrt(nonWT.df)
    ))
    
    select <- floor0(selectSam-selWTSam)
    nonselect <-floor0(nonselectSam-nonWTSam)
    phi <- (select+pseudoCount)/nonselect
    logPhi <- log10(phi)
    return(cbind(phi=phi,logPhi=logPhi))
  })
  sds <- apply(do.call(zbind,boot),c(1,2),function(x)sd(fin(x)))
  
  #estimate joint degrees of freedom
  jdf <- rowSums(msc[,c("select.df","nonselect.df","selWT.df","nonWT.df")])-4
  jdf <- sapply(jdf,max,1)
  
  return(data.frame(phi=phi,phi.se=sds[,"phi"],logPhi=logPhi,logPhi.se=sds[,"logPhi"],df=jdf))
  
}

#' Calculate enrichment ratios (phi) and perform error propagation
#'
#' @param msc data.frame containing the output of the 'mean.sd.count()' function
#' @return a data.frame containing phi, log(phi) and their respective stdevs
calcPhi <- function(msc) {
	select <- with(msc,floor0(select.mean-selWT.mean))
	nonselect <- with(msc,floor0(nonselect.mean-nonWT.mean))
	#determine pseudocounts
	smallestSelect <- unique(sort(na.omit(msc$select.mean)))[[2]]
	pseudoCount <- 10^floor(log10(smallestSelect))
	#calculate phi
	phi <- (select+pseudoCount)/nonselect

	#and its stdev
	phi.sd <- if ("select.sd.bayes" %in% colnames(msc)) {
		with(msc,ratio.sd(
			m1=select+pseudoCount,
			m2=nonselect,
			sd1=sqrt(select.sd.bayes^2+selWT.sd.bayes^2),
			sd2=sqrt(nonselect.sd.bayes^2+nonWT.sd.bayes^2)
		))
	} else rep(NA,length(phi))
	#and stderr via joint df
	jdf <- rowSums(msc[,c("select.df","nonselect.df","selWT.df","nonWT.df")])-4
	jdf <- sapply(jdf,max,1)
	phi.se <- phi.sd/sqrt(jdf)


	#as well as logPhi and its stdev
	logPhi <- log10(phi)
	logPhi.sd <- with(msc,log.sd(phi,phi.sd))

	#FIXME: for now, cap logPhi at 2*95% quantile to get around heuristic problems
	#ultimately, this should be fixed with bootstrapping
	lpsCap <- 2*quantile(logPhi.sd,0.95,na.rm=TRUE)
	logPhi.sd[logPhi.sd > lpsCap] <- lpsCap
	
	logPhi.se <- logPhi.sd/sqrt(jdf)

	return(data.frame(phi=phi,phi.se=phi.se,logPhi=logPhi,logPhi.se=logPhi.se,df=jdf))
}

#' Run read-frequecy bias correction on logPhi
#'
#' @param msc the enrichment table
#'
#' @return the bias-corrected enrichment (bce) and its stddev
#' @export
biasCorrection <- function(msc) {
  
  #exclude non-depleted nonsense variants
  nonsenseF <- msc[with(msc,type=="nonsense" & is.na(filter) & logPhi < 0),]
  synonymousF <- msc[with(msc,type=="synonymous" & is.na(filter)),]
  
  if (nrow(nonsenseF) == 0 || nrow(synonymousF)==0) {
    logWarn("Unable to perform bias correction, due to insufficient reference data!")
    return(data.frame(bce=rep(NA,nrow(msc)),bce.se=rep(NA,nrow(msc))))
  }
  
  #simple linear regression
  zns <- with(nonsenseF, lm(logPhi~log10(nonselect.mean)) )
  zsyn <- with(synonymousF, lm(logPhi~log10(nonselect.mean)) )
  
  #calculate pivot points for each variant based on their read frequency (=nonselect.mean)
  nsmean <- msc[,"nonselect.mean",drop=FALSE]
  msc$esyn <- predict.lm(zsyn,nsmean)
  msc$enon <- predict.lm(zns,nsmean)
  msc$ediff <- with(msc,floor0(esyn-enon))
  
  bces <- as.df(lapply(1:nrow(msc), function(i) with(msc[i,],{
    c(
      bce=(logPhi-enon)/ediff,
      bce.se=logPhi.se/abs(ediff)
    )
  })))
  return(bces)
}
