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

#' run scoring function
#' 
#' @param dataDir working data directory
#' @param paramFile input parameter file. defaults to <dataDir>/parameters.json
#' @return NULL. Results are written to file.
#' @export
scoring <- function(dataDir,paramFile=paste0(dataDir,"parameters.json"),logger=NULL,mc.cores=6,srOverride=FALSE) {

	op <- options(stringsAsFactors=FALSE)

	library(yogitools)
	library(hgvsParseR)
	library(pbmcapply)
	library(optimization)

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
	params <- withCallingHandlers(
		parseParameters(paramFile,srOverride=srOverride),
		warning=function(w)logWarn(conditionMessage(w))
	)
	

	logInfo("Scoring function uses the following parameters:")
	logInfo("countThreshold =",params$scoring$countThreshold)
	logInfo("pseudoReplicates (pseudo.n) =",params$scoring$pseudo.n)
	logInfo("sdThreshold =",params$scoring$sdThreshold)

	#find counts folder
	# subDirs <- list.dirs(dataDir,recursive=FALSE)
	# countDirs <- subDirs[grepl("_mut_call$",subDirs)]
	# if (length(countDirs) == 0) {
	# 	stop("No mutation call output found!")
	# }
	# latestCountDir <- sort(countDirs,decreasing=TRUE)[[1]]
	# logInfo("Selecting latest data directory: ",latestCountDir)
	# #extract time stamp
	# timeStamp <- extract.groups(latestCountDir,"/(\\d{4}-\\d{2}-\\d{2}-\\d{2}-\\d{2}-\\d{2})")

	#find counts folder
	latest <- latestSubDir(parentDir=dataDir,pattern="_mut_call$|mut_count$")

	#create a matching output directory
	outDir <- paste0(dataDir,latest[["label"]],latest[["timeStamp"]],"_scores/")
	dir.create(outDir,recursive=TRUE,showWarnings=FALSE)

	#identify selection conditions
	selectConds <- getSelects(params)

	#load the marginal counts	
	marginalCountFile <- paste0(latest[["dir"]],"/marginalCounts.csv")
	marginalCounts <- read.csv(marginalCountFile)

	#filter out frameshifts and indels
	toAA <- extract.groups(marginalCounts$aaChange,"\\d+(.*)$")
	indelIdx <- which(toAA=="-" | nchar(toAA) > 1)
	marginalCounts <- marginalCounts[-indelIdx,]

	#extract variant positions and assign to regions and tiles
	marginalCounts$position <- as.integer(extract.groups(marginalCounts$codonChange,"(\\d+)"))
	marginalCounts$region <- sapply(marginalCounts$position, function(pos) {
		which(params$regions[,"Start AA"] <= pos & params$regions[,"End AA"] >= pos)
	})
	marginalCounts$tile <- sapply(marginalCounts$position, function(pos) {
		which(params$tiles[,"Start AA"] <= pos & params$tiles[,"End AA"] >= pos)
	})

	#iterate over (possibly multiple different) selection conditons
	for (sCond in selectConds) {

		#pull up matching nonselect and WT controls
		nCond <- getNonselectFor(sCond,params)
		condQuad <- c(
			select=sCond,
			nonselect=nCond,
			selWT=getWTControlFor(sCond,params),
			nonWT=getWTControlFor(nCond,params)
		)
		condNames <- names(condQuad)

		#Workaround: Since R doesn't like variable names starting with numerals, 
		# we need to adjust any of those cases
		if (any(grepl("^\\d",condQuad))) {
			culprits <- which(grepl("^\\d",condQuad))
			#we add the prefix "X" to the name, just like the table parser does.
			condQuad[culprits] <- paste0("X",condQuad[culprits])
		}

		#iterate over time points
		for (tp in params$timepoints$`Time point name`) {

			logInfo("Processing selection",sCond,"; time point",tp)

			#variable for collecting error model parameters as they become available
			allModelParams <- NULL

			#iterate over mutagenesis regions and process separately.
			regions <- 1:nrow(params$regions)
			scoreTable <- do.call(rbind,lapply(regions, function(region) {

				#complain if there's no data for this region
				if (!any(marginalCounts$region == region)) {
					logWarn("No data for region ",region)
					return(NULL)
				} else {
					logInfo("Processing region",region)
				}

				#extract relevant count data for this region
				regionalCounts <- marginalCounts[which(marginalCounts$region == region),]

				#Calculate means, stdevs and average count for each condition
				msc <- do.call(cbind,lapply(condQuad, mean.sd.count, regionalCounts, tp, params))

				#error modeling only if more than one replicate exists
				if (params$numReplicates[[sCond]] > 1) {

					#fit tile-specific error models
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

					#Regularize stdev of raw frequencies
					logInfo("Performing error regularization")
					msc <- cbind(msc,
						regularizeRaw(msc, condNames, 
							tiles=regionalCounts$tile,
							n=params$numReplicates[condQuad], pseudo.n=params$scoring$pseudo.n, 
							modelFunctions=modelFunctions
						)
					)

					#Apply raw filter (count and frequency thresholds met)
					logInfo("Filtering...")
					msc$filter <- rawFilter(msc,params$scoring$countThreshold)

				} else {
					#quick-and-dirty filter for single-replicate data
					msc$filter <- sapply(msc$nonselect.count < params$scoring$countThreshold | msc$select.count < 1, ifelse, "count", NA) 
				}
				#Calculate enrichment ratios (phi) and propagate error
				logInfo("Scoring...")
				msc <- cbind(msc,calcPhi(msc))

				#if this is a negative assay, invert the scores
				if (params$assay[["selection"]] == "Negative") {
					msc$logPhi <- -msc$logPhi
				}

				#Normalize to synonymous and nonsense medians
				logInfo("Normalizing...")
				#check if overrides were provided for the distribution modes
				normOvr <- getNormOverrides(params,sCond,tp,region)
				if (any(!is.na(normOvr))) {
					logWarn("Applying normalization override!")
				}
				#and apply
				msc <- cbind(msc,normalizeScores(msc,regionalCounts$aaChange,params$scoring$sdThreshold,normOvr))

				#flooring and stderr calculation only where more than one replicate exists
				if (params$numReplicates[[sCond]] > 1) {
					#Flooring and SD adjustment
					msc <- cbind(msc,flooring(msc,params))

					#add stderr
					msc$se.floored <- msc$sd.floored/sqrt(params$numReplicates[[sCond]])
				}
				
				#attach labels and return
				return(cbind(regionalCounts[,1:4],msc))
			}))

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
			#export to file
			outFile <- paste0(outDir,sCond,"_t",tp,"_complete.csv")
			write.csv(scoreTable,outFile,row.names=FALSE)

			#simplified (MaveDB-compatible) format
			simpleCols <- if (srOverride) {
				c("hgvsc","hgvsp","score","score.sd","score.sd")
			} else {
				c("hgvsc","hgvsp","score.floored","sd.floored","se.floored")
			}
			simpleTable <- scoreTable[
				is.na(scoreTable$filter),
				simpleCols
			]
			colnames(simpleTable) <- c("hgvs_nt","hgvs_pro","score","sd","se")

			#export to file
			logInfo("Writing simplified table to file.")
			outFile <- paste0(outDir,sCond,"_t",tp,"_simple.csv")
			write.csv(simpleTable,outFile,row.names=FALSE)

			#collapse by amino acid change
			logInfo("Collapsing amino acid changes...")
			if (!srOverride) {
				aaTable <- as.df(with(simpleTable,tapply(1:length(hgvs_pro),hgvs_pro, function(is) {
					joint <- join.datapoints(
						score[is],
						sd[is],
						rep(params$numReplicates[[sCond]],length(is))
					)
					list(
						hgvs_pro=unique(hgvs_pro[is]),
						score=joint[["mj"]],
						sd=joint[["sj"]],
						df=joint[["dfj"]],
						se=joint[["sj"]]/sqrt(joint[["dfj"]])
					)
				})))
			} else {
				aaTable <- as.df(with(simpleTable,tapply(1:length(hgvs_pro),hgvs_pro, function(is) {
					list(
						hgvs_pro=unique(hgvs_pro[is]),
						score=mean(fin(score[is])),
						sd=NA,
						df=length(fin(score[is]))
					)
				})))
			}

			logInfo("Writing AA-centric table to file.")
			outFile <- paste0(outDir,sCond,"_t",tp,"_simple_aa.csv")
			write.csv(aaTable,outFile,row.names=FALSE)

		}

	}

	options(op)
	logInfo("Scoring complete.")
	return(NULL)
}

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
mean.sd.count <- function(cond,regionalCounts,tp,params) {
	nrep <- params$numReplicates[[cond]]
	freqs <- regionalCounts[,sprintf("%s.t%s.rep%d.frequency",cond,tp,1:params$numReplicates[[cond]])]
	counts <- regionalCounts[,sprintf("%s.t%s.rep%d.count",cond,tp,1:params$numReplicates[[cond]])]
	means <- if (nrep > 1) rowMeans(freqs,na.rm=TRUE) else freqs
	sds <- if (nrep > 1) apply(freqs,1,sd,na.rm=TRUE) else rep(NA,length(freqs))
	out <- data.frame(
		mean=means,
		sd=sds,
		cv=sds/means,
		count=if(nrep>1)rowMeans(counts,na.rm=TRUE) else counts
		# csd=apply(counts,1,sd,na.rm=TRUE)
	)
	return(out)
}

#' Function to apply filtering based on raw counts / frequency
#'
#' @param msc dataframe containing the output of the 'mean.sd.count()' function.
#' @param countThreshold the threshold of raw counts that must be met
#' @return a vector listing for each row in the table which (if any) filters apply, otherwise NA.
rawFilter <- function(msc,countThreshold) {
	#apply nonselect filter
	nsFilter <- with(msc, 
		nonselect.count < countThreshold | nonselect.mean <= nonWT.mean + 3*nonWT.sd.bayes
	)
	sFilter <- with(msc, 
		select.mean <= selWT.mean + 3*selWT.sd.bayes
	)
	mapply(function(ns,s) {
		if (ns) "frequency" else if (s) "bottleneck" else NA
	},nsFilter,sFilter)
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
			tryCatch(fit.cv.model(msc[i,],cond),error=function(e) {
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
fit.cv.model <- function(subscores,cond) {
	#poisson model of CV
	cv.poisson <- 1/sqrt(subscores[,paste0(cond,".count")])
	#filter out NAs and infinites
	filter <- which(is.na(cv.poisson) | is.infinite(cv.poisson) | 
		is.na(subscores[,paste0(cond,".cv")]) | is.infinite(subscores[,paste0(cond,".cv")]))
	#model function
	log.cv.model <- function(log.cv.poisson,static,additive,multiplicative) {
		sapply(log.cv.poisson, function(x) max(static, multiplicative*x + additive))
	}
	#objective function
	objective <- function(theta) {
		reference <- log10(subscores[-filter,paste0(cond,".cv")])
		prediction <- log.cv.model(log10(cv.poisson[-filter]),theta[[1]],theta[[2]],theta[[3]])
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
regularizeRaw <- function(msc,condNames,tiles,n,pseudo.n, modelFunctions) {
	#index replicates so we can look them up
	names(n) <- condNames
	#iterate over conditions and join as columns
	regul <- do.call(cbind,lapply(condNames, function(cond) {
		sd.prior <- sapply(1:nrow(msc), function(i) {
			count <- msc[i,paste0(cond,".count")]
			tile <- as.character(tiles[[i]])
			modelfun <- modelFunctions[[cond]][[tile]]
			if (is.function(modelfun)) {
				cv.prior <- modelfun(count)
			} else {
				#if model fitting failed, use simple poisson model
				cv.prior <- 1/sqrt(count)
			}
			return(cv.prior * msc[i,paste0(cond,".mean")])
		})
		sd.prior[which(msc[,paste0(cond,".mean")] == 0)] <- 0

		sd.empiric <- msc[,paste0(cond,".sd")]
		sd.bayes <- bnl(pseudo.n,n[[cond]],sd.prior,sd.empiric)
		out <- cbind(sd.prior=sd.prior,sd.bayes=sd.bayes)
		colnames(out) <- paste0(cond,c(".sd.prior",".sd.bayes"))
		out
	}))
	return(regul)
}


#' Calculate enrichment ratios (phi) and perform error propagation
#'
#' @param msc data.frame containing the output of the 'mean.sd.count()' function
#' @return a data.frame containing phi, log(phi) and their respective stdevs
calcPhi <- function(msc) {
	#calculate phi
	phi <- with(msc,{
		floor0(select.mean-selWT.mean)/floor0(nonselect.mean-nonWT.mean)
	})
	#and its stdev
	phi.sd <- if ("select.sd.bayes" %in% colnames(msc)) {
		with(msc,ratio.sd(
			m1=floor0(select.mean-selWT.mean),
			m2=floor0(nonselect.mean-nonWT.mean),
			sd1=sqrt(select.sd.bayes^2+selWT.sd.bayes^2),
			sd2=sqrt(nonselect.sd.bayes^2+nonWT.sd.bayes^2)
		))
	} else rep(NA,length(phi))
	#as well as logPhi and its stdev
	logPhi <- log10(phi)
	logPhi.sd <- with(msc,log.sd(phi,phi.sd))
	return(data.frame(phi=phi,phi.sd=phi.sd,logPhi=logPhi,logPhi.sd=logPhi.sd))
}

#' Helper function to extract normalization override values from the parameter sheet
#' if they exist
#' @param params the parameter sheet
#' @param sCond the current selective condition ID
#' @param tp the time point ID
#' @param region the region ID
#' @return a vector with syn and stop overrides, or \code{NA} if they weren't defined.
getNormOverrides <- function(params,sCond,tp,region) {
	pullValue <- function(type) {
		with(params$normalization,{
			row <- which(Condition == sCond & `Time point` == tp & Region == region & Type == type)
			if (length(row)==0) NA else Value[row]
		})
	}
	if (length(params$normalization) == 0 || nrow(params$normalization) == 0) {
		return(c(syn=NA,non=NA))
	} else {
		return(c(syn=pullValue("synonymous"),non=pullValue("nonsense")))
	}
}

#' Calculate scores by normalizing logPhi values to syn/nonsense medians
#' @param msc data.frame containing the enrichment ratios (phi) and log(phi)
#' @param aac the vector of corresponding amino acid changes
#' @param sdThreshold stdev threshold for finding the syn/stop means
#' @return data.frame with variant type, score and stdev
normalizeScores <- function(msc,aac,sdThreshold,overrides=c(syn=NA,non=NA)) {

	#determine mutation types
	fromAA <- substr(aac,1,1)
	toAA <- substr(aac,nchar(aac),nchar(aac))
	msc$type <- as.vector(mapply(function(from,to) {
		if (from==to) "synonymous" else if (to=="*") "nonsense" else "missense"
	},fromAA,toAA))

	#apply filter
	mscFiltered <- if (!all(is.na(msc$logPhi.sd))) {
		with(msc,msc[which(is.na(filter) & logPhi.sd < sdThreshold),])
	} else {
		with(msc,msc[which(is.na(filter)),])
	}

	#calculate medians
	synonymousMedian <- with(mscFiltered,median(
		logPhi[which(type == "synonymous")]
	,na.rm=TRUE))
	nonsenseMedian <- with(mscFiltered,median(
		logPhi[which(type == "nonsense")]
	,na.rm=TRUE))

	#apply any potential overrides
	if (!is.na(overrides[["syn"]])) {
		synonymousMedian <- overrides[["syn"]]
	}
	if (!is.na(overrides[["non"]])) {
		nonsenseMedian <- overrides[["non"]]
	}

	#use medians to normalize
	score <- (msc$logPhi - nonsenseMedian) / (synonymousMedian - nonsenseMedian)
	score.sd <- msc$logPhi.sd / (synonymousMedian - nonsenseMedian)

	return(data.frame(type=msc$type,score=score,score.sd=score.sd))
}

#' Apply flooring and SD adjustment to negative scores
#' 
#' This assumes that no cell can be "deader than dead", and that negative scores are just due to 
#' measurement error. From this follows that the more negative the score, the more certain we 
#' can be that the cell is dead. We can therefore adjust the SD to reflect this, by setting 
#' the score to zero and shifting the SD such that the implied probability of being above 0 
#' remains the same.
#' 
#' @param msc the score matrix
#' @param params the global parameter object
#' @return a matrix containing the floored scores and SDs
flooring <- function(msc,params) {
	#depending on the selection type of the assay, we're either flooring at the top or bottom.
	switch(params$assay[["selection"]],
		Positive={
			#the target null-like score towards which we will shift these values
			targetScore <- 0
			#the quantile for which we want to keep the p-value fixed
			quantile <- 1
			#the row numbers containing the cases to be fixed
			toFix <- which(msc$score < targetScore)
		},
		Negative={
			#if we're dealing with an inverse assay, we have to apply a ceiling instead of flooring
			#the target functional (but dead) score towards which we will shift these values
			targetScore <- 1
			#the quantile for which we want to keep the p-value fixed
			quantile <- 0
			#the row numbers containing the cases to be fixed
			toFix <- which(msc$score > targetScore)
		}
	)
	score.floored <- msc$score
	sd.floored <- msc$score.sd
	score.floored[toFix] <- targetScore
	#the equivalent sds of a normal distribution with the target mean based on the above area
	sd.floored[toFix] <- with(msc[toFix,], score.sd*(quantile-targetScore)/(quantile-score))

	return(cbind(score.floored=score.floored,sd.floored=sd.floored))
}

#' join multiple datapoints weighted by stdev
#' @param ms data means
#' @param sds data stdevs
#' @param degrees of freedom for each data point
#' @return a vector containing the joint mean and joint stdev
join.datapoints <- function(ms,sds,dfs) {
	#weights
	ws <- (1/sds)/sum(1/sds)
	#weighted mean
	mj <- sum(ws*ms)
	#weighted joint variance
	vj <- sum(ws*(sds^2+ms^2)) -mj^2
	#return output
	return(c(mj=mj,sj=sqrt(vj),dfj=sum(dfs)))
}
