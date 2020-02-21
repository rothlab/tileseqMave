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
scoring <- function(dataDir,paramFile=paste0(dataDir,"parameters.json"),logger=NULL,mc.cores=6, 
	countThreshold=10,pseudo.n=8,sdThreshold=0.3) {

	op <- options(stringsAsFactors=FALSE)

	library(yogitools)
	library(hgvsParseR)
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


	logInfo("Reading parameters")
	params <- parseParameters(paramFile)


	#find counts folder
	subDirs <- list.dirs(dataDir,recursive=FALSE)
	countDirs <- subDirs[grepl("_mut_call$",subDirs)]
	if (length(countDirs) == 0) {
		stop("No mutation call output found!")
	}
	latestCountDir <- sort(countDirs,decreasing=TRUE)[[1]]
	logInfo("Selecting latest data directory: ",latestCountDir)
	#extract time stamp
	timeStamp <- extract.groups(latestCountDir,"/(\\d{4}-\\d{2}-\\d{2}-\\d{2}-\\d{2}-\\d{2})")

	#create a matching output directory
	outDir <- paste0(dataDir,timeStamp,"_scores/")
	dir.create(outDir,recursive=TRUE,showWarnings=FALSE)

	#identify selection conditions
	selectConds <- getSelects(params)

	#load the marginal counts	
	marginalCountFile <- paste0(latestCountDir,"/marginalCounts.csv")
	marginalCounts <- read.csv(marginalCountFile)

	#filter out frameshifts and indels
	toAA <- extract.groups(marginalCounts$aaChange,"\\d+(.*)$")
	indelIdx <- which(toAA=="-" | nchar(toAA) > 1)
	marginalCounts <- marginalCounts[-indelIdx,]

	#extract variant positions and assign to regions
	marginalCounts$position <- as.integer(extract.groups(marginalCounts$codonChange,"(\\d+)"))
	marginalCounts$region <- sapply(marginalCounts$position, function(pos) {
		which(params$regions[,"Start AA"] <= pos & params$regions[,"End AA"] >= pos)
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

		#iterate over time points
		for (tp in params$timepoints$`Time point name`) {

			logInfo("Processing selection",sCond,"; time point",tp)

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
				
				#Regularize stdev of raw frequencies
				msc <- cbind(msc,
					regularizeRaw(msc,condNames=names(condQuad),n=params$numReplicates,pseudo.n=pseudo.n)
				)

				#Apply raw filter (count and frequency thresholds met)
				msc$filter <- rawFilter(msc,countThreshold)

				#Calculate enrichment ratios (phi) and propagate error
				msc <- cbind(msc,calcPhi(msc))

				#if this is a negative assay, invert the scores
				if (params$assay[["selection"]] == "Negative") {
					msc$logPhi <- -msc$logPhi
				}

				#Normalize to synonymous and nonsense medians
				msc <- cbind(msc,normalizeScores(msc,regionalCounts$aaChange,sdThreshold))

				#Flooring and SD adjustment
				msc <- cbind(msc,flooring(msc,params))

				#add stderr
				msc$se.floored <- msc$sd.floored/sqrt(params$numReplicates)
				
				#attach labels and return
				return(cbind(regionalCounts[,1:4],msc))
			}))

			#apply some vanity formatting
			type.idx <- which(colnames(scoreTable)=="type")
			filter.idx <- which(colnames(scoreTable)=="filter")
			new.order <- c(1:4,type.idx,filter.idx,setdiff(5:ncol(scoreTable),c(type.idx,filter.idx)))
			scoreTable <- scoreTable[,new.order]

			#export to file
			outFile <- paste0(outDir,sCond,"_t",tp,"_complete.csv")
			write.csv(scoreTable,outFile,row.names=FALSE)

			#simplified (MaveDB-compatible) format
			simpleTable <- scoreTable[
				is.na(scoreTable$filter),
				c("hgvsc","hgvsp","score.floored","sd.floored","se.floored")
			]
			colnames(simpleTable) <- c("hgvs_nt","hgvs_pro","score","sd","se")

			#export to file
			outFile <- paste0(outDir,sCond,"_t",tp,"_simple.csv")
			write.csv(simpleTable,outFile,row.names=FALSE)

			#collapse by amino acid change
			aaTable <- as.df(with(simpleTable,tapply(1:length(hgvs_pro),hgvs_pro, function(is) {
				joint <- join.datapoints(
					score[is],
					sd[is],
					rep(params$numReplicates,length(is))
				)
				list(
					hgvs_pro=unique(hgvs_pro[is]),
					score=joint[["mj"]],
					sd=joint[["sj"]],
					df=joint[["dfj"]],
					se=joint[["sj"]]/sqrt(joint[["dfj"]])
				)
			})))

			outFile <- paste0(outDir,sCond,"_t",tp,"_simple_aa.csv")
			write.csv(aaTable,outFile,row.names=FALSE)

		}

	}


	options(op)
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
	freqs <- regionalCounts[,sprintf("%s.t%s.rep%d.frequency",cond,tp,1:params$numReplicates)]
	counts <- regionalCounts[,sprintf("%s.t%s.rep%d.count",cond,tp,1:params$numReplicates)]
	means <- rowMeans(freqs,na.rm=TRUE)
	sds <- apply(freqs,1,sd,na.rm=TRUE)
	out <- data.frame(
		mean=means,
		sd=sds,
		cv=sds/means,
		count=rowMeans(counts,na.rm=TRUE)
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

#' Function to perform Baldi & Long's error regularization
#'
#' @param msc dataframe containing the output of the 'mean.sd.count()' function
#' @param condNames the names of the different count conditions (select, nonselect etc)
#' @param n the number of replicates present for each datapoint
#' @param pseudo.n the number of pseudo-replicates that detemine the weight of the prior
#'    in the regularization
#' @return a data.frame containing for each condition the possion prior, bayesian regularized, and
#'    pessimistic (maximum) standard deviation.
regularizeRaw <- function(msc,condNames,n,pseudo.n) {
	regul <- do.call(cbind,lapply(condNames, function(cond) {
		sd.empiric <- msc[,paste0(cond,".sd")]
		cv.poisson <- 1/sqrt(msc[,paste0(cond,".count")])
		cv.poisson[which(msc[,paste0(cond,".mean")] == 0)] <- 0
		sd.poisson <- cv.poisson*msc[,paste0(cond,".mean")]
		sd.bayes <- bnl(pseudo.n,n,sd.poisson,sd.empiric)
		sd.pessim <- mapply(max,sd.empiric,sd.poisson)
		out <- cbind(sd.poisson=sd.poisson,sd.bayes=sd.bayes,sd.pessim=sd.pessim)
		colnames(out) <- paste0(cond,c(".sd.poisson",".sd.bayes",".sd.pessim"))
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
	phi.sd <- with(msc,ratio.sd(
		m1=floor0(select.mean-selWT.mean),
		m2=floor0(nonselect.mean-nonWT.mean),
		sd1=sqrt(select.sd.bayes^2+selWT.sd.bayes^2),
		sd2=sqrt(nonselect.sd.bayes^2+nonWT.sd.bayes^2)
	))
	#as well as logPhi and its stdev
	logPhi <- log10(phi)
	logPhi.sd <- with(msc,log.sd(phi,phi.sd))
	return(data.frame(phi=phi,phi.sd=phi.sd,logPhi=logPhi,logPhi.sd=logPhi.sd))
}

#' Calculate scores by normalizing logPhi values to syn/nonsense medians
#' @param msc data.frame containing the enrichment ratios (phi) and log(phi)
#' @param aac the vector of corresponding amino acid changes
#' @param sdThreshold the maximum stdev of logPhi to filter for before calculating medians
#' @return data.frame with variant type, score and stdev
normalizeScores <- function(msc,aac,sdThreshold) {
	#determine mutation types
	fromAA <- substr(aac,1,1)
	toAA <- substr(aac,nchar(aac),nchar(aac))
	type <- as.vector(mapply(function(from,to) {
		if (from==to) "synonymous" else if (to=="*") "nonsense" else "missense"
	},fromAA,toAA))
	#calculate medians
	nonsenseMedian <- with(msc,median(
		logPhi[which(is.na(filter) & type == "nonsense" & logPhi.sd < sdThreshold)]
	,na.rm=TRUE))
	synonymousMedian <- with(msc,median(
		logPhi[which(is.na(filter) & type == "synonymous" & logPhi.sd < sdThreshold)]
	,na.rm=TRUE))
	#use medians to normalize
	score <- (msc$logPhi - nonsenseMedian) / (synonymousMedian - nonsenseMedian)
	score.sd <- msc$logPhi.sd / (synonymousMedian - nonsenseMedian)

	return(data.frame(type=type,score=score,score.sd=score.sd))
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
