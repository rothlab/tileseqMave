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
			regionalResults <- do.call(rbind,lapply(regions, function(region) {

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
				#(This should be moved into its own function)
				msc <- cbind(msc,calcPhi(msc))

				#TODO: Introduce flag for inverse selection assays!!

				#Normalize to synonymous and nonsense medians
				#(This should be moved into its own function)
				aac <- regionalCounts$aaChange
				msc <- cbind(msc,normalizeScores(msc,aac,sdThreshold))
				
				#TODO: Flooring and SD adjustment
				return(cbind(regionalCounts[,1:4],msc))

			}))

		}

	}


	options(op)
	return(NULL)
}

#Calculate enrichment ratios (phi) and perform error propagation
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

#calculate scores by normalizing logPhi values to syn/nonsense medians
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

	return(cbind(type=type,score=score,score.sd=score.sd))
}


#error propagation for ratios
ratio.sd <- function(m1,m2,sd1,sd2,cov12=0) {
	abs(m1/m2) * sqrt( (sd1/m1)^2 + (sd2/m2)^2 - 2*cov12/(m1*m2) )
}

#error propagation for logarithms
log.sd <- function(m1,sd1,base=10) {
	abs(sd1/(m1*log(base)))
}

floor0 <- function(xs,bottom=0) sapply(xs, function(x) if (x < bottom) bottom else x)

#Baldi & Long's formula
bnl <- function(pseudo.n,n,model.sd,empiric.sd) {
	sqrt((pseudo.n * model.sd^2 + (n - 1) * empiric.sd^2)/(pseudo.n + n - 2))
}

#' Helper function to form mean, stdev and average raw count for a given condition
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

#' Helper function to apply raw count / frequecy filtering
#' @param msc dataframe containing the output of the 'mean.sd.count()' function.
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



##########################
# SCRAPS BELOW!!!!!!!!!!
##########################

op <- par(mfrow=c(2,2))
with(msc,plot(nonselect.count,phi.sd,log="xy",pch=".",main="Unfiltered"))
with(msc[is.na(msc$filter),],plot(nonselect.count,phi.sd,log="xy",pch=".",main="Filtered"))
with(msc,plot(logPhi,logPhi.sd,log="y",pch=".",main="Unfiltered"))
with(msc[is.na(msc$filter),],plot(logPhi,logPhi.sd,log="y",pch=".",main="Filtered"))
par(op)



layout(cbind(1,2,3))
plot(sd.poisson,sd.empiric,log="xy",xlim=c(1e-8,1e-3),ylim=c(1e-8,1e-3),pch=".")
abline(0,1,col="gray",lty="dotted")
plot(sd.poisson,sd.bayes,log="xy",xlim=c(1e-8,1e-3),ylim=c(1e-8,1e-3),pch=".")
abline(0,1,col="gray",lty="dotted")
plot(sd.empiric,sd.bayes,log="xy",xlim=c(1e-8,1e-3),ylim=c(1e-8,1e-3),pch=".")
abline(0,1,col="gray",lty="dotted")






logcv <- log10(mscFiltered[,paste0(cond,".cv")])
# logmean <- log10(mscFiltered[,paste0(cond,".mean")])
logexpect <- log10(1/sqrt(mscFiltered[,paste0(cond,".count")]))
inverse <- 1/logexpect
z <- coefficients(lm(logcv~logexpect+inverse))
priorCV <- function(count) 10^(z[[1]] + z[["logexpect"]]*log10(1/sqrt(count)) + z[["inverse"]]/log10(1/sqrt(count)))
oneoversqr <- function(x)1/sqrt(x)
with(mscFiltered,plot(nonselect.count,nonselect.cv,log="xy"))
# curve(priorCV,from=.1,to=5000,col="red",add=TRUE)
curve(oneoversqr,from=.1,to=5000,col="gray",add=TRUE)
# runningMean <- with(mscFiltered,runningFunction(nonselect.count,nonselect.cv,10,100,fun=median,logScale=TRUE))
# lines(runningMean,col="red",lwd=2)

# m2cv <- function(m,n=600000) 1/sqrt(m*n)
# with(mscFiltered,plot(nonselect.mean,nonselect.cv,pch=".",log="xy"))
# curve(priorCV,add=TRUE,from=1e-6,to=1e-2)
# curve(m2cv,add=TRUE,from=1e-6,to=1e-2,col="red")
# x11()
# curve(priorCV,from=1e-6,to=1e-2)
# curve(m2cv,add=TRUE,from=1e-6,to=1e-2,col="red")
# x11()
with(mscFiltered,plot(1/sqrt(nonselect.count),nonselect.cv,log="xy",xlim=c(1e-4,1),ylim=c(1e-4,1)))
runningMean <- with(mscFiltered,runningFunction(1/sqrt(nonselect.count),nonselect.cv,0.02,100))
lines(runningMean,col="red",lwd=2)
abline(0,1,col="gray",lty="dotted")

abline(v=1/sqrt(600))

