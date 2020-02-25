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
selectionQC <- function(dataDir,paramFile=paste0(dataDir,"parameters.json"),logger=NULL,mc.cores=6) {


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


	#find scores folder
	subDirs <- list.dirs(dataDir,recursive=FALSE)
	scoreDirs <- subDirs[grepl("_scores$",subDirs)]
	if (length(scoreDirs) == 0) {
		stop("No mutation call output found!")
	}
	latestScoreDir <- sort(scoreDirs,decreasing=TRUE)[[1]]
	#find corresponding count directory
	latestCountDir <- sub("_scores$","_mut_call",latestScoreDir)
	#extract time stamp
	timeStamp <- extract.groups(latestScoreDir,"/(\\d{4}-\\d{2}-\\d{2}-\\d{2}-\\d{2}-\\d{2})")

	#create a matching output directory
	outDir <- paste0(dataDir,timeStamp,"_QC/")
	if (!dir.exists(outDir)) {
		dir.create(outDir,recursive=TRUE,showWarnings=FALSE)
	}

	logInfo("Reading count data")
	marginalCountFile <- paste0(latestCountDir,"/marginalCounts.csv")
	marginalCounts <- read.csv(marginalCountFile)

	#filter out frameshifts and indels
	toAA <- extract.groups(marginalCounts$aaChange,"\\d+(.*)$")
	indelIdx <- which(toAA=="-" | nchar(toAA) > 1)
	marginalCounts <- marginalCounts[-indelIdx,]


	#iterate over conditions
	for (sCond in getSelects(params)) {

		#iterate over timepoints
		for (tp in params$timepoints$`Time point name`) {

			logInfo("Processing condition",sCond, "; time",tp)

			#load score table for this condition
			scoreFile <- paste0(latestScoreDir,"/",sCond,"_t",tp,"_complete.csv")
			if (!file.exists(scoreFile)) {
				logWarn("No score file found! Skipping...")
				next
			}
			scores <- read.csv(scoreFile)

			#run replicate correlation analysis
			if (params$numReplicates > 1) {
				replicateCorrelation(scores, marginalCounts, params, sCond, tp, outDir)
			}

			#TODO: Regularization analysis

			#TODO: Score distributions & syn/non medians

			#TODO: Error profile

		}
	}

}


replicateCorrelation <- function(scores, marginalCounts, params, sCond, tp, outDir) {

	nrep <- params$numReplicates

	#pull up matching nonselect and WT controls
	nCond <- getNonselectFor(sCond,params)
	condQuad <- c(
		select=sCond,
		nonselect=nCond,
		selWT=getWTControlFor(sCond,params),
		nonWT=getWTControlFor(nCond,params)
	)

	#Workaround: Since R doesn't like variable names starting with numerals, 
	# we need to adjust any of those cases
	if (any(grepl("^\\d",condQuad))) {
		culprits <- which(grepl("^\\d",condQuad))
		#we add the prefix "X" to the name, just like the table parser does.
		condQuad[culprits] <- paste0("X",condQuad[culprits])
	}

	#check that labels match between tables
	if (!all(scores$hgvsp == marginalCounts$hgvsp)) {
		stop("Variants differ between counts and scores!")
	}

	#replicate correlations
	repMatrix <- do.call(cbind,lapply(1:nrep,
		function(repi) sprintf("%s.t%s.rep%d.frequency",condQuad,tp,repi)
	))
	rownames(repMatrix) <- names(condQuad)

	#extract replicate values for this condition
	repValues <- lapply(1:nrep, function(repi) {
		selFreq <-  floor0(
			marginalCounts[,repMatrix["select",repi]] - 
			marginalCounts[,repMatrix["selWT",repi]]
		)
		nonFreq <-  floor0(
			marginalCounts[,repMatrix["nonselect",repi]] - 
			marginalCounts[,repMatrix["nonWT",repi]]
		)
		logphi <- log10(selFreq / nonFreq)
		data.frame(select=selFreq,nonselect=nonFreq,logphi=logphi)
	})
	names(repValues) <- as.character(1:nrep)
	repValues <- do.call(cbind,repValues)
	#apply filter from scoring function
	repValues <- repValues[is.na(scores$filter),]


	if (nrep == 2) {
		outfile <- paste0(outDir,sCond,"_t",tp,"_replicates.pdf")
		pdf(outfile,10,5)
		layout(cbind(1,2))
		plot(
			repValues[,"1.nonselect"],repValues[,"2.nonselect"],
			xlab="Frequency Replicate 1", ylab="Frequency Replicate 2",
			main=sprintf(
				"non-select frequencies R = %.03f",
				cor(fin(log10(repValues[,sprintf("%d.nonselect",1:2)])))[1,2]
			),
			pch=".",
			log="xy"
		)
		plot(
			repValues[,"1.logphi"],repValues[,"2.logphi"],
			xlab="log(phi) Replicate 1", ylab="log(phi) Replicate 2",
			main=sprintf(
				"select / nonselect log-ratio R = %.03f",
				cor(fin(repValues[,sprintf("%d.logphi",1:2)]))[1,2]
			),pch="."
		)
		invisible(dev.off())
	} else {
		panel.cor <- function(x, y,...){
			usr <- par("usr"); on.exit(par(usr)); par(usr=c(0,1,0,1))
			r <- cor(fin(cbind(x,y)))[1,2]
			txt <- sprintf("R = %.02f",r)
			# cex.cor <- 0.8/strwidth(txt)
			text(0.5, 0.5, txt)
		}
		labels <- paste0("rep.",1:nrep)	
		imgSize <- max(4,nrep)

		outfile <- paste0(outDir,sCond,"_t",tp,"_ns_replicates.pdf")
		pdf(outfile,imgSize,imgSize)
		pairs(
			repValues[,sprintf("%d.nonselect",1:nrep)],
			lower.panel=panel.cor,pch=".",labels=labels,
			main="non-select frequencies"
		)
		invisible(dev.off())

		outfile <- paste0(outDir,sCond,"_t",tp,"_phi_replicates.pdf")
		pdf(outfile,imgSize,imgSize)
		pairs(
			repValues[,sprintf("%d.nonselect",1:nrep)],
			lower.panel=panel.cor,pch=".",labels=labels,
			main="select / nonselect log-ratios"
		)
		invisible(dev.off())
	}
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

