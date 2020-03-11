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
selectionQC <- function(dataDir,paramFile=paste0(dataDir,"parameters.json"),logger=NULL, sdCutoff=0.3,srOverride=TRUE) {


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
			if (params$numReplicates[[sCond]] > 1) {
				replicateCorrelation(scores, marginalCounts, params, sCond, tp, outDir)
			}

			#load error model file
			modelFile <- paste0(latestScoreDir,"/",sCond,"_t",tp,"_errorModel.csv")
			modelParams <- read.csv(modelFile)

			#Regularization analysis
			regularizationQC(scores,modelParams,params,sCond,tp,outDir)
			
			#Score distributions & syn/non medians
			scoreDistributions(scores,sCond,tp,outDir,sdCutoff)

			#Error profile
			errorProfile(scores,sCond,tp,outDir)

		}
	}

}


#' Draw replicate correlation plots
#' 
#' @param scores the score table
#' @param marginalCounts the marginal counts table
#' @param params the global parameter object
#' @param sCond the condition for which to draw the plot
#' @param tp the time point
#' @param outDir the output directory
#' @return NULL
replicateCorrelation <- function(scores, marginalCounts, params, sCond, tp, outDir) {

	# nrep <- params$numReplicates[[sCond]]

	#pull up matching nonselect and WT controls
	nCond <- getNonselectFor(sCond,params)
	sRep <- params$numReplicates[[sCond]]
	nsRep <- params$numReplicates[[msCond]]
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

	#replicate column name matrix
	repMatrix <- do.call(cbind,lapply(1:sRep,
		function(repi) sprintf("%s.t%s.rep%d.frequency",condQuad,tp,repi)
	))
	rownames(repMatrix) <- names(condQuad)

	#extract replicate values for this condition
	repValues <- lapply(1:sRep, function(repi) {
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
	names(repValues) <- as.character(1:sRep)
	repValues <- do.call(cbind,repValues)
	#apply filter from scoring function
	repValues <- repValues[is.na(scores$filter),]


	if (sRep == 2) {
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
		labels <- paste0("rep.",1:sRep)	
		imgSize <- max(4,sRep)

		outfile <- paste0(outDir,sCond,"_t",tp,"_ns_replicates.pdf")
		pdf(outfile,imgSize,imgSize)
		pairs(
			repValues[,sprintf("%d.nonselect",1:sRep)],
			lower.panel=panel.cor,pch=".",labels=labels,
			main="non-select frequencies"
		)
		invisible(dev.off())

		outfile <- paste0(outDir,sCond,"_t",tp,"_phi_replicates.pdf")
		pdf(outfile,imgSize,imgSize)
		pairs(
			repValues[,sprintf("%d.nonselect",1:sRep)],
			lower.panel=panel.cor,pch=".",labels=labels,
			main="select / nonselect log-ratios"
		)
		invisible(dev.off())
	}

	return(NULL)
}

#' Draw error regulariztion model QC plots
#' 
#' @param scores the score table
#' @param params the global parameter object
#' @param sCond the condition for which to draw the plot
#' @param tp the time point
#' @param outDir the output directory
#' @return NULL
regularizationQC <- function(scores,modelParams,params,sCond,tp,outDir) {

	#calculate tile assignments
	tileStarts <- params$tiles[,"Start AA"]
	positions <- as.integer(extract.groups(scores$codonChange,"(\\d+)")[,1])
	tiles <- sapply(positions,function(pos) max(which(tileStarts <= pos)))
	
	outfile <- paste0(outDir,sCond,"_t",tp,"_errorModel.pdf")
	pdf(outfile,8.5,11)
	opar <- par(mfrow=c(3,2))
	for (tile in params$tiles[,"Tile Number"]) {
		# model.fit <- tryCatch({
		# 	fit.cv.model(scores[which(tiles==tile),])
		# },error=function(e) {
		# 	NULL
		# })
		theta <- modelParams[tile,paste0("nonselect.",c("static","additive","multiplicative"))]
		cv.model <- function(count) {
			10^sapply(log10(1/sqrt(count)),function(x) max(theta[[1]], theta[[2]] + theta[[3]]*x))
		}

		with(scores[which(tiles==tile),],{

			plot(nonselect.count,nonselect.cv,log="xy",main=paste("Tile",tile))
			runningMean <- runningFunction(
				nonselect.count,nonselect.cv,nbins=20,logScale=TRUE
			)
			nsSamples <- seq(1,max(nonselect.count,na.rm=TRUE),length.out=100)
			lines(nsSamples,1/sqrt(nsSamples),col="chartreuse3",lty="dashed",lwd=2)
			lines(runningMean,col="firebrick3",lwd=2)
			
			lines(nsSamples,cv.model(nsSamples),col="blue",lwd=2)
			mtext(sprintf(
				"stat.=%.02f; add.=%.02f; mult.=%.02f",
				theta[[1]],theta[[2]],theta[[3]]
			))

			nonselect.sd.poisson <- 1/sqrt(nonselect.count)*nonselect.mean
			nonselect.sd.poisson[which(nonselect.mean==0)] <- 0
			plot(nonselect.sd.poisson, nonselect.sd,log="xy")
			runningMean <- runningFunction(
				nonselect.sd.poisson, nonselect.sd,nbins=20,logScale=TRUE
			)
			lines(runningMean,col="firebrick3",lwd=2)
			abline(0,1,col="chartreuse3",lty="dashed",lwd=2)
			depth <- mean(nonselect.count/nonselect.mean,na.rm=TRUE)
			lines(
				sqrt(nsSamples)/depth,
				cv.model(nsSamples)*(nsSamples/depth),
				col="blue",lwd=2
			)
		})
	}
	par(opar)
	invisible(dev.off())
}


# #' Create a model of the error in for a subset of marginal frequencies
# #' 
# #' @param  subcores a subset of the count table (for a given tile)
# #' @return a list containing the model function (mapping counts to expected CV), and the model parameters
# fit.cv.model <- function(subscores) {

# 	# subscores <- scores[which(tiles==tile),]

# 	#poisson model of CV
# 	cv.poisson <- 1/sqrt(subscores$nonselect.count)
# 	#filter out NAs and infinites
# 	filter <- which(is.na(cv.poisson) | is.infinite(cv.poisson) | 
# 		is.na(subscores$nonselect.cv) | is.infinite(subscores$nonselect.cv))

# 	#model function
# 	log.cv.model <- function(log.cv.poisson,static,additive,multiplicative) {
# 		sapply(log.cv.poisson, function(x) max(static, multiplicative*x + additive))
# 	}
# 	#objective function
# 	objective <- function(theta) {
# 		static <- theta[[1]]
# 		additive <- theta[[2]]
# 		multiplicative <- theta[[3]]
# 		reference <- log10(subscores$nonselect.cv[-filter])
# 		prediction <- log.cv.model(log10(cv.poisson[-filter]),static,additive,multiplicative)
# 		rmsd <- sum((prediction-reference)^2,na.rm=TRUE)
# 		return(rmsd)
# 	}
# 	#run optimization
# 	theta.start <- c(static=-2,additive=0,multiplicative=1)
# 	z <- optim_nm(objective,start=theta.start)
# 	theta.optim <- setNames(z$par,names(theta.start))

# 	cv.model <- function(count) {
# 		cv.poisson <- 1/sqrt(count)
# 		10^log.cv.model(
# 			log10(cv.poisson),
# 			theta.optim[["static"]],
# 			theta.optim[["additive"]],
# 			theta.optim[["multiplicative"]]
# 		)
# 	}

# 	# plot(cv.poisson,subscores$nonselect.cv,log="xy")
# 	# lines(runningFunction(cv.poisson,subscores$nonselect.cv,nbins=20,logScale=TRUE),col=2)
# 	# abline(0,1)
# 	# lines(seq(0,1,0.01), 10^log.cv.model(log10(seq(0,1,0.01)),z$par[[1]],z$par[[2]],z$par[[3]]),col="blue")

# 	return(c(list(cv.model=cv.model),theta.optim))
# }

#' Draw score distribution plots
#' 
#' @param scores the score table
#' @param sCond the condition for which to draw the plot
#' @param tp the time point
#' @param outDir the output directory
#' @return NULL
scoreDistributions <- function(scores,sCond,tp,outDir,sdCutoff) {
	outfile <- paste0(outDir,sCond,"_t",tp,"_logPhiDistribution.pdf")
	pdf(outfile,11,8.5)
	layout(rbind(1,2,3,4),heights=c(1.2,1,1.2,1))
	drawDistributions(scores,Inf)
	drawDistributions(scores,sdCutoff)
	invisible(dev.off())
}

#' Delegation function to draw score distribution plots with a given filter setting
#' 
#' @param scores the score table
#' @param sdCutoff the stdev cutoff to apply
#' @return NULL
drawDistributions <- function(scores,sdCutoff) {
	#extract filtered scores
	synScores <- with(scores,logPhi[grepl("=$",hgvsp) & is.na(filter) & logPhi.sd < sdCutoff ])
	stopScores <- with(scores,logPhi[grepl("Ter$",hgvsp) & is.na(filter)& logPhi.sd < sdCutoff])
	misScores <- with(scores,logPhi[!grepl("Ter$|=$",hgvsp) & is.na(filter) & logPhi.sd < sdCutoff])
	allScores <- c(synScores,stopScores,misScores)

	#calculate plot ranges to nearest integers
	left <- floor(quantile(allScores,0.01,na.rm=TRUE))
	right <- ceiling(quantile(allScores,0.99,na.rm=TRUE))
	farleft <- floor(min(c(0,allScores),na.rm=TRUE))
	farright <- ceiling(max(c(1,allScores),na.rm=TRUE))
	#set histogram bins based on range (with extra space on the right)
	breaks <- seq(farleft,farright+0.1,0.1)
	#fill bins
	synHist <- hist(synScores,breaks=breaks,plot=FALSE)
	stopHist <- hist(stopScores,breaks=breaks,plot=FALSE)
	misHist <- hist(misScores,breaks=breaks,plot=FALSE)

	#draw top plot for syn/stop
	op <- par(mar=c(2,4,2,1)+.1)
	xs <- barplot(
		rbind(synHist$density,stopHist$density),
		beside=TRUE,col=c("darkolivegreen3","firebrick3"),
		border=NA,ylab="density",space=c(0,0),
		main=if (is.infinite(sdCutoff)) "Unfiltered" else bquote(sigma < .(sdCutoff))
	)
	grid(NA,NULL)
	#add x-axis
	pips <- left:right
	pipx <- colMeans(xs[,which(breaks %in% pips)])
	axis(1,at=pipx,labels=pips)
	#draw legend
	legend("left",
		c("nonsense","synonymous","missense"),
		fill=c("firebrick3","darkolivegreen3","gray"),
		border=NA,bty="n"
	)
	#draw bottom plot (for missense)
	par(mar=c(1,4,0,1)+.1)
	xs <- barplot(
		-misHist$density,
		border=NA,ylab="density",space=0,
	)
	grid(NA,NULL)
	#establish user coordinate function on barplot
	x0 <- xs[which(breaks==0),1]
	x1 <- xs[which(breaks==1),1]
	uCoord <- function(x)  x0 + x*(x1-x0)
	#draw medians
	abline(v=uCoord(median(synScores)),col="darkolivegreen4",lwd=2,lty="dotted")
	abline(v=uCoord(median(stopScores)),col="firebrick3",lwd=2,lty="dotted")
	text(
		uCoord(median(synScores)),
		-max(misHist$density)/3,
		sprintf("Synonymous median\n%.03f",median(synScores)),
		col="darkolivegreen4"
	)
	text(
		uCoord(median(stopScores)),
		-2*max(misHist$density)/3,
		sprintf("Nonsense median\n%.03f",median(stopScores)),
		col="firebrick3"
	)
	par(op)

}

#' Draw error profile for the dataset
#' 
#' @param scores the score table
#' @param sCond the condition for which to draw the plot
#' @param tp the time point
#' @param outDir the output directory
#' @return NULL
errorProfile <- function(scores,sCond,tp,outDir) {

	outfile <- paste0(outDir,sCond,"_t",tp,"_errorProfile.pdf")
	pdf(outfile,5,5)

	sdRange <- range(log10(scores[is.na(scores$filter),"score.sd"]))
	scoreRange <- range(scores[is.na(scores$filter),"score"])

	layout(rbind(c(2,4),c(1,3)),widths=c(0.8,0.2),heights=c(0.2,0.8))
	op <- par(mar=c(5,4,0,0)+.1) 
	with(scores[is.na(scores$filter),],plot(score,score.sd,log="y",pch=".",ylab=expression(sigma)))
	par(mar=c(0,4,1,0)+.1) 
	breaks <- seq(scoreRange[[1]],scoreRange[[2]],length.out=50)
	scoreHist <- with(scores[is.na(scores$filter),], hist(score,breaks=breaks,plot=FALSE))
	barplot(scoreHist$density,border=NA,ylab="density",space=0)
	par(mar=c(5,0,0,0)+.1) 
	breaks <- seq(sdRange[[1]],sdRange[[2]],length.out=50)
	sdHist <- with(scores[is.na(scores$filter),], hist(log10(score.sd),breaks=breaks,plot=FALSE))
	barplot(sdHist$density,border=NA,horiz=TRUE,space=0,xlab="density")
	par(op)

	invisible(dev.off())
	
}
