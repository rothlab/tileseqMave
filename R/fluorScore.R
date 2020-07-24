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

#' score fluorescence data
#' 
#' @param dataDir working data directory
#' @param conditionOrder vector containting the order in which the conditions represent fluorescence bins
#' @param paramFile input parameter file. defaults to <dataDir>/parameters.json
#' @return NULL. Results are written to file.
#' @export
fluorScore <- function(dataDir,conditionOrder,paramFile=paste0(dataDir,"parameters.json"),srOverride=FALSE) {

	# conditionOrder <- c("lowest","lowmid","mid","midhigh","highest")

	op <- options(stringsAsFactors=FALSE)

	library(yogitools)
	library(hgvsParseR)
	library(pbmcapply)
	library(optimization)

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
	params <- parseParameters(paramFile,srOverride=srOverride)

	if (!all(conditionOrder %in% getSelects(params))) {
		stop("conditionOrder argument contains unknown conditions!")
	}


	#find counts and scores folder
	latestCount <- latestSubDir(parentDir=dataDir,pattern="_mut_call$|mut_count$")
	latestScore <- latestSubDir(parentDir=dataDir,pattern="_scores$")
	if (latestCount[["timeStamp"]] != latestScore[["timeStamp"]]) {
		stop("latest score folder does not match latest counts folder time stamp!")
	}

	#create a matching output directory
	outDir <- paste0(dataDir,latestCount[["label"]],latestCount[["timeStamp"]],"_fscores/")
	if (!dir.exists(outDir)) {
		dir.create(outDir,recursive=TRUE,showWarnings=FALSE)
	}

	logInfo("Reading count data")
	marginalCountFile <- paste0(latestCount[["dir"]],"/marginalCounts.csv")
	marginalCounts <- read.csv(marginalCountFile)
	rownames(marginalCounts) <- marginalCounts$hgvsc

	#filter out frameshifts and indels
	toAA <- extract.groups(marginalCounts$aaChange,"\\d+(.*)$")
	indelIdx <- which(toAA=="-" | nchar(toAA) > 1)
	marginalCounts <- marginalCounts[-indelIdx,]

	#annotate with positions, tiles and regions
	marginalCounts$position <- as.integer(extract.groups(marginalCounts$codonChange,"(\\d+)"))
	marginalCounts$region <- params$pos2reg(marginalCounts$position)
	marginalCounts$tile <- params$pos2tile(marginalCounts$position)

	
	#Calculate means, stdevs and average count for each condition
	names(conditionOrder) <- paste0("bin",seq_along(conditionOrder))
	conditions <- c(conditionOrder,getNonselectFor(conditionOrder[[1]],params),getWTControlFor(conditionOrder[[1]],params))
	names(conditions) <- c(paste0("bin",seq_along(conditionOrder)),"all","wt")

	msc <- do.call(cbind,lapply(conditions, mean.sd.count, marginalCounts, tp="1", params))





	foo <- msc[,paste0(names(conditionOrder),".mean")]

	subsample <- foo[sample(nrow(foo),100),]
	plot(NA,type="n",xlim=c(1,6),ylim=c(0,max(as.matrix(subsample))),axes=FALSE,xlab="bin",ylab="freq")
	axis(1,seq_along(conditionOrder),conditionOrder)
	axis(2)
	apply(subsample,1,function(x)lines(seq_along(x),x))


	tp <- "1"

	scoreTables <- lapply(names(conditionOrder), function(bini) {
		sCond <- conditionOrder[[bini]]
		scoreFile <- paste0(latestScore[["dir"]],"/",sCond,"_t",tp,"_complete.csv")
		if (!file.exists(scoreFile)) {
			logWarn("No scores found for",sCond)
			return(NA)
		}
		scores <- read.csv(scoreFile)
		#delete filtered values
		scores[!is.na(scores$filter),c("logPhi","logPhi.sd")] <- NA
		#extract relevant columns
		return(scores[,c(1:6,which(colnames(scores) %in% c("logPhi","logPhi.sd")))])
	})

	mutIDs <- lapply(scoreTables,function(x)x[,1])
	if (!all(sapply(mutIDs,function(x)identical(x,mutIDs[[1]])))) {
		stop("Score tables do not contain exactly the same mutants!")
	}

	allScores <- do.call(cbind,setNames(lapply(scoreTables,function(x) x[,c("logPhi","logPhi.sd")]),names(conditionOrder)))
	allScores <- cbind(scoreTables[[1]][,1:5],allScores)

	#only variants with at least one value that wasn't filtered out
	allScores <- allScores[apply(allScores[,-(1:5)],1,function(x)!all(is.na(x))),]

	`%~%` <- function(f,g) function(...) f(g(...))
	library(mclust)

	lphis <- allScores[,paste0(names(conditionOrder),".logPhi")]
	finlphis <- fin(lphis)
	mcfinlphis <- Mclust(finlphis)

	plotcols <- sapply(rainbow(mcfinlphis$G+1),yogitools::colAlpha,0.2)
	pdf("pairsClustering.pdf",7,7)
	pairs(finlphis,col=plotcols[mcfinlphis$classification],labels=conditionOrder,pch=20)
	dev.off()


	clsizes <- table(mcfinlphis$classification)


	pdf("clusterTrails.pdf",7,7)
	# opar <- par(mfrow=c(3,3),mar=c(5,4,1,0),bg="gray20",fg="gray90",col.lab="gray90",col.main="gray90",col.axis="gray90")
	opar <- par(mfrow=c(3,3),mar=c(5,4,1,0))
	lapply(names(clsizes), function(cli){
		cli <- as.integer(cli)
		subsample <- finlphis[mcfinlphis$classification==cli,][sample(clsizes[[cli]],50),]
		plot(NA,type="n",xlim=c(1,5),ylim=c(-2.5,1.5),axes=FALSE,
			xlab="bin",ylab="logPhi",
			main=sprintf("Cluster #%d (%d variants)",cli,clsizes[[cli]])
		)
		axis(1,seq_along(conditionOrder),conditionOrder)
		axis(2)
		invisible(lapply(1:nrow(subsample),function(i) lines(1:ncol(subsample),subsample[i,],col=plotcols[cli])))
	})
	par(opar)
	dev.off()


	#try again without lowest and highest
	finlphis <- finlphis[,-c(1,5)]
	mcfinlphis <- Mclust(finlphis)

	labels <- conditionOrder[-c(1,5)]
	plotcols <- sapply(rainbow(mcfinlphis$G+1),yogitools::colAlpha,0.2)

	pdf("pairsClustering.pdf",7,7)
	pairs(finlphis,col=plotcols[mcfinlphis$classification],labels=labels,pch=20)
	dev.off()


	clsizes <- table(mcfinlphis$classification)


	pdf("clusterTrails.pdf",7,7)
	# opar <- par(mfrow=c(3,3),mar=c(5,4,1,0),bg="gray20",fg="gray90",col.lab="gray90",col.main="gray90",col.axis="gray90")
	opar <- par(mfrow=c(3,3),mar=c(5,4,1,0),oma=c(1,1,1,1))
	lapply(names(clsizes), function(cli){
		cli <- as.integer(cli)
		subsample <- finlphis[mcfinlphis$classification==cli,][sample(clsizes[[cli]],50),]
		plot(NA,type="n",xlim=c(1,3),ylim=c(-2.5,1.5),axes=FALSE,
			xlab="bin",ylab="logPhi",
			main=sprintf("Cluster #%d (%d variants)",cli,clsizes[[cli]])
		)
		axis(1,seq_along(labels),labels)
		axis(2)
		invisible(lapply(1:nrow(subsample),function(i) lines(1:ncol(subsample),subsample[i,],col=plotcols[cli])))
	})
	par(opar)
	dev.off()




	prs <- read.csv("prs_patho.csv")
	nrs <- read.csv("nrs_beni.csv")

	pathoIdx <- do.call(c,lapply(prs$hgvs,function(h) which(allScores$hgvsp==h)))
	beniIdx <- do.call(c,lapply(nrs$hgvs,function(h) which(allScores$hgvsp==h)))

	opar <- par(mfrow=c(1,2))
	pathotracks <- lphis[pathoIdx,]
	plot(NA,type="n",xlim=c(1,5),ylim=c(-2.5,1.5),axes=FALSE,xlab="bin",ylab="logPhi",main="PRS (P/LP)")
	axis(1,seq_along(conditionOrder),conditionOrder)
	axis(2)
	invisible(lapply(1:nrow(pathotracks),function(i) lines(1:ncol(pathotracks),pathotracks[i,],col="firebrick3")))
	benitracks <- lphis[beniIdx,]
	plot(NA,type="n",xlim=c(1,5),ylim=c(-2.5,1.5),axes=FALSE,xlab="bin",ylab="logPhi",main="NRS (B/LB)")
	axis(1,seq_along(conditionOrder),conditionOrder)
	axis(2)
	invisible(lapply(1:nrow(benitracks),function(i) lines(1:ncol(benitracks),benitracks[i,],col="chartreuse3")))
	par(opar)

}