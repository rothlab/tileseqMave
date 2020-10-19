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
	# latestCount <- latestSubDir(parentDir=dataDir,pattern="_mut_call$|mut_count$")
	latestScore <- latestSubDir(parentDir=dataDir,pattern="_scores$")
	# if (latestCount[["timeStamp"]] != latestScore[["timeStamp"]]) {
	# 	stop("latest score folder does not match latest counts folder time stamp!")
	# }

	#create a matching output directory
	outDir <- paste0(dataDir,latestScore[["label"]],latestScore[["timeStamp"]],"_fscores/")
	if (!dir.exists(outDir)) {
		dir.create(outDir,recursive=TRUE,showWarnings=FALSE)
	}

	# logInfo("Reading count data")
	# marginalCountFile <- paste0(latestCount[["dir"]],"/marginalCounts.csv")
	# marginalCounts <- read.csv(marginalCountFile)
	# rownames(marginalCounts) <- marginalCounts$hgvsc

	# #filter out frameshifts and indels
	# toAA <- extract.groups(marginalCounts$aaChange,"\\d+(.*)$")
	# indelIdx <- which(toAA=="-" | nchar(toAA) > 1)
	# marginalCounts <- marginalCounts[-indelIdx,]

	# #annotate with positions, tiles and regions
	# marginalCounts$position <- as.integer(extract.groups(marginalCounts$codonChange,"(\\d+)"))
	# marginalCounts$region <- params$pos2reg(marginalCounts$position)
	# marginalCounts$tile <- params$pos2tile(marginalCounts$position)

	
	#Calculate means, stdevs and average count for each condition
	names(conditionOrder) <- paste0("bin",seq_along(conditionOrder))
	conditions <- c(conditionOrder,getNonselectFor(conditionOrder[[1]],params),getWTControlFor(conditionOrder[[1]],params))
	names(conditions) <- c(paste0("bin",seq_along(conditionOrder)),"all","wt")

	# msc <- do.call(cbind,lapply(conditions, mean.sd.count, marginalCounts, tp="1", params))





	# foo <- msc[,paste0(names(conditionOrder),".mean")]

	# subsample <- foo[sample(nrow(foo),100),]
	# plot(NA,type="n",xlim=c(1,6),ylim=c(0,max(as.matrix(subsample))),axes=FALSE,xlab="bin",ylab="freq")
	# axis(1,seq_along(conditionOrder),conditionOrder)
	# axis(2)
	# apply(subsample,1,function(x)lines(seq_along(x),x))


	tp <- "1"
	names(conditionOrder) <- paste0("bin",seq_along(conditionOrder))

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

	allScores <- do.call(cbind,setNames(lapply(scoreTables,
		function(x) x[,c("logPhi","logPhi.sd")]
	), names(conditionOrder)))
	allScores <- cbind(scoreTables[[1]][,1:5],allScores)

	#only variants with at least one value that wasn't filtered out
	allScores <- allScores[apply(allScores[,-(1:5)],1,function(x)!all(is.na(x))),]

	`%~%` <- function(f,g) function(...) f(g(...))
	library(mclust)

	sds <- allScores[,paste0(names(conditionOrder),".logPhi.sd")]
	qualFilter <- apply(sds,1,function(x)all(x < 1))

	lphis <- allScores[qualFilter,paste0(names(conditionOrder),".logPhi")]
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
	invisible(lapply(names(clsizes), function(cli){
		cli <- as.integer(cli)
		subsample <- finlphis[mcfinlphis$classification==cli,][sample(clsizes[[cli]],50),]
		plot(NA,type="n",xlim=c(1,5),ylim=c(-2.5,1.5),axes=FALSE,
			xlab="bin",ylab="logPhi",
			main=sprintf("Cluster #%d (%d variants)",cli,clsizes[[cli]])
		)
		axis(1,seq_along(conditionOrder),conditionOrder)
		axis(2)
		lapply(1:nrow(subsample),function(i) lines(1:ncol(subsample),subsample[i,],col=plotcols[cli]))
	}))
	par(opar)
	dev.off()


	model <- function(x,theta) {
		theta[[1]]+ x*theta[[2]]
	}
	xs <- seq(-2,2,1)

	fits <- as.df(pbmclapply(1:nrow(allScores),function(row) {
		values <- unlist(allScores[row,paste0(names(conditionOrder),".logPhi")])
		sds <- unlist(allScores[row,paste0(names(conditionOrder),".logPhi.sd")])
		filter <- which(!is.na(values) & !is.na(sds))
		if (length(filter) < 3) {
			return(c(intercept=NA,slope=NA,logLikeli=NA))
		}
		objective <- function(theta) {
			preds <- model(xs,theta)
			logLikeli <- sum(dnorm(preds[filter],values[filter],sds[filter],log=TRUE))
			return(logLikeli)
		}
		theta.start <- c(intercept=0,slope=0)
		z <- optim_nm(objective,start=theta.start,maximum=TRUE)
		theta.optim <- setNames(z$par,names(theta.start))
		c(theta.optim,logLikeli=z$function_value)
	},mc.cores=6))
	fits$logLikeli <- fits$logLikeli-max(fits$logLikeli,na.rm=TRUE)

	allScoresFits <- cbind(allScores,fits)

	pdf("modelFits.pdf",5,5)
	hist(fits$logLikeli,breaks=50,
		col="firebrick3",border=NA,
		xlab="Linear model log likelihood", 
		main="Model qualities"
	)
	dev.off()


	drawFit <- function(row,col="royalblue3") {
		values <- unlist(allScoresFits[
			row, paste0(names(conditionOrder),".logPhi")
		])
		sds <- unlist(allScoresFits[
			row, paste0(names(conditionOrder),".logPhi.sd")
		])
		filter <- which(!is.na(values) & !is.na(sds))
		plot(xs,values,
			ylim=c(-3,2),pch=20,
			main=allScoresFits[row,"hgvsp"],
			axes=FALSE,xlab="bins",ylab="log(phi)"
		)
		axis(1,at=xs,conditionOrder)
		axis(2)
		yogitools::errorBars(xs,values,sds)
		abline(h=0,col="gray",lty="dotted")
		if (all(is.finite(unlist(allScoresFits[row,c("intercept","slope")])))) {
			abline(allScoresFits[row,"intercept"],
				allScoresFits[row,"slope"],col=col,lwd=2
			)
		}
		mtext(sprintf("logL=%.01f",allScoresFits[row,"logLikeli"]),line=0,cex=0.7)
	}

	pdf("sample_fits.pdf",10,10)
	layout(matrix(1:25,nrow=5,ncol=5))
	for (i in sample(nrow(allScoresFits),25)) {
		drawFit(i)
	}
	dev.off()

	# op <- par(mfrow=c(1,2))
	# with(allScoresFits,plot(slope,logLikeli,pch=20,col=colAlpha(1,0.2)))
	# with(allScoresFits,plot(intercept,slope,pch=20,col=colAlpha(1,0.2)))
	# par(op)

	panel <- function(x,y,...) {
		points(x,y,...)
		abline(h=0,v=0,col="gray",lty="dashed")
	}
	pdf("modelParameters.pdf",5,5)
	pairs(fits,pch=".",col=colAlpha(1,0.2),lower.panel=panel,upper.panel=panel)
	dev.off()



	# #try again without lowest and highest
	# finlphis <- finlphis[,-c(1,5)]
	# mcfinlphis <- Mclust(finlphis)

	# labels <- conditionOrder[-c(1,5)]
	# plotcols <- sapply(rainbow(mcfinlphis$G+1),yogitools::colAlpha,0.2)

	# pdf("pairsClustering.pdf",7,7)
	# pairs(finlphis,col=plotcols[mcfinlphis$classification],labels=labels,pch=20)
	# dev.off()


	# clsizes <- table(mcfinlphis$classification)


	# pdf("clusterTrails.pdf",7,7)
	# # opar <- par(mfrow=c(3,3),mar=c(5,4,1,0),bg="gray20",fg="gray90",col.lab="gray90",col.main="gray90",col.axis="gray90")
	# opar <- par(mfrow=c(3,3),mar=c(5,4,1,0),oma=c(1,1,1,1))
	# lapply(names(clsizes), function(cli){
	# 	cli <- as.integer(cli)
	# 	subsample <- finlphis[mcfinlphis$classification==cli,][sample(clsizes[[cli]],50),]
	# 	plot(NA,type="n",xlim=c(1,3),ylim=c(-2.5,1.5),axes=FALSE,
	# 		xlab="bin",ylab="logPhi",
	# 		main=sprintf("Cluster #%d (%d variants)",cli,clsizes[[cli]])
	# 	)
	# 	axis(1,seq_along(labels),labels)
	# 	axis(2)
	# 	invisible(lapply(1:nrow(subsample),function(i) lines(1:ncol(subsample),subsample[i,],col=plotcols[cli])))
	# })
	# par(opar)
	# dev.off()




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



	# pathoIdx <- do.call(c,lapply(prs$hgvs,function(h) which(allScoresFits$hgvsp==h)))
	# beniIdx <- do.call(c,lapply(nrs$hgvs,function(h) which(allScoresFits$hgvsp==h)))
	pathoFilter <- allScoresFits$hgvsp %in% prs$hgvs
	beniFilter <- allScoresFits$hgvsp %in% nrs$hgvs

	pdf("prsSlopes.pdf",14,14)
	op <- par(mfrow=c(6,6))
	for (i in sample(which(pathoFilter & allScoresFits$logLikeli > -5),36)) drawFit(i,col="firebrick3")
	par(op)
	dev.off()

	pdf("nrsSlopes.pdf",6,6)
	op <- par(mfrow=c(2,2))
	for (i in which(beniFilter & allScoresFits$logLikeli > -5)) drawFit(i,col="chartreuse3")
	par(op)
	dev.off()

	library(beeswarm)

	pdf("prsNrsSlopes.pdf",5,5)
	prsNrsSlopes <- list(
		PRS=allScoresFits[pathoFilter & allScoresFits$logLikeli > -5,"slope"],
		NRS=allScoresFits[beniFilter & allScoresFits$logLikeli > -5,"slope"]
	)
	beeswarm(prsNrsSlopes,pch=20,
		ylab="model slope",
		col=c("firebrick3","chartreuse3")
	)
	bxplot(prsNrsSlopes,add=TRUE)
	dev.off()


	######
	# Export slope map
	neutral <- 0.1
	damaging <- -0.4
	fitness <- (allScoresFits$slope-damaging)/(neutral-damaging)
	error <- -allScoresFits$logLikeli/6
	mapOut <- data.frame(
		hgvs_pro = allScoresFits$hgvsp,
		score=fitness,
		se=error
	)
	write.csv(mapOut,"slopeMap.csv",row.names=FALSE)


	#Compare against PP2, Provean and SIFT
	builder <- new.hgvs.builder.p(3)

	provean <- read.delim("../../LDLR_provean.tsv")
	provean$hgvs <- sapply(1:nrow(provean), function(i) with(provean[i,],{
		if (RESIDUE_REF==RESIDUE_ALT) {
			builder$synonymous(POSITION,RESIDUE_REF)
		} else {
			builder$substitution(POSITION,RESIDUE_REF,RESIDUE_ALT)
		}
	}))
	rownames(provean) <- provean$hgvs
	pp2 <- read.delim("../../LDLR_pph2-fixed.txt")
	pp2$hgvs <- sapply(1:nrow(pp2), function(i) with(pp2[i,],{
		if (aa1==aa2) {
			builder$synonymous(o_pos,aa1)
		} else {
			builder$substitution(o_pos,aa1,aa2)
		}
	}))
	rownames(pp2) <- pp2$hgvs

	allScores$pp2 <- pp2[allScores$hgvsp,"pph2_prob"]
	allScores$provean <- provean[allScores$hgvsp,"SCORE"]
	allScores$sift <- provean[allScores$hgvsp,"SCORE.1"]

	plotAll <- function(sdCutoff=0.3,against="provean",label="PROVEAN") {
		invisible(lapply(1:5,function(bini) {
			cName <- conditionOrder[[bini]]
			mName <- sprintf("bin%d.logPhi",bini)
			sName <- sprintf("bin%d.logPhi.sd",bini)
			filtered <- allScores[which(allScores[,sName] < sdCutoff),]
			plot(filtered[,c(mName,against)],
				xlab=paste(cName,"log(phi)"),ylab=label,
				pch=".",col=colAlpha("black",0.3)
			)
			sc <- cor(fin(filtered[,c(mName,against)]),method="spearman")[1,2]
			mtext(sprintf("Spearman's rho = %.02f",sc))
		}))
	}

	pdf("predictor_comparison.pdf",11,8.5)
	op <- par(mfrow=c(3,5),oma=c(3,2,3,2))
	plotAll(against="provean",label="PROVEAN")
	plotAll(against="pp2",label="PolyPhen-2")
	plotAll(against="sift",label="SIFT")
	mtext("Filtered for replicate agreement",outer=TRUE)
	par(op)
	invisible(dev.off())

	pdf("predictor_comparison_unfiltered.pdf",11,8.5)
	op <- par(mfrow=c(3,5),oma=c(3,2,3,2))
	plotAll(sdCutoff=Inf,against="provean",label="PROVEAN")
	plotAll(sdCutoff=Inf,against="pp2",label="PolyPhen-2")
	plotAll(sdCutoff=Inf,against="sift",label="SIFT")
	mtext("Unfiltered",outer=TRUE)
	par(op)
	invisible(dev.off())

}