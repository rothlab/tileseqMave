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

#' QC plots for library coverage
#' @export
libraryQC <- function(dataDir,outdir,logger=NULL,mc.cores=8) {

	# library(hgvsParseR)
	# library(yogilog)
	library(yogitools)
	library(hash)
	library(pbmcapply)

	options(stringsAsFactors=FALSE)

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

	##############
	# Read and validate input data
	##############

	canRead <- function(filename) file.access(filename,mode=4) == 0
	
	#make sure data and out dir exist and ends with a "/"
	if (!grepl("/$",dataDir)) {
		dataDir <- paste0(dataDir,"/")
	}
	if (!grepl("/$",outdir)) {
		outdir <- paste0(outdir,"/")
	}
	if (!dir.exists(dataDir)) {
		#we don't use the logger here, assuming that whichever script wraps our function
		#catches the exception and writes to the logger (or uses the log error handler)
		stop("Data folder does not exist!")
	}
	if (!dir.exists(outdir)) {
		logWarn("Output folder does not exist, creating it now.")
		dir.create(outdir,recursive=TRUE)
	}

	mut2funcFile <- paste0(dataDir,"mut2func_info.csv")
	if (!canRead(mut2funcFile)) {
		stop("Unable to find or read 'mut2func_info.csv' file!")
	}
	seqDepthFile <- paste0(dataDir,"resultfile/sequencingDepth.csv")
	if (!canRead(seqDepthFile)) {
		stop("Unable to find or read 'sequencingDepth.csv' file!")
	}


	mut2func <- read.csv(mut2funcFile)
	tileRanges <- apply(extract.groups(colnames(mut2func)[-1],"X(\\d+)\\.(\\d+)"),2,as.integer)


	#Interpret sequencing depth file and obtain sample information
	seqDepthTable <- read.csv(seqDepthFile)
	rxGroups <- extract.groups(seqDepthTable$SampleID,"SampleID(\\d+)_(\\w+)(\\d+)")
	seqDepthTable$sample <- as.integer(rxGroups[,1])
	seqDepthTable$condition <- rxGroups[,2]
	seqDepthTable$replicate <- as.integer(rxGroups[,3])

	#Isolate first replicate of nonselect samples
	ns1.rows <- with(seqDepthTable,which(condition=="NS" & replicate == 1))
	ns1.depth <- seqDepthTable[ns1.rows,]

	logInfo("Compiling census...")

	#iterate over tiles/samples
	# censuses <- pbmclapply(1:nrow(ns1.depth), function(depth.row) {
	censuses <- lapply(1:nrow(ns1.depth), function(depth.row) {
		sample.id <- ns1.depth[depth.row,"sample"]
		sample.depth <- ns1.depth[depth.row,2]


		#read the file for mult-mutant CPMs (counts per million reads)
		multiCounts <- read.delim(
			paste0(dataDir,"mutationCallfile/",sample.id,"MultipleMut.txt"),
			header=FALSE
		)
		colnames(multiCounts) <- c("wtaa","pos","mutaa","wtcodon","mutcodon","cpm")

		wtaas <- strsplit(multiCounts$wtaa,"\\|")
		mutaas <- strsplit(multiCounts$mutaa,"\\|")
		positions <- strsplit(multiCounts$pos,"\\|")
		#create ids for multimutants
		multiVars <- mapply(function(w,p,m) paste0(w,p,m), wtaas, positions, mutaas)
		#number of aa-changes per multimutant
		nmuts <- sapply(multiVars,length)
		maxMuts <- max(nmuts)

		if (min(nmuts) < 2) {
			logWarn("MultiMuts contains single mutations!")
		}

		#create a census of the number of co-occurring mutations
		multiCensus <- tapply(multiCounts$cpm,nmuts,sum)
		#make sure no number is skipped!
		multiCensus <- multiCensus[as.character(2:maxMuts)]
		names(multiCensus) <- as.character(2:maxMuts)
		if (any(is.na(multiCensus))) {
			multiCensus[which(is.na(multiCensus))] <- 0
		}

		#infer expected single-mutant CPMS from multimutant CPMs
		multiTotals <- hash()
		#also infer complexity underpinning each aa-change (i.e. the number of unique multimutants that contain it)
		varComplexity <- hash()
		for (i in 1:length(multiVars)) {
			for (variant in multiVars[[i]]) {
				if (hash::has.key(variant,multiTotals)) {
					multiTotals[[variant]] <- multiTotals[[variant]] + multiCounts[i,"cpm"]
					varComplexity[[variant]] <- varComplexity[[variant]] + 1
				} else {
					multiTotals[[variant]] <- multiCounts[i,"cpm"]
					varComplexity[[variant]] <- 1
				}
			}
		}

		#Next, read the single-mutant file. Note that these are not true single-mutants, but rather
		#the total count regardless of in-cis context (i.e. mix of single and multimutants)
		singleCounts <- read.delim(
			paste0(dataDir,"mutationCallfile/",sample.id,"AAchange.txt"),
			header=FALSE
		)
		colnames(singleCounts) <- c("wtaa","pos","mutaa","wtcodon","mutcodon","cpm")

		#collapse by aa-change
		singleVariants <- with(singleCounts,mapply(function(w,p,m) paste0(w,p,m), wtaa, pos, mutaa))
		singleTotals <- hash()
		for (i in 1:length(singleVariants)) {
			v <- singleVariants[[i]]
			if (hash::has.key(v,singleTotals)) {
				singleTotals[[v]] <- singleTotals[[v]] + singleCounts[[i,"cpm"]]
			} else {
				singleTotals[[v]] <- singleCounts[[i,"cpm"]]
			}
		}

		#list of unique aa changes
		uniqueVars <- keys(singleTotals)

		#calculate the true single mutant CPMs by subtracting the corresponding multi-mutant CPMs
		singletonCounts <- sapply(uniqueVars, function(uvar) {
			if (has.key(uvar,multiTotals)) {
				singleTotals[[uvar]] - multiTotals[[uvar]]
			} else {
				singleTotals[[uvar]] 
			}
		})


		#Finally, read the deletion and inseration CPMs
		delCounts <- read.delim(
			paste0(dataDir,"mutationCallfile/",sample.id,"deletion.txt"),
			header=FALSE
		)
		insCounts <- read.delim(
			paste0(dataDir,"mutationCallfile/",sample.id,"insertion.txt"),
			header=FALSE
		)
		indelCPM <- sum(delCounts[,2])+sum(insCounts[,4])


		#Now, having compiled the required numbers, we can combine them into a census table.
		census <- c(
			indel = indelCPM,
			WT = 1e6 - (indelCPM + sum(singleCounts$cpm)),
			single = sum(singletonCounts),
			multiCensus
		) * 100 / 1e6


		#Fit a poisson-distribution to the census
		lambdas <- seq(0.01,4,0.01)
		dists <- sapply(lambdas, function(lambda) {
			sum(abs(census[-1] - dpois(0:maxMuts,lambda)*100))
		})
		# cat(min(dists),"\n")
		best.lambda <- lambdas[[which.min(dists)]]
		

		#Also generate a table that compares variant CPMs against variant complexity
		countVComplexity <- data.frame(
			cpm=singletonCounts[hash::keys(varComplexity)],
			complexity=hash::values(varComplexity)
		)

		return(list(census=census,lambda=best.lambda,complexity=countVComplexity))

	# },mc.cores=mc.cores)
	})


	censusRows <- lapply(censuses,`[[`,"census")
	N <- max(sapply(censusRows,length))
	censusTable <- do.call(rbind,lapply(censusRows,function(cr) c(cr,rep(0,N-length(cr)))))
	colnames(censusTable)[-(1:3)] <- 2:(N-2)

	lambdas <- sapply(censuses,`[[`,"lambda")

	complexities <- do.call(rbind,lapply(censuses,`[[`,"complexity"))

	###########################
	# Plot read counts vs complexities
	#############################

	outfile <- paste0(outdir,"complexities.pdf")
	pdf(outfile,5,5)
	layout(rbind(c(2,4),c(1,3)),heights=c(2,8),widths=c(8,2))
	op <- par(mar=c(5,4,0,0)+.1)
	plot(complexities,pch=".",log="xy",ylab="unique variants per AA-change")
	grid()
	par(mar=c(0,4,1,0)+.1)
	hist(
		log10(complexities$cpm),breaks=50,main="",axes=FALSE,
		col="gray",border=NA
	)
	# axis(1,c(-15,-11,-7,-3,1))
	axis(2)
	par(mar=c(5,0,0,1)+.1)
	yhist <- hist(log10(complexities$complexity),breaks=50,plot=FALSE)
	with(yhist,plot(
		NA,xlim=c(0,max(counts)),ylim=range(breaks),
		xlab="Frequency",ylab="",main="",axes=FALSE
	))
	axis(1)
	# axis(2,c(0,1,2))
	with(yhist,rect(0,breaks[-length(breaks)],counts,breaks[-1],col="gray",border=NA))
	par(op)
	invisible(dev.off())


	############################
	# Plot each tile's census
	#############################

	#function to find a good plot layout for a given number of tiles
	layoutFactors <- function(x) {
		sx <- sqrt(x)
		a <- ceiling(sx)
		b <- if (floor(sx)*a < x) a else floor(sx)
		c(a,b)
	}

	cr <- layoutFactors(nrow(censusTable))

	plotCensus <- function(i) {
		depth <- ns1.depth[i,2]
		barplot(censusTable[i,],
			ylim=c(0,100),border=NA,ylab="% reads",xlab="mutant classes",
			col=c("firebrick3","gray",rep("darkolivegreen3",N-2)),
			main=sprintf("Tile #%d: %.02fM reads",i,depth)
		)
		grid(NA,NULL)
		midx <- mean(par("usr")[1:2])
		text(midx,60,bquote(lambda == .(lambdas[[i]])))
		text(midx,40,sprintf("%.02f%% indels",censusTable[i,"indel"]))
	}

	outfile <- paste0(outdir,"tileCensus.pdf")
	pdf(outfile,2*cr[[1]],2.5*cr[[2]])
	layout(t(matrix(1:(cr[[1]]*cr[[2]]),ncol=cr[[1]])))
	op <- par(las=3,mar=c(6,4,3,1)+.1)
	invisible(lapply(1:nrow(censusTable),plotCensus))
	par(op)
	invisible(dev.off())


	###########################
	# Extrapolate overall variant census
	#############################


	#extrapolate overall indel rate
	# = 1 minus prob of zero indels in all tilesr
	# overallIndel <- 100*(1-dbinom(0,nrow(censusTable),mean(censusTable[,"indel"])/100))
	overallIndel <- 100*(1-(1-mean(censusTable[,"indel"])/100)^nrow(censusTable))
	#and overall #mut per clone
	# = mean of lambdas times number of tiles = sum of lambdas
	overallLambda <- sum(lambdas)

	overallCensus <- c(indel=overallIndel,setNames(dpois(0:8,overallLambda),0:8)*(100-overallIndel))

	outfile <- paste0(outdir,"extrapolatedCensus.pdf")
	pdf(outfile,5,5)
	op <- par(las=3,mar=c(6,4,3,1)+.1)
	barplot(overallCensus,
		ylim=c(0,100),border=NA,ylab="% reads",xlab="mutant classes",
		col=c("firebrick3","gray",rep("darkolivegreen3",length(overallCensus)-2)),
		main="Extrapolation for whole CDS"
	)
	grid(NA,NULL)
	midx <- mean(par("usr")[1:2])
	text(midx,60,bquote(lambda == .(overallLambda)))
	text(midx,40,sprintf("%.02f%% indels",overallCensus[["indel"]]))
	par(op)
	invisible(dev.off())

	## Plot distribution of lambdas
	# plot(function(x)dnorm(x,))
	# abline(v=jitter(lambdas),col="gray")


	###################
	# Coverage heatmap
	###################

	#parse all variant CPMs of all tiles and combine in one table
	allSingleCPMs <- do.call(rbind,lapply(1:nrow(ns1.depth), function(depth.row) {
		sample.id <- ns1.depth[depth.row,"sample"]
		sample.depth <- ns1.depth[depth.row,2]

		singleCounts <- read.delim(
			paste0(dataDir,"mutationCallfile/",sample.id,"AAchange.txt"),
			header=FALSE
		)
		colnames(singleCounts) <- c("wtaa","pos","mutaa","wtcodon","mutcodon","cpm")
		singleCounts
	}))

	#condense the table to unique amino acid changes
	mutcpms <- tapply(allSingleCPMs$cpm, with(allSingleCPMs,paste0(wtaa,pos,mutaa)), sum)
	allCPMs <- data.frame(
		mutname=names(mutcpms),
		wtaa=substr(names(mutcpms),1,1),
		pos=as.integer(substr(names(mutcpms),2,nchar(names(mutcpms))-1)),
		mutaa=substr(names(mutcpms),nchar(names(mutcpms)),nchar(names(mutcpms))),
		cpm=mutcpms
	)

	#calculate divider positions between tiles
	ntiles <- nrow(tileRanges)
	tileDividers <- tileRanges[-ntiles,2]


	#length protein (= highest AA position found)
	start <- min(allCPMs$pos)
	end <- max(allCPMs$pos)
	#prepare list of all possible amino acids
	aas <- toChars("AVILMFYWRHKDESTNQCGP_")

	#Prepare a color mapping function for CPM values
	# - find highest CPM value (to map to darkest color)
	maxPower <- ceiling(log10(max(allCPMs$cpm)))
	cmap <- yogitools::colmap(c(-1,2,maxPower),c("white","orange","firebrick3"))

	#Prepare heatmap coordinates and associated colors
	x <- allCPMs$pos
	y <- sapply(allCPMs$mutaa,function(aa) 22-which(aas==aa))
	colvals <- cmap(log10(allCPMs$cpm))

	#Determine plot size
	mapwidth <- end-start+2
	plotwidth <- mapwidth/10 + 3
	
	#Draw the heatmap plot
	outfile <- paste0(outdir,"coverageHeatmap.pdf")
	pdf(outfile,plotwidth*2/3,4)
	# layout(cbind(1,2),widths=c(mapwidth/10+1,3))
	layout(
		rbind(0:ntiles+3,c(rep(1,ntiles),2)),
		widths=c(rep(mapwidth/ntiles/10+1,ntiles),3),
		heights=c(1,1)
	)
	op <- par(xaxs="i",yaxs="i",mar=c(5,4,0,0)+.1)
	plot(
		NA,xlim=c(start-1,end+1),ylim=c(0,22),
		xlab="AA position",ylab="AA residue",
		axes=FALSE
	)
	axis(1)
	axis(2,at=21:1,aas,las=2,cex.axis=0.7)
	rect(start-0.5,0.5,end+.5,21.5,col="gray",border=NA)
	rect(x-0.5,y-0.5,x+0.5,y+0.5,col=colvals,border=NA)
	abline(v=tileDividers+.5)
	par(op)
	#legend
	op <- par(mar=c(5,0,0,6))
	plot(
		NA,xlim=c(0,1),ylim=c(-2.5,maxPower+.5),
		axes=FALSE,xlab="",ylab=""
	)
	rect(0,(-2:maxPower)-0.5,1,(-2:maxPower)+0.5,border=NA,
		col=cmap(c(NA,(-1:maxPower)))
	)
	axis(4,at=-2:maxPower,c("0","< 0.1","1",10^(1:(maxPower-1)),paste(">",10^maxPower)),las=2)
	mtext("Rel. reads (per mil.)",side=4,line=4)
	par(op)
	#census plots
	op <- par(las=3,mar=c(6,6,3,3)+.1)
	invisible(lapply(1:nrow(censusTable),plotCensus))
	par(op)
	invisible(dev.off())


	logInfo("libraryQC function completed successfully.")

}
