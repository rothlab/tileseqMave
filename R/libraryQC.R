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
libraryQC <- function(dataDir,paramFile=paste0(dataDir,"parameters.json"),logger=NULL,mc.cores=6,srOverride=FALSE) {

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
	params <- parseParameters(paramFile,srOverride=srOverride)

	#find counts folder
	subDirs <- list.dirs(dataDir,recursive=FALSE)
	countDirs <- subDirs[grepl("_mut_call$",subDirs)]
	if (length(countDirs) == 0) {
		stop("No mutation call output found!")
	}
	latestCountDir <- sort(countDirs,decreasing=TRUE)[[1]]
	#extract time stamp
	timeStamp <- extract.groups(latestCountDir,"/(\\d{4}-\\d{2}-\\d{2}-\\d{2}-\\d{2}-\\d{2})")

	#create a matching output directory
	outDir <- paste0(dataDir,timeStamp,"_QC/")
	dir.create(outDir,recursive=TRUE,showWarnings=FALSE)


	#identify nonselect conditions
	nsConditions <- getNonselects(params)
	# nsConditions <- unique(with(as.data.frame(params$conditions$definitions),{
	# 	`Condition 2`[which(Relationship == "is_selection_for")]
	# }))

	#if this is a pure QC run, there are no conditions declared as nonselect, so
	# we treat *all* conditions as nonselect.
	if (length(nsConditions) == 0) {
		if (nrow(params$conditions$definitions) > 0) {
			nsConditions <- params$conditions$definitions[,3]
		} else {
			nsConditions <-  params$conditions$names
		}
	}


	logInfo("Reading count data")

	allCountFile <- paste0(latestCountDir,"/allCounts.csv")
	allCounts <- read.csv(allCountFile)

	marginalCountFile <- paste0(latestCountDir,"/marginalCounts.csv")
	marginalCounts <- read.csv(marginalCountFile)

	depthTableFile <-  paste0(latestCountDir,"/sampleDepths.csv")
	depthTable <- read.csv(depthTableFile)


	logInfo("Interpreting variant descriptors. This may take some time...")

	#Extract affected positions and codon details
	splitCodonChanges <- function(ccs) {
		df <- as.data.frame(
			extract.groups(ccs,"(\\D*)(\\d+)(\\D*)"),
			stringsAsFactors=FALSE
		)
		colnames(df) <- c("from","pos","to")
		df$pos <- as.integer(df$pos)
		df
	}
	allSplitChanges <- pbmclapply(
		strsplit(allCounts$codonChanges,"\\|"), 
		splitCodonChanges, mc.cores=mc.cores
	)
	marginalSplitChanges <- splitCodonChanges(marginalCounts$codonChange)

	#Infer tile assignments
	logInfo("Calculating tile assignments. This may take some time...")
	tileStarts <- params$tiles[,"Start AA"]
	inferTiles <- function(cct) {
		sapply(cct$pos,function(pos) max(which(tileStarts <= pos)))
	}
	allTiles <- pbmclapply(allSplitChanges,inferTiles,mc.cores=mc.cores)	
	marginalTiles <- inferTiles(marginalSplitChanges)

	# #Quick Hack: Pull frameshift positions from initial HGVS strings
	# hackpos <- as.integer(extract.groups(allCounts$hgvsp,"(\\d+)fs")[,1])
	# allTiles[which(!is.na(hackpos))] <- sapply(
	# 	hackpos[which(!is.na(hackpos))],
	# 	function(pos) max(which(tileStarts <= pos))
	# )
	# hackpos <- as.integer(extract.groups(marginalCounts$hgvsp,"(\\d+)fs")[,1])
	# marginalTiles[which(!is.na(hackpos))] <- sapply(
	# 	hackpos[which(!is.na(hackpos))],
	# 	function(pos) max(which(tileStarts <= pos))
	# )
	

	#helper function to find the WT control condition for a given condition
	# findWTCtrl <- function(params,cond) {
	# 	defs <- params$conditions$definitions
	# 	defs[which(defs[,3] == cond & defs[,2] == "is_wt_control_for"),1]
	# }


	################################################################################
	# ITERATE OVER (POSSIBLY MULTIPLE) NONSELECT CONDITIONS AND ANALYZE SEPARATELY
	################################################################################

	for (nsCond in nsConditions) {

		logInfo("Processing",nsCond)

		#########################################
		# CALCULATE NONSELECT MEANS AND NORMALIZE
		#########################################
		#pull out nonselect condition and average over replicates
		nsReps <- sprintf("%s.t%s.rep%s.frequency",nsCond,params$timepoints[1,1],1:params$numReplicates[[nsCond]])
		if (params$numReplicates[[nsCond]] > 1) {
			nsMarginalMeans <- rowMeans(marginalCounts[,nsReps],na.rm=TRUE)
			nsAllMeans <- rowMeans(allCounts[,nsReps],na.rm=TRUE)
		} else {
			nsMarginalMeans <- marginalCounts[,nsReps]
			nsAllMeans <- allCounts[,nsReps]
		}

		#if WT controls are present, average over them as well and subtract from nonselect
		wtCond <- getWTControlFor(nsCond,params)
		# wtCond <- findWTCtrl(params,nsCond)
		if (length(wtCond) > 0) {
			wtReps <- sprintf("%s.t%s.rep%s.frequency",wtCond,params$timepoints[1,1],1:params$numReplicates[[wtCond]])

			if (params$numReplicates[[wtCond]] > 1) {
				wtMarginalMeans <- rowMeans(marginalCounts[,wtReps],na.rm=TRUE)
			} else {
				wtMarginalMeans <- marginalCounts[,wtReps]
			}
			nsMarginalMeans <- mapply(function(nsf,wtf) {
				max(0,nsf-wtf)
			},nsf=nsMarginalMeans,wtf=wtMarginalMeans)

			if (params$numReplicates[[wtCond]] > 1) {
				wtAllMeans <- rowMeans(allCounts[,wtReps],na.rm=TRUE)
			} else {
				wtAllMeans <- allCounts[,wtReps]
			}
			nsAllMeans <- mapply(function(nsf,wtf) {
				max(0,nsf-wtf)
			},nsf=nsAllMeans,wtf=wtAllMeans)

		} else {
			logWarn("No WT control defined for",nsCond,": Unable to perform error correction!")
		}

		#build a simplified marginal frequency table from the results
		simplifiedMarginal <- cbind(marginalSplitChanges,freq=nsMarginalMeans,tile=marginalTiles)


		#################################
		# RUN NUCLEOTIDE BIAS ANALYSIS
		#################################
		logInfo("Checking nucleotide distribution")
		pdf(paste0(outDir,nsCond,"_nucleotide_bias.pdf"),8.5,11)
		opar <- par(mfrow=c(6,1))
		nuclRates <- lapply(params$tiles[,1], function(tile) {
			if (any(simplifiedMarginal$tile == tile)){
				nucleotideBiasAnalysis(
					simplifiedMarginal[which(simplifiedMarginal$tile == tile),],
					tile
				)
			} else {
				plot.new()
				rect(0,0,1,1,col="gray80",border="gray30",lty="dotted")
				text(0.5,0.5,"no data")
				mtext(paste0("Tile #",tile),side=4)
				return(NA)
			}
		})
		par(opar)
		invisible(dev.off())

		#calculate global filters that can be combined with tile-specific filters below
		isFrameshift <- grepl("fs$",allCounts$aaChanges)
		isInFrameIndel <- sapply(allSplitChanges, function(changes) {
			any(changes[,3] != "indel" & nchar(changes[,3]) != 3)
		})
		numChanges <- sapply(allSplitChanges, nrow)
		maxChanges <- max(numChanges)
		
		###########################
		# CALCULATE VARIANT CENSUS
		###########################

		logInfo("Calculating variant census")
		#Census per tile
		tileCensi <- do.call(rbind,lapply(params$tiles[,1], function(tile) {
			#filter for entries in given tile
			isInTile <- sapply(allTiles,function(tiles) !any(is.na(tiles)) && tile %in% tiles)
			if (!any(isInTile)) {
				return(c(fs=NA,indel=NA,WT=NA,setNames(rep(NA,maxChanges),1:maxChanges)))
			}
			# WT frequency
			wtFreq <- 1 - sum(nsAllMeans[isInTile])
			#frameshift freq
			fsFreq <- sum(nsAllMeans[isInTile & isFrameshift])
			#indel frequency
			indelFreq <- sum(nsAllMeans[isInTile & !isFrameshift & isInFrameIndel])
			#substitution mutation census
			isSubst <- which(isInTile & !isFrameshift & !isInFrameIndel)
			census <- tapply(nsAllMeans[isSubst], numChanges[isSubst],sum)
			census <- c(
				fs=fsFreq,
				indel=indelFreq,
				WT=wtFreq,
				sapply(as.character(1:maxChanges),function(n) if (n %in% names(census)) census[[n]] else 0)
			)
			return(census)
		}))

		#calculate lambda and Poisson fits for each census
		tileLambdas <- do.call(rbind,lapply(1:nrow(tileCensi), function(tile) {
			freqs <- tileCensi[tile,-c(1,2)]
			if (any(is.na(freqs))) {
				return(c(lambda=NA,rmsd=NA))
			}
			ns <- 1:length(freqs)-1
			lambda <- sum(ns*(freqs/sum(freqs)))
			rmsd <- sqrt(mean((freqs-dpois(ns,lambda))^2))
			return(c(lambda=lambda,rmsd=rmsd))
		}))

		####################################################
		# Extrapolate overall census per mutagenesis region
		####################################################

		logInfo("Extrapolating variant censi for each region")
		#determine which tiles are in which regions
		tilesPerRegion <- tilesInRegions(params)
		pdf(paste0(outDir,nsCond,"_census.pdf"),8.5,11)
		opar <- par(mfrow=c(4,3))
		regionCensi <- lapply(1:length(tilesPerRegion), function(ri) {

			relevantTiles <- tilesPerRegion[[ri]]
			#if there is no data for any of those tiles
			if (all(is.na(tileCensi[relevantTiles,"WT"]))) {
				plot.new()
				rect(0,0,1,1,col="gray80",border="gray30",lty="dotted")
				text(0.5,0.5,"no data")
				mtext(paste0("Extrapolation for Region #",ri))
				return(c(fs=NA,indel=NA,WT=NA,setNames(rep(NA,maxChanges),1:maxChanges)))
			}

			#extrapolate overall frameshift and indel rates
			overallFS <- 1-(1-mean(tileCensi[relevantTiles,"fs"],na.rm=TRUE))^length(relevantTiles)
			overallIndel <- 1-(1-mean(tileCensi[relevantTiles,"indel"],na.rm=TRUE))^length(relevantTiles)
			#Neat: Potential length differences between tiles cancel out, as the lambdas are not length-normalized!
			overallLambda <- sum(tileLambdas[relevantTiles,"lambda"],na.rm=TRUE)
			overallCensus <- c(fs=overallFS,indel=overallIndel,
				setNames(dpois(0:maxChanges,overallLambda), 0:maxChanges)*(1-(overallIndel+overallFS))
			)

			#Plot overall census
			plotCensus(overallCensus,overallLambda,main=paste0("Extrapolation for Region #",ri))
			return(overallCensus)
		})
		par(opar)
		invisible(dev.off())


		######################
		# BUILD COVERAGE MAP #
		######################
		logInfo("Building coverage map.")
		#load translation table
		data(trtable)
		#translate marginals and join by translation
		aaMarginal <- simplifiedMarginal[with(simplifiedMarginal,which(!is.na(pos) & nchar(to)==3)),]
		aaMarginal$toaa <- sapply(aaMarginal$to,function(codon)trtable[[codon]])
		aaMarginal <- as.df(with(aaMarginal,tapply(1:nrow(aaMarginal),paste0(pos,toaa), function(is) {
			list(
				from=trtable[[unique(from[is])]],
				pos=unique(pos[is]),
				to=unique(toaa[is]),
				tile=unique(tile[is]),
				freq=sum(freq[is])
			)
		})))

		#plan page layout
		tpr <- 4 #tiles per row
		rpp <- 3 #rows per page
		ntiles <- nrow(params$tiles)
		#build layout plan
		tileLayout <- lapply(seq(1,ntiles,tpr),function(i)i:min(ntiles,i+tpr-1))
		nrows <- length(tileLayout)
		pageLayout <- lapply(seq(1,nrows,rpp), function(i) tileLayout[i:min(nrows,i+rpp-1)])

		#build a matrix that indicate the grid positions of plot in order of drawing
		plotIndex <- do.call(rbind,lapply(1:rpp, function(ri) {
			bottom <- (ri-1)*5+1
			#position map on bottom of the row, and census plots for the corresponding tiles above
			rbind(bottom+1:tpr, rep(bottom,tpr))
		}))

		#now do the actual plotting.
		pdf(paste0(outDir,nsCond,"_coverage.pdf"),8.5,11)
		layout(plotIndex)
		invisible(lapply(pageLayout, function(tileSets) {
			#For each row on the current page...
			lapply(tileSets, function(tiles) {
				#plot coverage map
				startPos <- min(params$tiles[tiles,"Start AA"])
				endPos <- max(params$tiles[tiles,"End AA"])
				seps <- params$tiles[tiles,"Start AA"][-1]
				coverageSubmap(startPos,endPos,aaMarginal,seps)
				#TODO: Obtain number of reads in each tile and add to plot
				#and plot the corresponding censi
				lapply(tiles, function(tile) {
					depths <- with(depthTable,depth[
						Condition==nsCond & Time.point==params$timepoints[1,1] & Tile.ID==tile
					])
					plotCensus(
						tileCensi[tile,],
						lambda=tileLambdas[tile,"lambda"],
						d=if(length(depths)>0) min(depths,na.rm=TRUE) else NULL,
						main=paste0("Tile #",tile)
					)
				})
			})
		}))
		invisible(dev.off())

		#######################
		# Complexity analysis #
		#######################
		logInfo("Running complexity analysis.")
		ccList <- strsplit(allCounts$codonChanges,"\\|")
		#count the number of combinations in which each codon change occurs (="complexity")
		ccComplex <- hash()
		for (ccs in ccList) {
			for (cc in ccs) {
				if (has.key(cc,ccComplex)) {
					ccComplex[[cc]] <- ccComplex[[cc]]+1
				} else {
					ccComplex[[cc]] <- 1
				}
			}
		}
		#tabulate "complexity" vs marginal frequency
		cplxPlot <- do.call(rbind,lapply(1:nrow(simplifiedMarginal), function(i) with(simplifiedMarginal[i,],{
			cplx <- ccComplex[[paste0(from,pos,to)]]
			c(cplx=if (is.null(cplx)) NA else cplx,freq=freq,tile=tile)
		})))

		#draw the plot for each tile
		pdf(paste0(outDir,nsCond,"_complexity.pdf"),8.5,11)
		opar <- par(mfrow=c(3,2))
		invisible(tapply(1:nrow(cplxPlot),cplxPlot[,"tile"],function(is) {
			if (length(is) > 1) {
				tile <- unique(cplxPlot[is,3])
				if (!all(cplxPlot[is,2] == 0)) {
					plot(cplxPlot[is,1:2],log="xy",
						xlab="#unique contexts in tile",ylab="marginal frequency",
						main=paste0("Tile #",tile)
					)
				}
			}
		}))
		par(opar)
		invisible(dev.off())
		
		##############################
		# Well-measuredness analysis #
		##############################
		logInfo("Running Well-measuredness analysis")
		reachable <- reachableChanges(params)

		data(trtable)
		#add translations to marginal table
		simplifiedMarginal$toaa <- sapply(simplifiedMarginal$to, function(codon){
			if (codon=="indel") {
				"fs"
			} else if (nchar(codon) < 3) {
				"del"
			} else if (nchar(codon) == 3) {
				trtable[[codon]]
			} else {
				"ins"
			}
		})

		#filter out frameshifts and indels
		codonMarginal <- simplifiedMarginal[nchar(simplifiedMarginal$toaa)==1,]
		#add flag for SNVs
		codonMarginal$isSNV <- sapply(1:nrow(codonMarginal), function(i) with(codonMarginal[i,], {
			sum(sapply(1:3,function(j) substr(from,j,j) != substr(to,j,j)))==1
		}))
		#and collapse by aa change
		aaMarginal <- as.df(with(codonMarginal,tapply(1:length(pos),paste0(pos,toaa),function(is) {
			list(pos=unique(pos[is]),toaa=unique(toaa[is]),freq=sum(freq[is]),tile=unique(tile[is]))
		})))
		#add index for quick lookup
		aaMarginal$index <- with(aaMarginal,paste0(pos,toaa))

		#list of frequency thresholds to test
		thresholds <- 10^seq(-7,-2,0.1)

		#iterate over regions and analyze separately
		pdf(paste0(outDir,nsCond,"_wellmeasured.pdf"),8.5,11)
		opar <- par(mfrow=c(3,2))
		coverageCurves <- lapply(1:length(tilesPerRegion), function(ri) {

			tiles <- tilesPerRegion[[ri]]
			rStart <- min(params$tiles[tiles,"Start AA"])
			rEnd <- max(params$tiles[tiles,"End AA"])
			rLength <- rEnd-rStart+1

			#if there's no data, return an empty plot
			if (!any(aaMarginal$tile %in% tiles)) {
				plot.new()
				rect(0,0,1,1,col="gray80",border="gray30",lty="dotted")
				text(0.5,0.5,"no data")
				mtext(paste0("Region #",ri))
				return(NULL)
			}

			regionReachable <- reachable[with(reachable,pos >= rStart & pos <= rEnd),]
			regionReachable$index <- with(regionReachable,paste0(pos,mutaa))

			nreachableAA <- nrow(regionReachable)
			npossibleAA <- rLength * 21 #includes syn + stop options
			npossibleSNV <- rLength * 9 # = three possible changes at three positions per codon
			npossibleCodon <- rLength * 63 # = 4*4*4-1

			#calculate coverage curves
			coverageCurves <- as.df(lapply(thresholds, function(thr) {
				list(
					freqCutoff = thr,
					fracPossibleCodon = with(codonMarginal,sum(freq > thr & tile %in% tiles))/npossibleCodon,
					fracSNV=with(codonMarginal,sum(freq > thr & tile %in% tiles & isSNV))/npossibleSNV,
					fracPossibleAA = with(aaMarginal,sum(freq > thr & tile %in% tiles))/npossibleAA,
					fracReachableAA = with(aaMarginal,
						sum(freq > thr & tile %in% tiles & index %in% regionReachable$index)
					)/nreachableAA
				)
			}))

			#draw the plot
			plotcolors <- c(
				`codon changes`="black",`SNVs`="gray",
				`AA changes`="firebrick3",`SNV-reachable AA changes`="firebrick2"
			)
			plot(NA,type="n",
				log="x",xlim=range(thresholds),ylim=c(0,1),
				xlab="Marginal read frequency cutoff",
				ylab="Fraction of possible variants measured",
				main=paste0("Region #",ri)
			)
			for (i in 1:4) {
				lines(coverageCurves[,c(1,i+1)],col=plotcolors[[i]])
			}
			grid()
			legend("bottomleft",names(plotcolors),col=plotcolors,lty=1)

			return(coverageCurves)

		})
		par(opar)
		invisible(dev.off())
		
		#################################################
		# Stacked barplots for missense/syn/stop/indel/frameshift (for each tile)
		##################################################

		logInfo("Plotting mutation type breakdown per tile")
		#FIXME: Frameshifts have been accidentally joined into one entry
		simplifiedMarginal$fromaa <- sapply(simplifiedMarginal$from, 
			function(codon) if (is.na(codon)) "" else trtable[[codon]]
		)
		mutbreakdown <- do.call(cbind,lapply(1:nrow(params$tiles), function(tile) {
			if (!any(simplifiedMarginal$tile == tile)) {
				return(setNames(rep(NA,5),c("fs","indel","stop","syn","mis")))
			}
			marginalSubset <- simplifiedMarginal[simplifiedMarginal$tile == tile,]
			freqsums <- with(marginalSubset,{
				c(
					fs=sum(freq[toaa=="fs"]),
					indel=sum(freq[which(toaa !="fs" & nchar(to) != 3)]),
					stop=sum(freq[toaa=="*"]),
					syn=sum(freq[fromaa==toaa])
				)
			})
			freqsums[["mis"]] <- sum(marginalSubset$freq)-sum(freqsums)
			return(freqsums)
		}))
		colnames(mutbreakdown) <- params$tiles[,1]
		
		plotcols <- c(
			frameshift="firebrick3",indel="orange",nonsense="gold",
			synonymous="steelblue3",missense="darkolivegreen3"
		)
		pdf(paste0(outDir,nsCond,"_mutationtypes.pdf"),8.5,11)
		barplot(mutbreakdown,col=plotcols,border=NA,horiz=TRUE,
			xlab="sum of marginal frequencies",ylab="Tile",main="Mutation types"
		)
		legend("right",names(plotcols),fill=plotcols)
		invisible(dev.off())

	}

	options(op)
	logInfo("QC analysis complete.")
	return(NULL)
}

#helper function to draw a subsection of the coverage map.
coverageSubmap <- function(startPos,endPos,aaMarginal,seps=NULL,thresholds=c(1e-6,5e-5)) {
	mapwidth <- endPos-startPos+2
	aas <- toChars("AVILMFYWRHKDESTNQCGP*")
	#TODO: Adjust color thresholds based on sequencing depth
	cmap <- yogitools::colmap(c(log10(thresholds[[1]]),log10(thresholds[[2]]),0),c("white","orange","firebrick3"))

	#set a filter for positions in the selected part of the map
	is <- with(aaMarginal,which(pos >= startPos & pos <= endPos))
	#if there is no data, return an empty plot
	if (length(is) == 0) {
		plot.new()
		rect(0,0,1,1,col="gray80",border="gray30",lty="dotted")
		text(0.5,0.5,"no data")
		invisible(return(NULL))	
	}

	plotData <- as.df(lapply(is, function(i) with(aaMarginal[i,],list(
		x=pos,
		y=22-which(aas==to),
		colval=cmap(log10(freq))
	))))

	#prepare canvas
	op <- par(xaxs="i",yaxs="i",mar=c(5,4,0,1)+.1)
	plot(
		NA,xlim=c(startPos-1,endPos+1),ylim=c(0,22),
		xlab="AA position",ylab="AA residue",
		axes=FALSE
	)
	#add axis labels
	axis(1)
	axis(2,at=21:1,aas,las=2,cex.axis=0.7)
	#gray background
	rect(startPos-0.5,0.5,endPos+.5,21.5,col="gray",border=NA)
	#heatmap squares
	with(plotData,rect(x-0.5,y-0.5,x+0.5,y+0.5,col=colval,border=NA))
	#tile separator lines
	if (!is.null(seps)) {
		abline(v=seps-.5)
	}
	par(op)

	invisible(return(NULL))	
}

#helper function to find the tiles that belong to each region
tilesInRegions <- function(params) {
	lapply(1:nrow(params$regions), function(ri) {
		rs <- params$regions[ri,"Start AA"]
		re <- params$regions[ri,"End AA"]
		which(sapply(params$tiles[,"Start AA"], function(ts){
			ts >=rs && ts < re
		}))
	})
}

#helper function to plot a census dataset
plotCensus <- function(census,lambda,d=NULL,main="") {
	if (is.null(census) || any(is.na(census))) {
		plot.new()
		rect(0,0,1,1,col="gray80",border="gray30",lty="dotted")
		text(0.5,0.5,"no data")
		mtext(main)
		invisible(return(NULL))
	}

	op <- par(las=2)
	barplot(
		census*100, ylim=c(0,100),
		col=c(rep("firebrick2",2),"gray",rep("darkolivegreen3",length(census)-3)),
		border=NA, main=main,
		xlab="mutant classes",ylab="% reads"
	)
	grid(NA,NULL)
	u <- par("usr")
	midx <- mean(u[1:2])
	midy <- mean(u[3:4])
	uppermid <- 0.4*u[[3]]+0.6*u[[4]]
	lowermid <- 0.6*u[[3]]+0.4*u[[4]]
	text(midx,lowermid,bquote(lambda == .(round(lambda,digits=3))))
	if (!is.null(d)) {
		text(midx,uppermid,paste("#reads =",d))
	}
	par(op)
	invisible(return(NULL))
}

#helper function to perform nucleotide bias analysis and draw the corresponding plot
nucleotideBiasAnalysis <- function(simplifiedMarginal,tile,draw=TRUE) {
	
	counters <- replicate(2,array(0, dim=c(4,4,3),
		dimnames=list(toChars("ACGT"), toChars("ACGT"), 1:3)
	), simplify=FALSE)
	names(counters) <- c("single","multi")

	#iterate over each variant entry
	for (i in 1:nrow(simplifiedMarginal)) {
		#skip frameshifts and deletions
		if (any(is.na(simplifiedMarginal[i,])) 
			|| nchar(simplifiedMarginal[i,"to"]) < 3
			|| simplifiedMarginal[i,"to"] == "indel") {
			next
		}

		#calculate number of nucleotide changes
		ndiff <- sum(sapply(1:3,function(j) {
			substr(simplifiedMarginal[i,"from"],j,j) != substr(simplifiedMarginal[i,"to"],j,j)
		}))

		#skip insertions, otherwise classify as SNV or MNV
		if (ndiff == 0) {
			next
		} else if (ndiff == 1) {
			type <- "single"
		} else {
			type <- "multi"
		}

		#add to counters
		for (j in 1:3) {
			from <- substr(simplifiedMarginal[i,"from"],j,j)
			to <- substr(simplifiedMarginal[i,"to"],j,j)
			freq <- simplifiedMarginal[i,"freq"]
			if (from != to) {
				counters[[type]][from,to,j] <- counters[[type]][from,to,j] + freq
			}
		}
	}

	#normalize columns
	rates <- lapply(c("single","multi"),function(type){
		apply(counters[[type]],c(1,3),function(xs) xs/sum(xs))
	})
	names(rates) <- c("single","multi")

	if (draw) {
		#calculate plot coordinates and colors
		plotVals <- expand.grid(
			from=toChars("ACGT"),to=toChars("ACGT"),
			pos=1:3,type=c("single","multi"),stringsAsFactors=FALSE
		)
		cm <- colmap(c(0,1),c("white","steelblue3"))
		plotVals <- cbind(plotVals,as.df(lapply(1:nrow(plotVals),function(i) with(plotVals[i,],{
			list(
				x = which(toChars("ACGT")==to) + 5*(pos-1) + switch(type,single=0,multi=16),
				y = 5-which(toChars("ACGT")==from),
				val = rates[[type]][from,to,pos],
				color = cm(rates[[type]][from,to,pos])
			)
		}))))

		op <- par(mar=c(5,4,0,1)+.1)
		plot(NA,type="n",xlim=c(0,30),ylim=c(0,5),axes=FALSE,xlab="from",ylab="to")
		with(plotVals,rect(x-1,y-1,x,y,col=color,border=NA))
		with(plotVals,text(
			x-.5, y-.5, 
			sapply(round(val*100),function(v) {
				if (is.na(v)) "n.d." else paste0(v,"%")
			}), 
			cex=0.8
		))
		text(c(2,7,12,18,23,28),4.33,rep(c("First","Second","Third"),2),cex=0.9)
		text(c(7,23),4.66,c("SNV","MNV"),cex=0.9)
		axis(2,at=(4:1)-.5,toChars("ACGT"))
		axis(1,at=c(1:4,6:9,11:14,17:20,22:25,27:30)-.5,rep(toChars("ACGT"),6))
		mtext(paste0("Tile #",tile),side=4)
		par(op)
	}

	return(rates)
}

