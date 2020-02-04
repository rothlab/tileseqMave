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
libraryQC <- function(dataDir,paramFile=paste0(dataDir,"parameters.json"),logger=NULL,mc.cores=6) {

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
	#extract time stamp
	timeStamp <- extract.groups(latestCountDir,"/(\\d{4}-\\d{2}-\\d{2}-\\d{2}-\\d{2}-\\d{2})")

	#create a matching output directory
	outDir <- paste0(dataDir,timeStamp,"_QC/")
	dir.create(outDir,recursive=TRUE,showWarnings=FALSE)


	#identify nonselect conditions
	nsConditions <- unique(with(as.data.frame(params$conditions$definitions),{
		`Condition 2`[which(Relationship == "is_selection_for")]
	}))


	logInfo("Reading count data")

	allCountFile <- paste0(latestCountDir,"/allCounts.csv")
	allCounts <- read.csv(allCountFile)

	marginalCountFile <- paste0(latestCountDir,"/marginalCounts.csv")
	marginalCounts <- read.csv(marginalCountFile)

	logInfo("Interpreting variant descriptors")



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
	tileStarts <- params$tiles[,"Start AA"]
	inferTiles <- function(cct) {
		sapply(cct$pos,function(pos) max(which(tileStarts <= pos)))
	}
	allTiles <- pbmclapply(allSplitChanges,inferTiles,mc.cores=mc.cores)	
	marginalTiles <- inferTiles(marginalSplitChanges)

	#Quick Hack: Pull frameshift positions from initial HGVS strings
	hackpos <- as.integer(extract.groups(allCounts$hgvsp,"(\\d+)fs")[,1])
	allTiles[which(!is.na(hackpos))] <- sapply(
		hackpos[which(!is.na(hackpos))],
		function(pos) max(which(tileStarts <= pos))
	)

	#TODO: 
	# * Iterate over (possibly multiple) nonselects and analyze:
	#   * Coverage map and census for each tile
	#   * Extrapolated census
	#   * Complexity analysis
	#   * Stacked barplots for missense/syn/stop/indel/frameshift (for each tile)

	#helper function to find the WT control condition for a given condition
	findWTCtrl <- function(params,cond) {
		defs <- params$conditions$definitions
		defs[which(defs[,3] == cond & defs[,2] == "is_wt_control_for"),1]
	}


	#Iterate over (possibly multiple) nonselects and analyze
	for (nsCond in nsConditions) {

		#pull out nonselect condition and average over replicates
		nsReps <- sprintf("%s.t%s.rep%s.frequency",nsCond,params$timepoints[1,1],1:params$numReplicates)
		nsMarginalMeans <- rowMeans(marginalCounts[,nsReps],na.rm=TRUE)
		nsAllMeans <- rowMeans(allCounts[,nsReps],na.rm=TRUE)

		#if WT controls are present, average over them as well and subtract from nonselect
		wtCond <- findWTCtrl(params,nsCond)
		if (length(wtcond) > 0) {
			wtReps <- sprintf("%s.t%s.rep%s.frequency",wtCond,params$timepoints[1,1],1:params$numReplicates)

			wtMarginalMeans <- rowMeans(marginalCounts[,wtReps],na.rm=TRUE)
			nsMarginalMeans <- mapply(function(nsf,wtf) {
				max(0,nsf-wtf)
			},nsf=nsMarginalMeans,wtf=wtMarginalMeans)

			wtAllMeans <- rowMeans(allCounts[,wtReps],na.rm=TRUE)
			nsAllMeans <- mapply(function(nsf,wtf) {
				max(0,nsf-wtf)
			},nsf=nsAllMeans,wtf=wtAllMeans)

		} else {
			logWarn("No WT control defined for",nsCond,": Unable to perform error correction!")
		}

		#build a simplified marginal frequency table from the results
		simplifiedMarginal <- cbind(marginalSplitChanges,freq=nsMarginalMeans,tile=marginalTiles)

		#generate nuclotide bias analysis output
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
		dev.off()

		#calculate global filters that can be combined with tile-specific filters below
		isFrameshift <- allCounts$aaChanges=="fs"
		isInFrameIndel <- sapply(allSplitChanges, function(changes) {
			any(nchar(changes[,3]) != 3)
		})
		numChanges <- sapply(allSplitChanges, nrow)
		maxChanges <- max(numChanges)
		
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

		#calculate best fitting lambda per census
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

		#FIXME: THIS NEEDS TO BE DONE PER REGION! NOT PER WHOLE CDS!!!
		#Extrapolate overall census
		#overall frameshift/indel rates
		# = 1 minus prob of zero frameshifts/indels in all tiles
		overallFS <- 1-(1-mean(tileCensi[,"fs"],na.rm=TRUE))^nrow(tileCensi)
		overallIndel <- 1-(1-mean(tileCensi[,"indel"],na.rm=TRUE))^nrow(tileCensi)
		#and overall #mut per clone
		# = mean of lambdas times number of tiles = sum of lambdas
		# This assumes that tiles are of equal size. Weighted sum over relative tile size would be better
		overallLambda <- sum(tileLambdas[,"lambda"],na.rm=TRUE)

		overallCensus <- c(fs=overallFS,indel=overallIndel,
			setNames(dpois(0:maxChanges,overallLambda), 0:maxChanges)*(1-(overallIndel+overallFS))
		)

		#Plot overall census
		pdf(paste0(outDir,nsCond,"_census.pdf"),5,5)
		plotCensus(overallCensus,overallLambda,main="Extrapolation for whole CDS")
		dev.off()

	}

	options(op)
	return(NULL)
}



plotCensus <- function(census,lambda,main="") {
	op <- par(las=2)
	barplot(
		census*100, ylim=c(0,100),
		col=c(rep("firebrick2",2),"gray",rep("darkolivegreen3",length(census)-3)),
		border=NA, main=main,
		xlab="mutant classes",ylab="% reads"
	)
	midx <- mean(par("usr")[1:2])
	midy <- mean(par("usr")[3:4])
	text(midx,midy,bquote(lambda == .(round(lambda,digits=3))))
	par(op)
}

nucleotideBiasAnalysis <- function(simplifiedMarginal,tile,draw=TRUE) {
	
	counters <- replicate(2,array(0, dim=c(4,4,3),
		dimnames=list(toChars("ACGT"), toChars("ACGT"), 1:3)
	), simplify=FALSE)
	names(counters) <- c("single","multi")

	#iterate over each variant entry
	for (i in 1:nrow(simplifiedMarginal)) {
		#skip frameshifts and deletions
		if (any(is.na(simplifiedMarginal[i,])) || nchar(simplifiedMarginal[i,"to"]) < 3) {
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

