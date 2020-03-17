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

#' converts legacy tileseq count data into the new count file format
#'
#' @param dataDir path to the existing data main directory
#' @param outdir path to the desired output directory
#' @param mut2funcFile path to the mut2func file. Defaults to <dataDir>/mut2func_info.csv
#' @param countFolder path the mutationCallfile folder. Defaults to <dataDir>/mutationCallfile/
#' @param logger An optional yogilogger object. Defaults to NULL, which just prints to stdout.
#' @return NULL. Results are written to file. 
#' @export
legacyCountConverter <- function(dataDir, outdir, 
	mut2funcFile=paste0(dataDir,"mut2func_info.csv"), 
	countFolder=paste0(dataDir,"mutationCallfile/"),
	logger=NULL) {
	
	library(hgvsParseR)
	# library(yogilog)
	library(yogitools)

	op <- options(stringsAsFactors=FALSE)

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
	logErr <- function(...) {
		if (!is.null(logger)) {
			logger$error(...)
		} else {
			do.call(cat,c("ERROR:",list(...),"\n"))
		}
	}


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

	# mut2funcFile <- paste0(dataDir,"mut2func_info.csv")
	if (!canRead(mut2funcFile)) {
		stop("Unable to find or read 'mut2func_info.csv' file!")
	}
	# seqDepthFile <- paste0(dataDir,"resultfile/sequencingDepth.csv")
	# if (!canRead(seqDepthFile)) {
	# 	stop("Unable to find or read 'sequencingDepth.csv' file!")
	# }

	parseReport <- function(reportFile) {
		lines <- scan(reportFile,what="character",sep="\n",quiet=TRUE)
		mreads <- as.numeric(extract.groups(lines[[7]],"is: (.+) million")[1,1])
		return(mreads*1e6)
	}


	#inverse lookup of sample IDs
	mut2func <- read.csv(mut2funcFile)
	tileRanges <- apply(extract.groups(colnames(mut2func)[-1],"X(\\d+)\\.(\\d+)"),2,as.integer)
	sampleIDs <- unique(sort(do.call(c,mut2func[,-1])))
	condRep <- extract.groups(mut2func[,1],"^(\\w+)(\\d+)$")
	# rownames(condRep) <- mut2func[,1]
	sampleMetadata <- as.df(lapply(sampleIDs,function(sid) {
		idx <- which(mut2func==sid,arr.ind=TRUE)
		reportFile <- paste0(countFolder,sid,"report.txt")
		depth <- parseReport(reportFile)
		list(
			sample=sid, tile=unique(idx[,2]), 
			condition=paste(unique(condRep[idx[,1],1]),collapse=","),
			replicate=paste(unique(condRep[idx[,1],2]),collapse=","),
			timepoint=1,depth=depth
		)
	}))

	processCountFile <- function(countFile,depth) {
		counts <- read.delim(countFile,header=FALSE)

		cbuilder <- new.hgvs.builder.c()
		hgvsc <- apply(counts,1,function(row) {
			wt <- strsplit(row[[4]],"\\|")[[1]]
			pos <- as.integer(strsplit(row[[2]],"\\|")[[1]])
			mut <- strsplit(row[[5]],"\\|")[[1]]
			#codon start indices
			cstart <- pos*3-2
			hgvscs <- mapply(function(wt,pos,mut,cstart) {
				#calculate differences between codons
				diffs <- sapply(1:3,function(i)substr(wt,i,i)!=substr(mut,i,i))
				ndiff <- sum(diffs)
				if (ndiff == 1) { #one difference is a SNV
					offset <- which(diffs)
					wtbase <- substr(wt,offset,offset)
					mutbase <- substr(mut,offset,offset)
					snvpos <- cstart+offset-1
					cbuilder$substitution(snvpos,wtbase,mutbase)
				} else if (ndiff > 1) { #multiple differences is a delIns
					cbuilder$delins(cstart,cstart+2,mut)
				} else {
					stop("mutation must differ from wt!")
				}
			},wt,pos,mut,cstart)
			if (length(hgvscs) > 1) {
				do.call(cbuilder$cis,as.list(hgvscs))
			} else {
				hgvscs
			}
		})
		data.frame(HGVS=hgvsc,count=round(counts[,6]*depth/1e6))
	}


	# #Sample: 14
	# #Tile: 3
	# #Condition: nonselect
	# #Replicate: 1
	# #Timepoint: 1
	# #Read-depth: 30156
	# HGVS,count
	# c.21A>C,34
	# c.[13_15delinsTTG;29G>C],10
	# c.[8T>G;18_19del],8
	# c.123_125delinsCCG,42
	# ...

	invisible(apply(sampleMetadata,1,function(row) {
		sid <- as.integer(row[1])

		multi <- processCountFile(paste0(countFolder,sid,"MultipleMut.txt"),as.integer(row["depth"]))
		single <- processCountFile(paste0(countFolder,sid,"AAchange.txt"),as.integer(row["depth"]))
		header <- mapply(
			function(name,value) sprintf("#%s: %s",name,value), 
			name=c("Sample","Tile","Condition","Replicate","Timepoint","Read-depth"),
			value=row
		)

		con <- file(paste0(outdir,"counts_sample",sid,".csv"),open="w")
		writeLines(header,con)
		write.csv(rbind(multi,single),con,row.names=FALSE,quote=FALSE)
		close(con)

	}))

	options(op)

	return(NULL)
}

#' analyze tileseq counts from legacy pipeline
#'
#' This analysis function performs the following steps for each mutagenesis region:
#' 1. Construction of HGVS variant descriptor strings.
#' 2. Collapsing equivalent codons into amino acic change counts.
#' 3. Error regularization at the level of pre- and post-selection counts.
#' 4. Quality-based filtering filtering based on "Song's rule".
#' 5. Fitness score calculation and error propagation.
#' 6. Secondary error regularization at the level of fitness scores.
#' 7. Determination of synonymous and nonsense medians and re-scaling of fitness scores.
#' 8. Flooring of negative scores and adjustment of associated error.
#' 9. Output in MaveDB format.
#' 
#' @param countfile the path to the "rawData.txt" file produced by the legacy pipeline.
#' @param regionfile the path to a csv file describing the mutagenesis regions. Must contain columns 
#'  'region', start', 'end', 'syn', 'stop', i.e. the region id, the start position, end position, and
#'  and optional synonymous and stopm mean overrides.
#' @param outdir path to desired output directory
#' @param logger a yogilogger object to be used for logging (or NULL for simple printing)
#' @param inverseAssay a boolean flag to indicate that the experiment was done with an inverse assay
#'       i.e. protein function leading to decreased fitness. Defaults to FALSE
#' @param pseudoObservations The number of pseudoObservations to use for the Baldi&Long regularization.
#'       Defaults to 2.
#' @param conservativeMode Boolean flag. When turned on, pseudoObservations are not counted towards 
#'       standard error and the first round of regularization uses pessimistic error estimates.
#' @return nothing. output is written to various files in the output directory
#' @export
analyzeLegacyTileseqCounts <- function(countfile,regionfile,outdir,logger=NULL,
	inverseAssay=FALSE,pseudoObservations=2,conservativeMode=TRUE) {

	library(hgvsParseR)
	# library(yogilog)
	library(yogitools)

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
	logErr <- function(...) {
		if (!is.null(logger)) {
			logger$error(...)
		} else {
			do.call(cat,c("ERROR:",list(...),"\n"))
		}
	}

	##############
	# Read and validate input data
	##############

	canRead <- function(filename) file.access(filename,mode=4) == 0
	stopifnot(
		canRead(countfile),
		canRead(regionfile)
	)

	#parse the input files
	rawCounts <- read.delim(countfile)
	regions <- read.delim(regionfile, sep=",")

	#check that all required columns exist and contain the correct data types.
	stopifnot(
		c(
			"wt_codon","pos","mut_codon","wt_aa","mut_aa",
			"nonselect1","nonselect2","select1","select2",
			"controlNS1","controlNS2","controlS1","controlS2"
		) %in% colnames(rawCounts),
		all(apply(regions[,1:3],2,class)=="integer"),
		all(apply(regions[,4:5],2,class) %in% c("integer","numeric","logical")),
		c("region","start","end","syn","stop") %in% colnames(regions)
	)

	#make sure output directory path (outdir) ends with a "/"
	if (!grepl("/$",outdir)) {
		outdir <- paste0(outdir,"/")
	}
	#and if it doesn't exist, create it
	if (!dir.exists(outdir)) {
		dir.create(outdir,recursive=TRUE)
	}


	##################
	# Build HGVS variant descriptor strings
	##################

	#for nucleotide level descriptors
	cbuilder <- new.hgvs.builder.c()
	#for protein-level descriptors
	pbuilder <- new.hgvs.builder.p(aacode=3)

	#for nucleotide level descriptors
	hgvsc <- apply(rawCounts,1,function(row) {
		pos <- as.numeric(row["pos"])
		#codon start indices
		cstart <- pos*3-2
		wt <- row["wt_codon"]
		mut <- row["mut_codon"]
		#calculate differences between codons
		diffs <- sapply(1:3,function(i)substr(wt,i,i)!=substr(mut,i,i))
		ndiff <- sum(diffs)
		if (ndiff == 1) { #one difference is a SNV
			offset <- which(diffs)
			wtbase <- substr(wt,offset,offset)
			mutbase <- substr(mut,offset,offset)
			snvpos <- cstart+offset-1
			return(cbuilder$substitution(snvpos,wtbase,mutbase))
		} else if (ndiff > 1) { #multiple differences is a delIns
			return(cbuilder$delins(cstart,cstart+2,mut))
		} else {
			stop("mutation must differ from wt!")
		}
	})

	hgvsp <- apply(rawCounts,1,function(row) {
		pos <- as.numeric(row["pos"])
		wt <- row["wt_aa"]
		mut <- row["mut_aa"]
		if (mut=="_") mut <- "*" #correct stop character
		if (wt == mut) {
			return(pbuilder$synonymous(pos,wt))
		} else {
			return(pbuilder$substitution(pos,wt,mut))
		}
	})

	logInfo(sprintf(
		"Parsed data for %d codon variants (%d protein variants),\n %d (%d) of which are missense.",
		length(hgvsc),length(unique(hgvsp)), 
		sum(rawCounts$annotation == "NONSYN"), length(unique(hgvsp[!grepl("(Ter|=)$",hgvsp)]))
	))

	#####################
	# Pre-processing and formatting input
	####################

	#extract and examine condition names and replicates
	conditions <- colnames(rawCounts)[-c(1:6)]
	condStruc <- extract.groups(conditions,"^(\\w+)(\\d+)$")
	condNames <- unique(condStruc[,1])
	repNames <- unique(condStruc[,2])
	condMatrix <- do.call(rbind,lapply(condNames,paste0,repNames))
	dimnames(condMatrix) <- list(condNames,repNames)

	#TODO: Throw error if conditions are missing

	logInfo(sprintf(
		"Detected %d replicates for %d conditions: %s",
		length(repNames), length(unique(condNames)), paste(condNames,collapse=", ")
	))


	#apply pre-filter on variants
	rawmsd <- do.call(cbind,lapply(condNames,function(cond) {
		msd <- t(apply(rawCounts[,condMatrix[cond,]],1,function(xs){
			c(mean(xs,na.rm=TRUE), sd(xs,na.rm=TRUE))
		}))
		colnames(msd) <- paste0(cond,c(".mean",".sd"))
		msd
	}))
	flagged1 <- with(as.data.frame(rawmsd), {
		controlNS.mean + 3*controlNS.sd >= nonselect.mean
	})
	flagged2 <- with(as.data.frame(rawmsd), {
		controlS.mean + 3*controlS.sd >= select.mean
	})

	logInfo(sprintf(
"Filtering out %d variants (=%.02f%%):
%d (=%.02f%%) due to likely sequencing error.
%d (=%.02f%%) due to likely bottlenecking.",
		sum(flagged1|flagged2), 100*sum(flagged1|flagged2)/length(flagged1|flagged2),
		sum(flagged1), 100*sum(flagged1)/length(flagged1),
		sum(flagged2)-sum(flagged1), 100*(sum(flagged2)-sum(flagged1))/length(flagged2)
	))
	rawCountsFiltered <- rawCounts[!(flagged1|flagged2),]
	hgvsc <- hgvsc[!(flagged1|flagged2)]
	hgvsp <- hgvsp[!(flagged1|flagged2)]

	logInfo(sprintf(
		"Data remains for %d codon variants (%d protein variants),\n %d (%d) of which are missense.",
		length(hgvsc),length(unique(hgvsp)), 
		sum(rawCountsFiltered$annotation == "NONSYN"), length(unique(hgvsp[!grepl("(Ter|=)$",hgvsp)]))
	))

	# logInfo(sprintf(
	# 	"Data remains for for %d variants covering %d amino acid changes",
	# 	length(hgvsc),length(unique(hgvsp))
	# ))

	

	#collapse codons into unique AA changes
	logInfo("Collapsing variants by outcome...")
	combiCounts <- as.df(tapply(1:nrow(rawCountsFiltered),hgvsp,function(is) {
		# mut <- unique(hgvsp[is])
		hgvsp <- hgvsp[is[[1]]]
		hgvsc <- paste(unique(hgvsc[is]),collapse=" ")
		c(list(hgvsp=hgvsp,hgvsc=hgvsc),colSums(rawCountsFiltered[is,conditions],na.rm=TRUE))
	},simplify=FALSE))

	logInfo("Parsing variant strings...")
	combiCountMuts <- parseHGVS(combiCounts$hgvsp)
	combiCountMuts$type[which(combiCountMuts$variant=="Ter")] <- "nonsense"
	logInfo("Parsing complete.")



	###############################
	# Plot replicate correlations
	###############################

	logInfo("Plotting replicate correlations...")

	panel.cor <- function(x, y,...){
		usr <- par("usr"); on.exit(par(usr)); par(usr=c(0,1,0,1))
		r <- cor(fin(cbind(x,y)))[1,2]
		txt <- sprintf("R = %.02f",r)
		# cex.cor <- 0.8/strwidth(txt)
		text(0.5, 0.5, txt)
	}

	labels <- paste0("rep.",1:ncol(condMatrix))	
	imgSize <- max(4,ncol(condMatrix))

	pdfFile <- paste0(outdir,"replicateCorrelations_nonselect.pdf")
	pdf(pdfFile,imgSize,imgSize)
	# Plot non-select
	if (ncol(condMatrix) > 2) {
		pairs(
			combiCounts[,condMatrix["nonselect",]],
			lower.panel=panel.cor,pch=".",labels=labels,
			main="Non-select counts"
		)
	} else {
		plotTitle <- sprintf(
			"Non-select counts R = %.03f",
			cor(fin(combiCounts[,condMatrix["nonselect",]]))[1,2]
		)
		plot(
			combiCounts[,condMatrix["nonselect",]],
			main=plotTitle
		)
	}
	invisible(dev.off())

	# Plot phi after collapsing codons
	lfcRepsCombi <- log10((combiCounts[,condMatrix["select",]] - combiCounts[,condMatrix["controlS",]]) /
		(combiCounts[,condMatrix["nonselect",]] - combiCounts[,condMatrix["controlNS",]]))

	pdfFile <- paste0(outdir,"replicateCorrelations_logPhi.pdf")
	pdf(pdfFile,imgSize,imgSize)
	if (ncol(condMatrix) > 2) {
		pairs(
			lfcRepsCombi,pch=".",lower.panel=panel.cor,labels=labels,
			main="Select/Non-select Log-ratios"
		)
	} else {
		plotTitle <- sprintf(
			"Select/Non-select Log-ratios R = %.03f",
			cor(fin(lfcRepsCombi))[1,2]
		)
		plot(
			lfcRepsCombi, main=plotTitle,
			xlab=expression(log(phi[1])), ylab=expression(log(phi[2]))
		)
	}
	
	invisible(dev.off())

	##############
	# Iterate over regions
	##############
	for (region.i in 1:nrow(regions)) {

		regionStart <- regions[region.i,"start"]
		regionEnd <- regions[region.i,"end"]
		region.rows <- with(combiCountMuts,which(start >= regionStart & start <= regionEnd))
		region.syn <- regions[region.i,"syn"]
		region.stop <- regions[region.i,"stop"]

		if (length(region.rows) < 1) {
			logWarn(sprintf("\n\nNo variants found for region %d! Skipping...",region.i))
			next
		} else {
			logInfo(sprintf("\n\nProcessing region %d. (pos. %d-%d)",region.i,regionStart,regionEnd))
		}

		localCombiCounts <- combiCounts[region.rows,]
		localCombiCountMuts <- combiCountMuts[region.rows,]

		#########
		# Regularize SD of individual counts
		#########

		#Baldi & Long's formula
		bnl <- function(pseudo.n,n,model.sd,empiric.sd) {
			sqrt((pseudo.n * model.sd^2 + (n - 1) * empiric.sd^2)/(pseudo.n + n - 2))
		}

		#a helper function to construct an appropriately sized useful pseudocount 
		# pseudocount <- function(xs) if (all(xs==0)) 1e-4 else min(xs[xs!=0],na.rm=TRUE)/10
		ps <- 1e-4

		#calculate replicate means and coefficient of variation (CV) for each condition
		# (we use CV here instead of stdev, because that's what anti-correlates with read depth
		#  therefore, this will be the subject of our regression below)
		mcv <- do.call(cbind,lapply(condNames,function(cond) {
			mcv <- t(apply(localCombiCounts[,condMatrix[cond,]],1,function(xs){
				m <- mean(xs,na.rm=TRUE)
				# ps <- pseudocount(m)
				m <- m+ps
				c(m, (ps+sd(xs,na.rm=TRUE))/m)
			}))
			colnames(mcv) <- paste0(cond,c(".mean",".cv"))
			mcv
		}))

		#draw scatterplots for means vs CV
		pdfFile <- paste0(outdir,"region",region.i,"_regularizationInput.pdf")
		pdf(pdfFile,6,9)
		op <- par(mfrow=c(3,2))
		with(as.data.frame(mcv),{
			hist(log10(select.mean),breaks=100,col="gray",border=NA,main="")
			hist(log10(nonselect.mean),breaks=100,col="gray",border=NA,main="")
			plot(nonselect.mean,nonselect.cv,log="x",pch=".")
			plot(nonselect.mean,select.cv,log="x",pch=".")
			plot(nonselect.mean,controlNS.cv,log="x",pch=".")
			plot(controlNS.mean,controlNS.cv,log="x",pch=".")
		})
		par(op)
		invisible(dev.off())

		exportCoefficients <- function(z) {
			coef <- coefficients(z)
			coef <- coef[order(abs(coef),decreasing=TRUE)]
			paste(mapply(function(x,y)sprintf("%s = %.02f",x,y),x=names(coef),y=coef),collapse="; ")
		}

		#build regression input matrix by transforming to logspace and adding pseudocounts
		model.in <- log10(mcv)
		model.in <- as.data.frame(cbind(
			model.in,
			log.ratio=model.in[,"select.mean"]-model.in[,"nonselect.mean"]
		))

		#for each condition, perform the regularization
		regul <- do.call(cbind,lapply(condNames,function(cond) {
			#construct the regression formula
			input.column <- paste0(cond,".cv")
			mean.column <- paste0(cond,".mean")
			#extract empirical cv
			empiric.cv <- mcv[,input.column]
			input.column.id <- which(colnames(model.in)==input.column)
			regr.formula <- as.formula(paste0(input.column,"~."))
			#identify duplicated input columns, so they can be excluded
			dupli <- setdiff(which(apply(model.in,2,function(x) {
				all(x == model.in[,input.column])
			})),input.column.id)
			.model.in <- if (length(dupli) > 0) model.in[,-dupli] else model.in
			#perform regression
			z <- lm(regr.formula,data=.model.in)
			logInfo(paste("Regression coefficents for",cond,":",exportCoefficients(z)))
			#calculate prior (log(cv)) from regressor
			prior.lcv <- predict(z)
			logInfo(sprintf("Prior PCC=%.02f",cor(10^prior.lcv,empiric.cv)))
			#calculate bayesian regularization
			observations <- length(repNames)
			bayes.cv <- bnl(pseudoObservations,observations,10^prior.lcv,empiric.cv)
			#calculate pessimistic regularization
			pessim.cv <- mapply(max,10^prior.lcv,empiric.cv)
			#calculate the 1% quantile stdev (excluding the pseudocounts) as a minimum floor
			minCut <- quantile(empiric.cv[empiric.cv > 1e-5],0.01)
			#apply that minimum flooring to the pessimisting value
			pessim.cv <- sapply(pessim.cv,function(x) if (x < minCut) minCut else x)
			#return the result
			means <- mcv[,mean.column]
			out <- cbind(
				10^prior.lcv, (10^prior.lcv)*means,
				bayes.cv, bayes.cv*means,
				pessim.cv, pessim.cv*means
			)
			colnames(out) <- paste0(cond,c(
				".prior.cv",".prior.sd",
				".bayes.cv",".bayes.sd",
				".pessim.cv",".pessim.sd"
			))
			return(out)
		}))

		combiCountsReg <- cbind(localCombiCounts,mcv,regul)


		###################
		#Plot regularization results
		###################
		pdfFile <- paste0(outdir,"region",region.i,"_regularizationRound1.pdf")
		pdf(pdfFile,6,12)
		op <- par(mfrow=c(4,2))
		for (cond in condNames) {
			topoScatter(
				log10(mcv[,paste0(cond,".cv")]),
				log10(regul[,paste0(cond,".prior.cv")]),
				resolution=80,
				xlab=expression("Empiric"~log[10](sigma/mu)),
				ylab=expression("Prior"~log[10](sigma/mu)),
				main=cond
			)
			abline(0,1,col="gray",lty="dashed")
			topoScatter(
				log10(mcv[,paste0(cond,".cv")]),
				log10(regul[,paste0(cond,".bayes.cv")]),
				resolution=80,
				xlab=expression("Empiric"~log[10](sigma/mu)),
				ylab=expression("Regularized"~log[10](sigma/mu))
			)
			abline(0,1,col="gray",lty="dashed")
		}
		par(op)
		invisible(dev.off())

		logInfo("First round of regularization complete.")


		#####################
		# Quality checks
		####################

		#TODO: plot spatial distribution and histogram of sequencing error
		#


		#############
		#Filter based on high sequencing error
		# In this version, we just try to catch any straggles that weren't identified before regularization
		############
		#Song's rule: Filter out anything where the nonselect count is smaller than the WT control plus three SDs.
		#(That is, where the nonselect count could be explained by sequencing error)
		if (conservativeMode) {
			flagged <- with(combiCountsReg, controlNS.mean + 3*controlNS.pessim.sd >= nonselect.mean)
		} else {
			flagged <- with(combiCountsReg, controlNS.mean + 3*controlNS.bayes.sd >= nonselect.mean)
		}
		logInfo(sprintf(
			"Post-regularization: Filtering out an additional %d variants (=%.02f%%) due to likely sequencing error.",
			sum(flagged), 100*sum(flagged)/nrow(combiCountsReg)
		))
		combiCountsFiltered <- combiCountsReg[!flagged,]


		############
		#Calculate raw scores
		############

		#Helper function: Taylor approximation for variance of ratio
		approx.ratio.var <- function(mr,ms,vr,vs,covrs) {
			(mr^2)/(ms^2) * (vr/(mr^2) - 2*(covrs/(mr*ms)) + (vs/ms^2))
		}

		#pseudocounts
		c.ps <- 1e-4
		phivar.ps <- 1e-4

		rawScores <- as.df(lapply(1:nrow(combiCountsFiltered),function(i) {
			#mean numerator
			mnum <- combiCountsFiltered[i,"select.mean"]-combiCountsFiltered[i,"controlS.mean"]
			if (mnum < c.ps) mnum <- c.ps #apply flooring, so the wt control can't result in negatives counts
			#mean denominator
			mden <- combiCountsFiltered[i,"nonselect.mean"]-combiCountsFiltered[i,"controlNS.mean"]
			#variances
			if (conservativeMode) {
				#variance numerator
				vnum <- combiCountsFiltered[i,"select.pessim.sd"]^2+combiCountsFiltered[i,"controlS.pessim.sd"]^2
				#variance denominator
				vden <- combiCountsFiltered[i,"nonselect.pessim.sd"]^2+combiCountsFiltered[i,"controlNS.pessim.sd"]^2
			} else {
				vnum <- combiCountsFiltered[i,"select.bayes.sd"]^2+combiCountsFiltered[i,"controlS.bayes.sd"]^2
				vden <- combiCountsFiltered[i,"nonselect.bayes.sd"]^2+combiCountsFiltered[i,"controlNS.bayes.sd"]^2
			}
			#covariance of numerator and denominator
			covnumden <- cov(
				unlist(combiCountsFiltered[i,condMatrix["select",]])
				-unlist(combiCountsFiltered[i,condMatrix["controlS",]]),
				unlist(combiCountsFiltered[i,condMatrix["nonselect",]])
				-unlist(combiCountsFiltered[i,condMatrix["controlNS",]])
			)
			#Use helper function to estimate variance for phi
			phivar <- approx.ratio.var(mnum,mden,vnum,vden,covnumden)
			#As this is based on a Taylor approximation, we can sometimes get freaky negative values
			if (phivar < phivar.ps) phivar <- phivar.ps

			phi <- mnum/mden
			phisd <- sqrt(phivar)
			list(
				hgvsp=combiCountsFiltered[i,"hgvsp"],
				hgvsc=combiCountsFiltered[i,"hgvsc"],
				phiPrime=min(combiCountsFiltered[i,c("nonselect1","nonselect2")]),
				mean.c=mnum,
				mean.phi=phi,
				sd.phi=phisd,
				mean.lphi=log10(phi),
				sd.lphi=abs(phisd/(log(10)*phi))
			)
		}))

		#Filter out for selection bottlenecking
		flagged <- rawScores$mean.c <= c.ps
		logInfo(sprintf(
			"Filtering out %d variants (=%.02f%%) due to likely selection bottlenecking.",
			sum(flagged), 100*sum(flagged)/nrow(rawScores)
		))
		rawScores <- rawScores[!flagged,]


		logInfo(sprintf(
			"Data remains for %d protein-level variants, %d of which are missense.",
			nrow(rawScores),sum(!grepl("Ter|=$",rawScores$hgvsp))
		))


		##############
		# 2nd round of regularization: This time for SD of scores
		##############

		# #regression input matrix
		# splinemat <- with(rawScores,data.frame(
		# 	logsd=log10(sd.lphi),
		# 	logminbc=log10(phiPrime+1),
		# 	lphi=mean.lphi
		# ))
		# #run regression
		# z <- lm(logsd ~.,data=splinemat)
		#calculate prior
		# priorSD <- 10^predict(z)

		logInfo("Starting second round of regularization.")

		#regression input matrix
		splinemat <- with(rawScores,data.frame(
			logcv=log10(sd.phi/mean.phi),
			logminbc=log10(phiPrime+.1)#,
			# lphi=mean.lphi
		))
		#run regression
		z <- lm(logcv ~.,data=splinemat)
		#calculate prior
		priorSD <- 10^predict(z)*rawScores$mean.phi

		#apply regularization
		observations <- length(repNames)
		bayesSD <- bnl(pseudoObservations,observations,priorSD,rawScores$sd.phi)
		rawScores$bsd.phi=bayesSD
		rawScores$bsd.lphi=with(rawScores,abs(bsd.phi/(log(10)*mean.phi)))
		rawScores$df=if(conservativeMode) observations else observations+pseudoObservations

		#####################
		# Filter out broken values if any
		#####################
		broken <- which(is.na(rawScores$mean.lphi) | is.infinite(rawScores$mean.lphi))
		if (length(broken) > 0) {
			logWarn(sprintf("Values for %d variants could not be calculated and were removed!",
				length(broken)
			))
			rawScores <- rawScores[-broken,]
		}

		################
		# Plot results of regularization
		################
		pdfFile <- paste0(outdir,"region",region.i,"_regularizationRound2.pdf")
		pdf(pdfFile,7,7)
		op <- par(mfrow=c(2,2))
		#Plot phiPrime vs SD
		with(rawScores[rawScores$sd.lphi < 1,],topoScatter(phiPrime+1,sd.lphi,log="x",maxFreq=35,thresh=3,
			resolution=40, xlab="Non-select count (per M.)", ylab=expression(sigma)
		))
		#Plot score vs SD
		with(rawScores,topoScatter(mean.lphi,sd.lphi,log="y",pch=20,resolution=40, 
			xlab="Fitness score",ylab=expression(sigma),maxFreq=35,thresh=3
		))
		#Plot phiPrime vs regSD
		if (sum(rawScores$bsd.lphi < 1) > 1) {
			with(rawScores[rawScores$bsd.lphi < 1,],topoScatter(phiPrime+1,bsd.lphi,log="x",maxFreq=35,thresh=3,
				resolution=40, xlab="Non-select count (per M.)", 
				ylab=expression("Bayesian Regularized"~sigma)
			))
		}
		# abline(0,1,col="gray",lty="dashed")
		#Plot Empiric vs regularized
		with(rawScores,topoScatter(sd.lphi,bsd.lphi,resolution=60,maxFreq=30,log="xy",
			xlab=expression("Empiric"~sigma),ylab=expression("Bayesian Regularized"~sigma)
		))
		abline(0,1,col="gray",lty="dashed")
		par(op)
		invisible(dev.off())


		#################
		# If the functional assay is based on inverse fitness, invert the scores
		#################
		if (inverseAssay) {
			rawScores$mean.lphi <- -rawScores$mean.lphi
			#stdev remains unchanged as sqrt(((-1)^2)) = 1
		}

		#################
		# Estimate Modes of synonymous and stop
		#################

		sdCutoff <- 0.3

		####################3
		#TODO: Require minimum amount of filter passes rather than just 1
		#TODO: Try gaussian mixture models with two underlying distributions?
		#########################
		#if we can't find any syn/stop below the cutoff, increase the cutoff to 1
		if (with(rawScores,!any(grepl("Ter$",hgvsp) & bsd.lphi < sdCutoff) ||
			!any(grepl("=$",hgvsp) & bsd.lphi < sdCutoff) )) {
			sdCutoff <- 1
			#if we still can't find any below the new cutoff, get rid of it altogether
			if (with(rawScores,!any(grepl("Ter$",hgvsp) & bsd.lphi < sdCutoff) ||
				!any(grepl("=$",hgvsp) & bsd.lphi < sdCutoff) )) {
				sdCutoff <- Inf
			}
		}

		modes <- with(rawScores,{	
			stops <- mean.lphi[which(grepl("Ter$",hgvsp) & bsd.lphi < sdCutoff)]
			syns <- mean.lphi[which(grepl("=$",hgvsp) & bsd.lphi < sdCutoff)]
			c(stop=median(stops,na.rm=TRUE),syn=median(syns,na.rm=TRUE))
		})


		#################
		#Plot Syn vs stop
		#################

		plotStopSyn <- function(data,title,modes) {
			with(data,{
				stop.is <- which(grepl("Ter$",hgvsp))
				syn.is <- which(grepl("=$",hgvsp))
				miss.is <- setdiff(1:nrow(data),c(stop.is,syn.is))
				br <- seq(floor(min(mean.lphi)),ceiling(max(mean.lphi)),.1)

				hist(mean.lphi[miss.is],col="gray",breaks=br,xlab=expression(log(phi)),border=NA,main=title)
				hist(mean.lphi[stop.is],add=TRUE,col=colAlpha("firebrick3",.5),breaks=br)
				hist(mean.lphi[syn.is],add=TRUE,col=colAlpha("darkolivegreen3",.5),breaks=br)
				abline(v=modes,col=c("firebrick3","darkolivegreen3"),lty="dashed")
			})
		}

		pdfFile <- paste0(outdir,"scalingQC_region",region.i,".pdf")
		pdf(pdfFile)
		op <- par(mfrow=c(2,1))
		plotStopSyn(rawScores,"Unfiltered",modes)
		legend("topright",c("missense","synonymous","stop"),fill=c("gray","darkolivegreen3","firebrick3"))
		if (any(rawScores$bsd.lphi < sdCutoff)) {
			plotStopSyn(rawScores[rawScores$bsd.lphi < sdCutoff,],
				bquote(sigma["regularized"] < .(sdCutoff)),
				modes
			)
		}
		par(op)
		invisible(dev.off())


		#################
		# Use syn/stop modes to scale scores
		#################

		logInfo(sprintf(
			"Scaling to synonymous (log(phi)=%.02f) and nonsense (log(phi)=%.02f) medians.",modes[["syn"]],modes[["stop"]]
		))

		#if manual overrides for the synonymous and stop modes were provided, use them
		if (!is.na(region.syn)) {
			logInfo(sprintf(
				"Using manual override (=%.02f) for synonmous mode instead of automatically determined value (=%.02f).",region.syn,modes[["syn"]]
			))
			modes[["syn"]] <- region.syn
		}
		if (!is.na(region.stop)) {
			logInfo(sprintf(
				"Using manual override (=%.02f) for stop mode instead of automatically determined value (=%.02f).",region.stop,modes[["stop"]]
			))
			modes[["stop"]] <- region.stop
		}

		if (modes[["stop"]] >= modes[["syn"]]) {
			logErr("Stop mode is not below synonymous mode! Cannot normalize scores!")
			next
		}

		#apply the scaling
		denom <- modes[["syn"]]-modes[["stop"]]
		scoreMat <- with(rawScores,{
			sd <- bsd.lphi/denom
			score <- (mean.lphi - modes[["stop"]])/denom
			cbind(
				score=score,sd=sd,se=sd/sqrt(df)
			)
		})
		scores <- cbind(rawScores,scoreMat)


		##################
		# Floor negatives and fix their excessive variances
		##################

		if (!inverseAssay) {
			#the target null-like score towards which we will shift these values
			targetScore <- 0
			#the quantile for which we want to keep the p-value fixed
			quantile <- 1
			#the row numbers containing the cases to be fixed
			toFix <- which(scores$score < targetScore)
		} else {
			#if we're dealing with an inverse assay, we have to apply a ceiling instead of flooring
			#the target functional (but dead) score towards which we will shift these values
			targetScore <- 1
			#the quantile for which we want to keep the p-value fixed
			quantile <- 0
			#the row numbers containing the cases to be fixed
			toFix <- which(scores$score > targetScore)
		}
		#the equivalent sds of a normal distribution with the target mean based on the above area
		equivalent.sds <- with(scores[toFix,], sd*(quantile-targetScore)/(quantile-score))
		#apply the fixed values to the table
		scores$score.unfloored <- scores$score
		scores$sd.unfloored <- scores$sd
		scores$se.unfloored <- scores$se
		scores$score[toFix] <- targetScore
		scores$sd[toFix] <- equivalent.sds
		scores$se[toFix] <- equivalent.sds/sqrt(scores$df[toFix])
		
		logInfo(sprintf(
			"Flooring adjusted the values of %d variant scores",length(toFix)
		))


		##############################
		# Draw a histogram of standard errors in the dataset
		#############################
		pdfFile <- paste0(outdir,"errorProfile_region",region.i,".pdf")
		pdf(pdfFile,5,5)
		hist(
			log10(scores$se),col="gray",border=NA,breaks=50,
			main="Standard Error distribution",
			xlab=expression(log[10](sigma[bar(x)]))
		)
		dev.off()

		
		#################
		# Write output to file
		#################

		#detailed output of all intermediate results
		outfile <- paste0(outdir,"detailed_scores_region",region.i,".csv")
		write.csv(scores,outfile,row.names=FALSE)

		#protein-level MaveDB output
		mavedbProtScores <- scores[,c("hgvsp","score","sd","se")]
		colnames(mavedbProtScores) <- c("hgvs_pro","score","sd","se")
		outfile <- paste0(outdir,"mavedb_scores_perAA_region",region.i,".csv")
		write.csv(mavedbProtScores,outfile,row.names=FALSE)

		logInfo(sprintf(
			"Exported data for %d protein-level variants, %d of which are missense.",
			nrow(scores),sum(!grepl("Ter|=$",scores$hgvsp))
		))

		#nucleotide-level MaveDB output
		mavedbNuclScores <- do.call(rbind,lapply(1:nrow(scores),function(i) {
			hgvsc <- strsplit(scores[i,"hgvsc"]," ")[[1]]
			data.frame(
				hgvs_nt=hgvsc,
				hgvs_pro = scores[i,"hgvsp"],
				score = scores[i,"score"],
				sd = scores[i,"sd"],
				se = scores[i,"se"]
			)
		}))
		outfile <- paste0(outdir,"mavedb_scores_perNt_region",region.i,".csv")
		write.csv(mavedbNuclScores,outfile,row.names=FALSE)

		logInfo(sprintf(
			"Exported data for %d codon-level variants, %d of which are missense.",
			nrow(mavedbNuclScores),sum(!grepl("Ter|=$",mavedbNuclScores$hgvs_pro))
		))


	}

	###################
	# Join the results into a single files
	###################

	joinFiles <- function(filenameBase) {
		joined <- do.call(rbind,lapply(
			paste0(outdir,filenameBase,"_region",regions$region,".csv"),
			function(rfile) {
				if (file.exists(rfile)) {
					return(read.csv(rfile))
				} else {
					return(NULL)
				}
			}
		))
		outfile <- paste0(outdir,filenameBase,".csv")
		write.csv(joined,outfile,row.names=FALSE)
	}

	joinFiles("detailed_scores")
	joinFiles("mavedb_scores_perNt")
	joinFiles("mavedb_scores_perAA")

}



#' QC plots for library coverage
#' @param dataDir path to the existing data main directory
#' @param outdir path to the output directory
#' @param mut2funcFile path to the mut2func file. Defaults to <dataDir>/mut2func_info.csv
#' @param logger optional yogilogger object. Defaults to NULL, in which case messages will go to stdout.
#' @param mc.cores number of CPU cores to use in parallel
#' @return NULL. Results are written to file.
#' @export
libraryQCLegacy <- function(dataDir,outdir,
	mut2funcFile=paste0(dataDir,"mut2func_info.csv"),
	logger=NULL,mc.cores=8) {

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

	# mut2funcFile <- paste0(dataDir,"mut2func_info.csv")
	if (!canRead(mut2funcFile)) {
		stop("Unable to find or read 'mut2func_info.csv' file!")
	}
	# seqDepthFile <- paste0(dataDir,"resultfile/sequencingDepth.csv")
	# if (!canRead(seqDepthFile)) {
	# 	stop("Unable to find or read 'sequencingDepth.csv' file!")
	# }


	mut2func <- read.csv(mut2funcFile)
	tileRanges <- apply(extract.groups(colnames(mut2func)[-1],"X(\\d+)\\.(\\d+)"),2,as.integer)


	#Interpret sequencing depth file and obtain sample information
	# seqDepthTable <- read.csv(seqDepthFile)
	# rxGroups <- extract.groups(seqDepthTable$SampleID,"SampleID(\\d+)_(\\w+)(\\d+)")
	# seqDepthTable$sample <- as.integer(rxGroups[,1])
	# seqDepthTable$condition <- rxGroups[,2]
	# seqDepthTable$replicate <- as.integer(rxGroups[,3])

	# #Isolate first replicate of nonselect samples
	# ns1.rows <- with(seqDepthTable,which(condition=="NS" & replicate == 1))
	# ns1.depth <- seqDepthTable[ns1.rows,]

	nsSamples <- unlist(mut2func[which(mut2func$tiles=="nonselect1"),-1])
	ns1.depth <- as.df(lapply(nsSamples, function(sid) {
		reportFile <- paste0(dataDir,"mutationCallfile/",sid,"report.txt")
		reportLines <- scan(reportFile,what="character",sep="\n",quiet=TRUE)
		depth <- as.numeric(extract.groups(reportLines[grep("sequencing depth",reportLines)],"is: (.+) million")[,1])
		list(sample=sid,depth=depth)
	}))


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
		multiVars <- mapply(function(w,p,m) paste0(w,p,m), wtaas, positions, mutaas, SIMPLIFY=FALSE)
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

		#Next, read the AAchange file. Note that these are not true single-mutants, but rather
		#the marginal counts for each amino acid change
		marginalCounts <- read.delim(
			paste0(dataDir,"mutationCallfile/",sample.id,"AAchange.txt"),
			header=FALSE
		)
		colnames(marginalCounts) <- c("wtaa","pos","mutaa","wtcodon","mutcodon","cpm")

		#collapse by aa-change
		marginalVariants <- with(marginalCounts,mapply(function(w,p,m) paste0(w,p,m), wtaa, pos, mutaa))
		marginalTotals <- hash()
		for (i in 1:length(marginalVariants)) {
			v <- marginalVariants[[i]]
			if (hash::has.key(v,marginalTotals)) {
				marginalTotals[[v]] <- marginalTotals[[v]] + marginalCounts[[i,"cpm"]]
			} else {
				marginalTotals[[v]] <- marginalCounts[[i,"cpm"]]
			}
		}

		#list of unique aa changes
		uniqueVars <- keys(marginalTotals)

		#calculate the true single mutant CPMs by subtracting the multi-mutant marginal CPMs
		# from the total marginal CPMs
		singletonCounts <- sapply(uniqueVars, function(uvar) {
			if (has.key(uvar,multiTotals)) {
				marginalTotals[[uvar]] - multiTotals[[uvar]]
			} else {
				marginalTotals[[uvar]] 
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
		#Unfortunately, we have to treat the indels separately, as we have no information on the 
		# cross-talk between indels and missense/nonsense variants.
		census <- c(
			indel = indelCPM,
			# WT = 1e6 - (indelCPM + sum(marginalCounts$cpm)),
			WT = 1e6 - (sum(marginalCounts$cpm)),
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
	if (all(cr==c(2,1))) {
		pdf(outfile,5,10)
	} else {
		pdf(outfile,2*cr[[1]],2.5*cr[[2]])
	}
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


#legacy converter for mut2func_info
convertLegacyMut2Func <- function(infile,outfile) {

	library(yogitools)

	op <- options(stringsAsFactors=FALSE); on.exit(options(op))

	m2f <- read.csv(infile,header=FALSE)

	tiles <- cbind(1:(ncol(m2f)-1),do.call(rbind,strsplit(unlist(m2f[1,-1]),"-")))
	colnames(tiles) <- c("id","start","end")

	condRep <- extract.groups(m2f[-1,1],"^(\\w+)(\\d+)$")

	samples <- do.call(rbind,lapply(2:nrow(m2f), function(i) {
		cbind(id=unlist(m2f[i,-1]),tile=1:(ncol(m2f)-1),cond=condRep[i-1,1],tp=1,rep=condRep[i-1,2])
	}))

	con <- file(outfile,open="w")
	write.table(tiles,con,row.names=FALSE,quote=FALSE,sep="\t")
	writeLines("",con)
	write.table(samples,con,row.names=FALSE,quote=FALSE,sep="\t")
	close(con)

}

