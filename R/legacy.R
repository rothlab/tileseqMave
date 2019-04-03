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
#' @param logger a yogilogger object to be used for logging
#' @return nothing. output is written to various files in the output directory
#' @export
analyzeLegacyTileseqCounts <- function(countfile,regionfile,outdir,logger=NULL) {

	library(hgvsParseR)
	library(yogilog)
	library(yogitools)

	options(stringsAsFactors=FALSE)

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

	#countfile <- "/home/jweile/projects/ccbr2hgvs/HMGCR_S_resultfile/rawData.txt"
	#regionfile <- "/home/jweile/projects/ccbr2hgvs/HMGCR_S_resultfile/regions.txt"
	canRead <- function(filename) file.access(filename,mode=4) == 0
	stopifnot(
		canRead(countfile),
		canRead(regionfile)
	)

	rawCounts <- read.delim(countfile)
	regions <- read.delim(regionfile)

	stopifnot(
		c(
			"wt_codon","pos","mut_codon","wt_aa","mut_aa",
			"nonselect1","nonselect2","select1","select2",
			"controlNS1","controlNS2","controlS1","controlS2"
		) %in% colnames(rawCounts),
		all(apply(regions[,1:3],2,class)=="integer"),
		all(apply(regions[,4:5],2,class)=="numeric"),
		c("region","start","end","syn","stop") %in% colnames(regions)
	)

	#make sure outdir ends with a "/"
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
		"Parsed data for %d variants covering %d amino acid changes",
		length(hgvsc),length(unique(hgvsp))
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

	#collapse codons into unique AA changes
	logInfo("Collapsing variants by outcome...")
	combiCounts <- as.df(tapply(1:nrow(rawCounts),hgvsp,function(is) {
		# mut <- unique(hgvsp[is])
		hgvsp <- hgvsp[is[[1]]]
		hgvsc <- paste(unique(hgvsc[is]),collapse=" ")
		c(list(hgvsp=hgvsp,hgvsc=hgvsc),colSums(rawCounts[is,conditions]))
	},simplify=FALSE))

	logInfo("Parsing variant strings...")
	combiCountMuts <- parseHGVS(combiCounts$hgvsp)
	combiCountMuts$type[which(combiCountMuts$variant=="Ter")] <- "nonsense"
	logInfo("Parsing complete.")

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
			logWarn(sprintf("No variants found for region %d! Skipping...",region.i))
			next
		} else {
			logInfo(sprintf("Processing region %d. (pos. %d-%d)",region.i,regionStart,regionEnd))
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
		pseudocount <- function(xs) if (all(xs==0)) 1e-4 else min(xs[xs!=0],na.rm=TRUE)/10

		#calculate replicate means and coefficient of variation (CV) for each condition
		# (we use CV here instead of stdev, because that's what anti-correlates with read depth
		#  so, this will be the subject of our regression below)
		mcv <- do.call(cbind,lapply(condNames,function(cond) {
			mcv <- t(apply(localCombiCounts[,condMatrix[cond,]],1,function(xs){
				m <- mean(xs,na.rm=TRUE)
				ps <- pseudocount(m)
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
			hist(nonselect.mean,breaks=100,col="gray",border=NA,main="")
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
			bayes.cv <- bnl(2,2,10^prior.lcv,empiric.cv)
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


		#####################
		# Quality checks
		####################

		#TODO: plot spatial distribution and histogram of sequencing error
		#


		#############
		#Filter based on high sequencing error
		############
		#Song's rule: Filter out anything where the nonselect count is smaller than the WT control plus three SDs.
		#(That is, where the nonselect count could be explained by sequencing error)
		flagged <- with(combiCountsReg, controlNS.mean + 3*controlNS.pessim.sd >= nonselect.mean)
		logInfo(sprintf(
			"Filtering out %d variants (=%.02f%%) due to likely sequencing error.",
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

		rawScores <- as.df(lapply(1:nrow(combiCountsFiltered),function(i) {
			#mean numerator
			mnum <- combiCountsFiltered[i,"select.mean"]-combiCountsFiltered[i,"controlS.mean"]
			if (mnum < 1e-3) mnum <- 1e-3 #apply flooring, so the wt control can't result in negatives counts
			#mean denominator
			mden <- combiCountsFiltered[i,"nonselect.mean"]-combiCountsFiltered[i,"controlNS.mean"]
			#variance numerator
			vnum <- combiCountsFiltered[i,"select.pessim.sd"]^2+combiCountsFiltered[i,"controlS.pessim.sd"]^2
			#variance denominator
			vden <- combiCountsFiltered[i,"nonselect.pessim.sd"]^2+combiCountsFiltered[i,"controlNS.pessim.sd"]^2
			#covariance of numerator and denominator
			covnumden <- cov(
				unlist(combiCountsFiltered[i,c("select1","select2")])
				-unlist(combiCountsFiltered[i,c("controlNS1","controlNS2")]),
				unlist(combiCountsFiltered[i,c("nonselect1","nonselect2")])
				-unlist(combiCountsFiltered[i,c("controlS1","controlS2")])
			)
			#Use helper function to estimate variance for phi
			phivar <- approx.ratio.var(mnum,mden,vnum,vden,covnumden)
			#As this is based on a Taylor approximation, we can sometimes get freaky negative values
			if (phivar < 1e-4) phivar <- 1e-4

			phi <- mnum/mden
			phisd <- sqrt(phivar)
			list(
				hgvsp=combiCountsFiltered[i,"hgvsp"],
				hgvsc=combiCountsFiltered[i,"hgvsc"],
				phiPrime=min(combiCountsFiltered[i,c("nonselect1","nonselect2")]),
				mean.phi=phi,
				sd.phi=phisd,
				mean.lphi=log10(phi),
				sd.lphi=abs(phisd/(log(10)*phi))
			)
		}))

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
		bayesSD <- bnl(2,2,priorSD,rawScores$sd.phi)
		rawScores$bsd.phi=bayesSD
		rawScores$bsd.lphi=with(rawScores,abs(bsd.phi/(log(10)*mean.phi)))
		rawScores$df=4

		#####################
		# Filter out broken values if any
		#####################
		broken <- which(is.na(rawScores$mean.lphi) | is.infinite(rawScores$mean.lphi))
		if (length(broken) > 0) {
			logWarn(sprintf("Values for %d variants could not be calculated!",
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
		if (any(rawScores$bsd.lphi < 1)) {
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
		# Estimate Modes of synonymous and stop
		#################

		sdCutoff <- 0.3

		if (with(rawScores,!any(grepl("Ter$",hgvsp) & bsd.lphi < sdCutoff) ||
			!any(grepl("=$",hgvsp) & bsd.lphi < sdCutoff) )) {
			sdCutoff <- 1
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

		#the target null-like score towards which we will shift these values
		targetScore <- -0.01
		#the row numbers containing the cases to be fixed
		toFix <- which(scores$sd > .3 & scores$score < targetScore)
		#the area under the part of the normal distribution that exceeds zero
		ps <- with(scores,pnorm(0,score[toFix],se[toFix],lower.tail=FALSE))
		#the equivalent sds of a normal distribution with the target mean based on the above area
		equivalent.sds <- (-targetScore)/(qnorm(1-ps))
		#apply the fixed values to the table
		scores$score[toFix] <- targetScore
		scores$sd[toFix] <- equivalent.sds
		scores$se[toFix] <- equivalent.sds/sqrt(scores$df[toFix])
		
		logInfo(sprintf(
			"Flooring adjusted the values of %d variant scores",length(toFix)
		))

		
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


