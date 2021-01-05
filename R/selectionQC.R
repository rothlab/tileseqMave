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
#' @param countDir input directory for counts, defaults to subdirectory with latest timestamp ending in _mut_count.
#' @param scoreDir input directory for scores, defaults to subdirectory with latest timestamp ending in _scores.
#' @param outDir output directory, defaults to name of input directory with _QC tag attached.
#' @param paramFile input parameter file. defaults to <dataDir>/parameters.json
#' @return NULL. Results are written to file.
#' @export
selectionQC <- function(dataDir,countDir=NA, scoreDir=NA, outDir=NA, 
                        paramFile=paste0(dataDir,"parameters.json"),
                        srOverride=FALSE) {


	op <- options(stringsAsFactors=FALSE)

	library(yogitools)
	library(hgvsParseR)
	library(pbmcapply)
	# library(optimization)

	#make sure data exists and ends with a "/"
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


	logInfo("Reading parameters from",normalizePath(paramFile))
	params <- withCallingHandlers(
		parseParameters(paramFile,srOverride=srOverride),
		warning=function(w)logWarn(conditionMessage(w))
	)
	
	#find counts folder
	if (is.na(countDir)) {
	  latest <- latestSubDir(parentDir=dataDir,pattern="_mut_call$|mut_count$")
	  countDir <- latest[["dir"]]
	  timeStamp <- latest[["timeStamp"]]
	  runLabel <- latest[["label"]]
	} else { #if custom input dir was provided
	  #make sure it exists
	  if (!dir.exists(countDir)) {
	    stop("Input count folder ",countDir," does not exist!")
	  }
	  #try to extract a timestamp and label
	  lt <- extract.groups(countDir,"([^/]+_)?(\\d{4}-\\d{2}-\\d{2}-\\d{2}-\\d{2}-\\d{2})")[1,]
	  if (!any(is.na(lt))) {
	    runLabel <- lt[[1]]
	    timeStamp <- lt[[2]]
	  } else {
	    #if none can be extracted, use current time and no tag
	    timeStamp <- format(Sys.time(), "%Y-%m-%d-%H-%M-%S")
	    runLabel <- ""
	  }
	}
	#make sure it ends in "/"
	if (!grepl("/$",countDir)) {
	  countDir <- paste0(countDir,"/")
	}
	
	#if no score directory was provided
	if (is.na(scoreDir)) {
	  #try to guess its name
	  if (grepl("_mut_count/$",countDir)) {
	    scoreDir <- sub("_mut_count/$","_scores/",countDir)
	  } else {
	    scoreDir <- sub("/$","_scores/",countDir)
	  }
	  if (!dir.exists(scoreDir)) {
	    stop("No matching score directory found for ",countDir,"!\n",
	         "(Expecting ",scoreDir,")\n",
	         "Use --scores option to define explicitly."
	    )
	  }
	} else {
	  if (!dir.exists(scoreDir)) {
	    stop("Score directory ",scoreDir," does not exist!")
	  }
	}
	#make sure it ends in "/"
	if (!grepl("/$",scoreDir)) {
	  scoreDir <- paste0(scoreDir,"/")
	}
	
	#if not output directory was defined
	if (is.na(outDir)) {
	  #derive one from the input
	  if (grepl("_mut_count/$",countDir)) {
	    outDir <- sub("_mut_count/$","_QC/",countDir)
	  } else {
	    outDir <- sub("/$","_QC/",countDir)
	  }
	} 
	#make sure it ends in "/"
	if (!grepl("/$",outDir)) {
	  outDir <- paste0(outDir,"/")
	}
	
	#make sure outdir exists
	dir.create(outDir,recursive=TRUE,showWarnings=FALSE)
	
	logInfo("Using count directory",countDir,
	        "score directory", scoreDir,
	        "and output directory",outDir
	)
	#create PDF tag
	pdftag <- with(params,sprintf("%s (%s): %s%s",project,template$geneName,runLabel,timeStamp))
	params$pdftagbase <- pdftag

	
	logInfo("Reading count data")
	marginalCountFile <- paste0(countDir,"/marginalCounts.csv")
	if (!file.exists(marginalCountFile)) {
	  stop("Invalid counts directory ",countDir,"! Must contain marginalCounts.csv!")
	}
	marginalCounts <- read.csv(marginalCountFile,comment.char="#")
	rownames(marginalCounts) <- marginalCounts$hgvsc

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
			scoreFile <- paste0(scoreDir,"/",sCond,"_t",tp,"_enrichment.csv")
			if (!file.exists(scoreFile)) {
				logWarn("No score file found! Skipping...")
				next
			}
			scores <- read.csv(scoreFile,comment.char="#")

			#ordering should match scores
			rownames(scores) <- scores$hgvsc
			# scores <- scores[marginalCounts$hgvsc,]
			marginalSubset <- marginalCounts[scores$hgvsc,]

			#Score distributions & syn/non medians
			scoreDistributions(scores,sCond,tp,outDir,params,srOverride)
			
			#plot logPhi bias vs marginal frequency
			logPhiBias(scores,params,sCond,tp,outDir)

			#filter progression graph
			filterProgression(scores,sCond,tp,params,outDir)
			
			#examine codon agreement for same amino acids
			codonAgreement(scores,sCond,tp,params,outDir,srOverride)
			
			#running mean of synonymous, stop and their difference
			synNonDelta(scores,sCond,tp,params,outDir,srOverride)

			#all of these analyses require more than one replicate
			# if (params$numReplicates[[sCond]] > 1) {
			if (!srOverride) {

				#run replicate correlation analysis
				replicateCorrelation(scores, marginalSubset, params, sCond, tp, outDir)

				#load error model file
				modelFile <- paste0(scoreDir,"/",sCond,"_t",tp,"_errorModel.csv")
				if (!file.exists(modelFile)) {
					logWarn("No error model file found. Skipping regularization QC.")
				} else {
					#Regularization analysis
					modelParams <- read.csv(modelFile,row.names=1)
					regularizationQC(scores,modelParams,params,sCond,tp,outDir)
				}
			
				#If scores could not be assigned due to synonymous-nonsense median failure
				#then we can't run an error profile analysis
				if (!all(is.na(scores$logPhi)) && !any(scores$logPhi.sd < 0,na.rm=TRUE)) {
					#Error profile
					errorProfile(scores,sCond,tp,outDir,params)
				} else {
					logWarn("Cannot plot error profiles: logPhi values are not available!")
				}
			}

		}
	}

}

#' Draw replicate correlation plots
#' 
#' @param scores the score table
#' @param params the global parameter object
#' @return NULL
filterProgression <- function(scores,sCond,tp,params,outDir) {

  regSubsets <- data.frame(
    set=c("all", params$regions$`Region Number`),
    start=c(2,params$regions$`Start AA`),
    end=c(params$template$proteinLength,params$regions$`End AA`),
    len=c(params$template$proteinLength-1,params$regions$`End AA`-params$regions$`Start AA`+1)
  )
  
  
	#get reachable AA changes
	#(this includes stop, but not synonymous)
	reachable <- reachableChanges(params)
	
	#prep plot
	outfile <- paste0(outDir,sCond,"_t",tp,"_filtering.pdf")
	pdf(outfile,11,8.5)
	tagger <- pdftagger(paste(params$pdftagbase,"; selection condition:",sCond),cpp=4)
	opar <- par(oma=c(2,2,2,2),mar=c(1,1,4,1)+.1,mfrow=c(2,2))
	
	for (ri in 1:nrow(regSubsets)) {
	  
	  #and extract reachable codon changes
	  reachableCCs <- do.call(c,with(reachable,mapply(function(w,p,ms){
	    if (p >= regSubsets$start[[ri]] && p <= regSubsets$end[[ri]]) {
	      paste0(w,p,ms)
	    } else NULL
	  },wtcodon,pos,strsplit(mutcodons,"\\|"))))
	  
	  reachableAAs <- reachable[with(reachable,pos >= regSubsets$start[[ri]] & pos <= regSubsets$end[[ri]]),]
	  
	  #make filtered subsets of the score table
	  scorePos <- as.integer(gsub("\\D","",scores$aaChange))
	  localScores <- scores[which(scorePos >= regSubsets$start[[ri]] & 
	                                scorePos <= regSubsets$end[[ri]]),]
	  filteredScores <- scores[which(is.na(scores$filter) & scorePos >= regSubsets$start[[ri]] & 
	                             scorePos <= regSubsets$end[[ri]]),]
	  hqScores <- filteredScores[filteredScores$bce.sd < params$scoring$sdThreshold,]
	  
  	#calculate filter census
  	census <- rbind(
  		possible = c(
  		  AllCCs = regSubsets$len[[ri]] * (4^3-1),
  			ReachCCs = length(reachableCCs),
  			AllAACs = regSubsets$len[[ri]] * (19+1),
  			ReachAACs = nrow(reachableAAs)
  		),
  		found = c(
  			AllCCs = nrow(localScores),
  			ReachCCs = length(intersect(localScores$codonChange, reachableCCs)),
  			AllAACs = length(unique(localScores$hgvsp[scores$type != "synonymous"])),
  			ReachAACs = length(intersect(unique(localScores$hgvsp),reachable$hgvsp))
  		),
  		filtered = c(
  			AllCCs = nrow(filteredScores),
  			ReachCCs = length(intersect(filteredScores$codonChange, reachableCCs)),
  			AllAACs = length(unique(filteredScores$hgvsp[filteredScores$type != "synonymous"])),
  			ReachAACs = length(intersect(unique(filteredScores$hgvsp),reachable$hgvsp))
  		),
  		hiQual = c(
  			AllCCs = nrow(hqScores),
  			ReachCCs = length(intersect(hqScores$codonChange, reachableCCs)),
  			AllAACs = length(unique(hqScores$hgvsp[hqScores$type != "synonymous"])),
  			ReachAACs = length(intersect(unique(hqScores$hgvsp),reachable$hgvsp))
  		)
  	)
  
  	#draw plot
  	widths <- census/max(census)
  	percentages <- apply(census,2,function(xs)xs/xs[[1]])*100
  	
  	plotCols <- c("steelblue2","steelblue3","gold2","gold3")
  	ylabels <- c("All possible","Detected","Passed filter","High Quality")
  	xlabels <- c("All","SNV-reachable","All","SNV-reachable")
  	toplabels <- c("Codon changes","AA changes")
  	title <- if (regSubsets$set[[ri]]=="all") "Entire Construct" else {
  	  paste0("Region #",regSubsets$set[[ri]])
  	}
  	cex <- 0.8
  	
  	plot(NA,type="n",xlim=c(-.5,4.5),ylim=c(1,5),xlab="",ylab="",axes=FALSE,main=title)
  	abline(h=1:4,col="gray",lty="dotted")
  	text(-0.5,4:1,ylabels,pos=4,cex=cex)
  	text(1:4,4.4,xlabels,cex=cex)
  	text(c(1.5,3.5),4.8,toplabels,cex=cex)
  	invisible(lapply(1:4, function(cati) {
  		polygon(
  			c(cati-widths[,cati]/2,rev(cati+widths[,cati]/2)),
  			c(4:1,1:4),col=plotCols[[cati]], border=NA
  		)
  		text(cati,4:1,sprintf("%d (%.02f%%)",census[,cati],percentages[,cati]),cex=cex)
  	}))
  	tagger$cycle()
	
	}
	par(opar)
	invisible(dev.off())

}


#' Helper function to find the quad of conditions and time points that belong 
#' with the given selection condition
#'
#' @param params the parameter sheet object
#' @param sCond the current selection condition
#' @param tp the current time point of the selection condition
#'
#' @return a list containing the condition quad ('condQuad'), 
#' replicate matrix ('repMatrix') and timepoint quad ('tpQuad')
#' @export
#'
findQuads <- function(params,sCond,tp) {
  
  # nrep <- params$numReplicates[[sCond]]
  #this is to make up for missing WT conditions
  null2na <- function(x) if (length(x) == 0) NA else x
  
  #pull up matching nonselect and WT controls
  nCond <- getNonselectFor(sCond,params)
  condQuad <- c(
    select=sCond,
    nonselect=nCond,
    selWT=null2na(getWTControlFor(sCond,params)),
    nonWT=null2na(getWTControlFor(nCond,params))
  )
  sRep <- params$numReplicates[[sCond]]
  #if a condition is missing entirely, we give it one "pseudoreplicate with 0 scores"
  repQuad <- sapply(condQuad,function(con) if(is.na(con)) 1 else params$numReplicates[[con]])
  if (!all(repQuad == sRep)) {
    logWarn(paste(
      "Number of replicates in conditions is not balanced or WT is missing!",
      " => Correlation plot may be distorted due to recycled replicates!!",
      sep="\n"
    ))
  }
  #Workaround: Since R doesn't like variable names starting with numerals, 
  # we need to adjust any of those cases
  if (any(grepl("^\\d",condQuad))) {
    culprits <- which(grepl("^\\d",condQuad))
    #we add the prefix "X" to the name, just like the table parser does.
    condQuad[culprits] <- paste0("X",condQuad[culprits])
  }
  
  #deal with time points
  tpQuad <- sapply(condQuad, function(con) {
    if (is.na(con)) {
      NA
    } else if (tp %in% getTimepointsFor(con,params)) {
      tp
    } else {
      #hopefully there's only one other as you are only allowed to either have 
      #the same timepoints everywhere or just one
      getTimepointsFor(con,params)[[1]]
    }
  })
  
  #replicate column name matrix
  repMatrix <- do.call(rbind,lapply(names(condQuad),function(con) {
    sapply(1:repQuad[[con]], 
           function(repi) {
             if (is.na(condQuad[[con]])) {
               NA #again, compensating for potentially missing WT condition
             } else {
               sprintf("%s.t%s.rep%d.frequency",condQuad[[con]],tpQuad[[con]],repi)
             }
           }
    )
  }))
  rownames(repMatrix) <- names(condQuad)
  
  return(list(condQuad=condQuad,repMatrix=repMatrix,tpQuad=tpQuad))
}

#' Plot the logPhi bias vs the nonselect marginal frequency
#'
#' @param scores the score table
#' @param params the parameter sheet object
#' @param sCond the active selection condition
#' @param tp the active time point
#' @param outDir the output directory
#'
#' @return nothing. writes a pdf file to the output dir
#' @export
logPhiBias <- function(scores,params,sCond,tp,outDir) {
  
  scores$position <- as.integer(extract.groups(scores$codonChange,"(\\d+)"))
  scores$region <- params$pos2reg(scores$position)
  # scores$tile <- params$pos2tile(scores$position)
  # 
  # nonsense <- scores[scores$type=="nonsense",]
  # nonsenseF <- scores[scores$type=="nonsense" & is.na(scores$filter),]
  # zns <- with(nonsenseF, lm(logPhi~log10(nonselect.mean)) )
  # trns <- with(nonsenseF, yogitools::runningFunction(log10(nonselect.mean),logPhi,nbins=30))
  # 
  # synonymous <- scores[scores$type=="synonymous",]
  # synonymousF <- scores[scores$type=="synonymous" & is.na(scores$filter),]
  # zsyn <- with(synonymousF, lm(logPhi~log10(nonselect.mean)) )
  # trsyn <- with(synonymousF, yogitools::runningFunction(log10(nonselect.mean),logPhi,nbins=30))
  # 
  outfile <- paste0(outDir,sCond,"_t",tp,"_logPhiBias.pdf")
  pdf(outfile,8.5,11)
  tagger <- pdftagger(paste(params$pdftagbase,"; selection condition:",sCond),cpp=6)
  opar <- par(oma=c(15,2,2,2),mfrow=c(3,2))
  
  for (currRegion in params$regions[,1]) {
    
    nonsenseF <- scores[with(scores,type=="nonsense" & is.na(filter) & logPhi < 0 & region==currRegion),]
    synonymousF <- scores[with(scores,type=="synonymous" & is.na(filter) & region==currRegion),]
    
    zns <- with(nonsenseF, lm(logPhi~log10(nonselect.mean)) )
    zsyn <- with(synonymousF, lm(logPhi~log10(nonselect.mean)) )
    
    trns <- with(nonsenseF, yogitools::runningFunction(log10(nonselect.mean),logPhi,nbins=20))
    trsyn <- with(synonymousF, yogitools::runningFunction(log10(nonselect.mean),logPhi,nbins=20))
    with(rbind(nonsenseF,synonymousF), plot(
      log10(nonselect.mean),logPhi, type="n", main=paste("Region",currRegion)
    ))
    with(nonsenseF, {
      points(log10(nonselect.mean),logPhi,col="firebrick3",pch=20)
      yogitools::errorBars(log10(nonselect.mean),logPhi,logPhi.sd,col="firebrick3")
    })
    abline(zns,col="firebrick2",lty="dashed")
    # abline(a=thetaNon[[1]],b=thetaNon[[2]],col="firebrick2")
    lines(trns,col="firebrick2",lwd=2)
    
    with(synonymousF, {
      points(log10(nonselect.mean),logPhi,col="chartreuse3",pch=20)
      yogitools::errorBars(log10(nonselect.mean),logPhi,logPhi.sd,col="chartreuse3")
    })
    abline(zsyn,col="chartreuse2",lty="dashed")
    # abline(a=thetaSyn[[1]],b=thetaSyn[[2]],col="chartreuse2")
    lines(trsyn,col="chartreuse2",lwd=2)
    
    
    # plotcols <- sapply(c("gray","firebrick3"),yogitools::colAlpha,0.5)
    # with(nonsense, plot(
    #   log10(nonselect.mean),logPhi,
    #   col=plotcols[1+is.na(filter)],
    #   pch=20
    # ))
    # # abline(h=0:1,lty="dashed",col=2:3)
    # abline(zns,col="firebrick2")
    # lines(trns,col="firebrick2",lwd=2)
    # plotcols <- sapply(c("gray","chartreuse3"),yogitools::colAlpha,0.5)
    # with(synonymous, points(
    #   log10(nonselect.mean),logPhi,
    #   col=plotcols[1+is.na(filter)],
    #   pch=20
    # ))
    # # abline(h=0:1,lty="dashed",col=2:3)
    # abline(zsyn,col="chartreuse2")
    # lines(trsyn,col="chartreuse2",lwd=2)
    # plotcols <- sapply(c("chartreuse3","firebrick3","gray"),yogitools::colAlpha,0.5)
    # legend("bottomleft",c("synonymous","nonsense","filtered out"),col=plotcols,pch=20)
    tagger$cycle()
  }
  
  par(opar)
  invisible(dev.off())
  
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
  
  sRep <- params$numReplicates[[sCond]]
  quads <- findQuads(params,sCond,tp)
  repMatrix <- quads$repMatrix
  
	#check that labels match between tables
	if (!all(scores$hgvsp == marginalCounts$hgvsp)) {
	  stop("scores and marginal count files mismatch.")
	}
	
	#extract replicate values for this condition
	repValues <- lapply(1:sRep, function(repi) {
	  selRaw <- marginalCounts[,repMatrix["select",repi]]
	  selWTraw <- if (!is.na(repMatrix["selWT",repi])) {
	    marginalCounts[,repMatrix["selWT",repi]]
	  } else 0
		selFreq <-floor0(selRaw - selWTraw)
		
		nonRaw <- marginalCounts[,repMatrix["nonselect",repi]]
		nonWTraw <- if (!is.na(repMatrix["nonWT",repi])) {
		  marginalCounts[,repMatrix["nonWT",repi]]
		} else 0
		nonFreq <-  floor0(nonRaw - nonWTraw)
		
		smallestSelect <- unique(sort(na.omit(selFreq)))[[2]]
		pseudoCount <- 10^floor(log10(smallestSelect))
		
		logphi <- log10((selFreq+pseudoCount) / nonFreq)
		data.frame(select=selFreq,nonselect=nonFreq,logphi=logphi)
	})
	names(repValues) <- as.character(1:sRep)
	repValues <- do.call(cbind,repValues)
	#apply filter from scoring function
	repValues <- repValues[is.na(scores$filter),]


	if (sRep == 2) {
		outfile <- paste0(outDir,sCond,"_t",tp,"_replicates.pdf")
		pdf(outfile,11,8.5)
		tagger <- pdftagger(paste(params$pdftagbase,"; selection condition:",sCond),cpp=2)
		layout(cbind(1,2))
		opar <- par(oma=c(2,2,2,2),mar=c(15,4,4,1)+.1)
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
		tagger$cycle()
		plot(
			repValues[,"1.logphi"],repValues[,"2.logphi"],
			xlab="log(phi) Replicate 1", ylab="log(phi) Replicate 2",
			main=sprintf(
				"select / nonselect log-ratio R = %.03f",
				cor(fin(repValues[,sprintf("%d.logphi",1:2)]))[1,2]
			),pch="."
		)
		par(opar)
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
		# imgSize <- max(4,sRep)

		#TODO: This needs to be tested!
		outfile <- paste0(outDir,sCond,"_t",tp,"_ns_replicates.pdf")
		# pdf(outfile,imgSize,imgSize)
		pdf(outfile,8.5,11)
		# tagger <- pdftagger(paste(params$pdftagbase,"; selection condition:",sCond),cpp=2)
		pairs(
			repValues[,sprintf("%d.nonselect",1:sRep)],
			lower.panel=panel.cor,pch=".",labels=labels,
			main="non-select frequencies"
		)
		op <- par(oma=c(2,2,2,2))
		mtext(paste(params$pdftagbase,"; selection condition:",sCond),side=1,outer=TRUE,line=0,cex=0.5)
		par(op)
		invisible(dev.off())

		outfile <- paste0(outDir,sCond,"_t",tp,"_phi_replicates.pdf")
		# pdf(outfile,imgSize,imgSize)
		pdf(outfile,8.5,11)
		pairs(
			repValues[,sprintf("%d.nonselect",1:sRep)],
			lower.panel=panel.cor,pch=".",labels=labels,
			main="select / nonselect log-ratios"
		)
		op <- par(oma=c(2,2,2,2))
		mtext(paste(params$pdftagbase,"; selection condition:",sCond),side=1,outer=TRUE,line=0,cex=0.5)
		par(op)
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
	# tileStarts <- params$tiles[,"Start AA"]
	positions <- as.integer(extract.groups(scores$codonChange,"(\\d+)")[,1])
	# tiles <- sapply(positions,function(pos) max(which(tileStarts <= pos)))
	tiles <- params$pos2tile(positions)
	
	outfile <- paste0(outDir,sCond,"_t",tp,"_errorModel.pdf")
	pdf(outfile,8.5,11)
	tagger <- pdftagger(paste(params$pdftagbase,"; selection condition:",sCond),cpp=6)
	opar <- par(mfrow=c(3,2),oma=c(2,2,2,2))
	for (tile in params$tiles[,"Tile Number"]) {
		# model.fit <- tryCatch({
		# 	fit.cv.model(scores[which(tiles==tile),])
		# },error=function(e) {
		# 	NULL
		# })
		if (!(tile %in% tiles)) {
			plot.new()
			rect(0,0,1,1,col="gray80",border="gray30",lty="dotted")
			text(0.5,0.5,"no data")
			mtext(paste0("Tile ",tile),side=3)
			tagger$cycle()
			next
		}

		with(scores[which(tiles==tile),],{

			theta <- modelParams[as.character(tile),paste0("nonselect.",c("static","additive","multiplicative"))]
			cv.model <- function(count) {
				10^sapply(log10(1/sqrt(count)),function(x) max(theta[[1]], theta[[2]] + theta[[3]]*x))
			}

			plot(nonselect.count,nonselect.cv,log="xy",main=paste("Tile",tile,"non-select"))
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
			tagger$cycle()

			theta <- modelParams[as.character(tile),paste0("select.",c("static","additive","multiplicative"))]
			cv.model <- function(count) {
				10^sapply(log10(1/sqrt(count)),function(x) max(theta[[1]], theta[[2]] + theta[[3]]*x))
			}

			plot(select.count,select.cv,log="xy",main=paste("Tile",tile,"select"))
			runningMean <- runningFunction(
				select.count,select.cv,nbins=20,logScale=TRUE
			)
			sSamples <- seq(1,max(select.count,na.rm=TRUE),length.out=100)
			lines(sSamples,1/sqrt(sSamples),col="chartreuse3",lty="dashed",lwd=2)
			lines(runningMean,col="firebrick3",lwd=2)
			
			lines(sSamples,cv.model(sSamples),col="blue",lwd=2)
			mtext(sprintf(
				"stat.=%.02f; add.=%.02f; mult.=%.02f",
				theta[[1]],theta[[2]],theta[[3]]
			))
			tagger$cycle()

			# nonselect.sd.poisson <- 1/sqrt(nonselect.count)*nonselect.mean
			# nonselect.sd.poisson[which(nonselect.mean==0)] <- 0
			# plot(nonselect.sd.poisson, nonselect.sd,log="xy")
			# runningMean <- runningFunction(
			# 	nonselect.sd.poisson, nonselect.sd,nbins=20,logScale=TRUE
			# )
			# lines(runningMean,col="firebrick3",lwd=2)
			# abline(0,1,col="chartreuse3",lty="dashed",lwd=2)
			# depth <- mean(nonselect.count/nonselect.mean,na.rm=TRUE)
			# lines(
			# 	sqrt(nsSamples)/depth,
			# 	cv.model(nsSamples)*(nsSamples/depth),
			# 	col="blue",lwd=2
			# )
		})
	}
	par(opar)
	invisible(dev.off())
}


#' Draw score distribution plots
#' 
#' @param scores the score table
#' @param sCond the condition for which to draw the plot
#' @param tp the time point
#' @param outDir the output directory
#' @return NULL
scoreDistributions <- function(scores,sCond,tp,outDir,params,srOverride) {
  
  filteredScores <- scores[is.na(scores$filter),]
  if (!srOverride) {
    resErr <- residualError(filteredScores$bce,filteredScores$bce.sd)
    nerr <- resErr - min(resErr,na.rm=TRUE) + min(abs(resErr),na.rm=TRUE)
  }
  
	#collapse by amino acid consequence and associate with regions
	aaScores <- as.df(with(filteredScores,tapply(1:length(hgvsp),hgvsp, function(is) {
		if (!srOverride && !any(is.na(bce.sd[is]))) {
			joint <- join.datapoints(
			  ms=bce[is],
			  sds=bce.sd[is],
				dfs=rep(params$numReplicates[[sCond]],length(is)),
				ws=(1/nerr[is])/sum(1/nerr[is])
			)
		} else {
			#this is the case if srOverride is turned on and only
			#one replicate was available
			joint <- c(
				mj=mean(bce[is],na.rm=TRUE),sj=NA,
				dfj=params$numReplicates[[sCond]]*length(is)
			)
		}
	  #extract AA position and derive region assignment
		p <- unique(as.integer(extract.groups(aaChange[is],"(\\d+)")[,1]))
		mutregion <- params$pos2reg(p)
		list(
			hgvsp=unique(hgvsp[is]),
			bce=joint[["mj"]],
			sd=joint[["sj"]],
			df=joint[["dfj"]],
			se=joint[["sj"]]/sqrt(joint[["dfj"]]),
			pos=p,
			region=mutregion
		)
	})))
	
	sdCutoff <- params$scoring$sdThreshold

	#subdivide data into regions
	# mutpos <- as.integer(extract.groups(scores$aaChange,"(\\d+)")[,1])
	# mutregion <- with(as.data.frame(params$regions),sapply(mutpos,function(p) {
	# 	which(p >= `Start AA` & p <= `End AA`)
	# }))

	outfile <- paste0(outDir,sCond,"_t",tp,"_logPhiDistribution.pdf")
	pdf(outfile,11,8.5)
	opar <- par(oma=c(2,2,2,2))
	tagger <- pdftagger(paste(params$pdftagbase,"; selection condition:",sCond),cpp=2)
	layout(rbind(1,2,3,4),heights=c(1.2,1,1.2,1))
	invisible(tapply(1:nrow(aaScores),aaScores$region, function(is) {
		reg <- unique(aaScores$region[is])
		drawDistributions(aaScores[is,],Inf,reg)
		tagger$cycle()
		if (!srOverride && !any(is.na(aaScores[is,"se"]))) {
			drawDistributions(aaScores[is,],sdCutoff,reg)
			tagger$cycle()
		} else {
		  plot.new()
		  rect(0,0,1,1,col="gray80",border="gray30",lty="dotted")
		  text(0.5,0.5,"unable to perform\nquality filtering")
		  plot.new()
		  tagger$cycle()
		}
		return(NULL)
	}))
	drawDistributions(aaScores,Inf,"all")
	tagger$cycle()
	if (!any(is.na(aaScores[,"se"]))) {
		drawDistributions(aaScores,sdCutoff,"all")
		tagger$cycle()
	}
	par(opar)
	invisible(dev.off())
}

#' Delegation function to draw score distribution plots with a given filter setting
#' 
#' @param aaScores the score table
#' @param seCutoff the stderr cutoff to apply
#' @return NULL
drawDistributions <- function(aaScores,seCutoff=Inf,reg=NA) {

  if (seCutoff < Inf) {
    #check that the cutoff leaves at least 10 nonsense variants to work with
    #otherwise, make it more permissive
  	numNsSurvive <- with(aaScores,sum(se < seCutoff & grepl("Ter$",hgvsp),na.rm=TRUE))
  	if (numNsSurvive < 10) {
  		nonsenseSDs <- with(aaScores,se[grepl("Ter$",hgvsp)])
  		r10Threshold <- sort(nonsenseSDs)[[min(10,length(nonsenseSDs))]]
  		logWarn(sprintf("sdThreshold %.03f is too restrictive! Using se < %.03f instead.",seCutoff,r10Threshold))
  		seCutoff <- r10Threshold
  	}
  }

	#extract filtered scores
	if (!all(is.na(aaScores$se))) {
		synScores <- with(aaScores,bce[grepl("=$",hgvsp) & se < seCutoff ])
		stopScores <- with(aaScores,bce[grepl("Ter$",hgvsp) & se < seCutoff])
		misScores <- with(aaScores,bce[!grepl("Ter$|=$",hgvsp) & se < seCutoff])
	} else {
		synScores <- with(aaScores,bce[grepl("=$",hgvsp)])
		stopScores <- with(aaScores,bce[grepl("Ter$",hgvsp)])
		misScores <- with(aaScores,bce[!grepl("Ter$|=$",hgvsp)])
	}
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
		main=if (is.infinite(seCutoff)) {
			paste("Region",reg,"; Unfiltered")
		} else {
			bquote("Region"~.(reg)~";"~sigma < .(seCutoff))
		}
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
	#draw medians and percentile
	abline(v=uCoord(median(synScores)),col="darkolivegreen4",lwd=2,lty="dotted")
	abline(v=uCoord(median(stopScores)),col="firebrick3",lwd=2,lty="dotted")
	abline(v=uCoord(quantile(misScores,0.1)),col="gray30",lwd=2,lty="dotted")
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
	text(
		uCoord(quantile(misScores,0.1)),
		-max(misHist$density)/2,
		sprintf("10th percentile\n%.03f",quantile(misScores,0.1)),
		col="gray30"
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
errorProfile <- function(scores,sCond,tp,outDir,params) {
  
  plotEP <- function(scores,varname,xlab) {
    sdname <- paste0(varname,".sd")
    sdRange <- range(log10(scores[is.na(scores$filter),sdname]),finite=TRUE)
    lphiRange <- range(scores[is.na(scores$filter),varname],finite=TRUE)
    runningMedian <- with(scores[is.na(scores$filter),],yogitools::runningFunction(
      get(varname),log10(get(sdname)),fun=median,nbins=20
    ))
    runningMedian[,2] <- 10^runningMedian[,2]
    expectation <- with(scores[is.na(scores$filter),],expectedError(get(varname),get(sdname)))
    par(mar=c(5,4,0,0)+.1)
    with(scores[is.na(scores$filter),],plot(
      get(varname),get(sdname),log="y",pch=".",
      xlab=xlab,ylab=expression(sigma)
    ))
    lines(runningMedian,col="steelblue3",lwd=2)
    curve(expectation,add=TRUE,col="firebrick3",lwd=2)
    par(mar=c(0,4,1,0)+.1) 
    breaks <- seq(lphiRange[[1]],lphiRange[[2]],length.out=50)
    lphiHist <- with(scores[is.na(scores$filter),], hist(get(varname),breaks=breaks,plot=FALSE))
    barplot(lphiHist$density,border=NA,ylab="density",space=0)
    par(mar=c(5,0,0,0)+.1) 
    breaks <- seq(sdRange[[1]],sdRange[[2]],length.out=50)
    sdHist <- with(scores[is.na(scores$filter),], hist(log10(get(sdname)),breaks=breaks,plot=FALSE))
    barplot(sdHist$density,border=NA,horiz=TRUE,space=0,xlab="density")
    tagger$cycle()
  }
  
	outfile <- paste0(outDir,sCond,"_t",tp,"_errorProfile.pdf")
	pdf(outfile,8.5,11)
	layout(rbind(c(2,0),c(1,3)),widths=c(0.8,0.2),heights=c(0.2,0.8))
	tagger <- pdftagger(paste(params$pdftagbase,"; selection condition:",sCond),cpp=4)
	opar <- par(oma=c(24,6,2,6)) 
	plotEP(scores,"logPhi",expression(log[10](phi)))
	plotEP(scores,"bce","bias-corrected enrichment")
	par(opar)

	invisible(dev.off())
	
}



#' Draw plots to analyze agreement between equivalent codons
#'
#' @param scores the score table
#' @param sCond currently active selection condition
#' @param tp currently active time point
#' @param params parameter sheet object
#' @param outdir output directory
#' @return nothing. writes plots to output directory
codonAgreement <- function(scores,sCond,tp,params,outDir,srOverride) {
  codonGroups <- tapply(1:nrow(scores),scores$hgvsp,c)
  nCodons <- sapply(codonGroups,length)
  
  #create pairwise combinations
  pairings <- t(do.call(cbind,lapply(codonGroups[nCodons > 1],combn,2)))
  pairScores <- as.df(apply(pairings,1,function(idxs) {
    filter <- any(!is.na(scores[idxs,"filter"]))
    lphi <- scores[idxs,"logPhi"]
    lpsd <- scores[idxs,"logPhi.sd"]
    c(list(filter=filter),c(logPhi=lphi,sd=lpsd))
  }))
  
  alphaCol <- colAlpha("black",0.2)
  
  outfile <- paste0(outDir,sCond,"_t",tp,"_codonCorr.pdf")
  tagger <- pdftagger(paste(params$pdftagbase,"; selection condition:",sCond),cpp=1)
  pdf(outfile,8.5,11)
  
  layout(cbind(1:3),heights=c(1,2,2))
  op <- par(mar=c(5,4,.1,1),oma=c(2,15,2,15))
  barplot(
    table(nCodons),
    xlab="#equivalent codons per AA",
    ylab="Frequency", col="steelblue3",border=NA
  )
  
  with(pairScores[!pairScores$filter,],{
    
    #plot codon scores against each other
    par(mar=c(5,4,2,1))
    plot(
      logPhi1,logPhi2,type="n",
      xlab=expression("codon 1"~log(phi)),
      ylab=expression("codon 2"~log(phi))
    )
    yogitools::errorBars(
      logPhi1,logPhi2,sd2,col=alphaCol
    )
    yogitools::errorBars(
      logPhi2,logPhi1,sd1,vertical=FALSE,col=alphaCol
    )
    abline(0,1,col="gray",lty="dashed")
    mtext(sprintf("PCC=%.02f",cor(logPhi1,logPhi2)))
   
    #plot score difference against max(sd)
    par(mar=c(5,4,2,1))
    if (!srOverride) {
      plot(
        abs(logPhi1-logPhi2),mapply(max,sd1,sd2),log="y",
        pch=20,col=alphaCol,
        xlab=expression(Delta~log(phi)),ylab=expression(max(sigma))
      )
    } else {
      plot.new()
      rect(0,0,1,1,col="gray80",border="gray30",lty="dotted")
      text(0.5,0.5,"no replicate data")
    }
  })
  tagger$cycle()
  par(op)
  invisible(dev.off())
  
}


#' Draws a plot showing running means of synonymous and nonsense variant scores
#' across the length of the protein, as well as difference track between the two.
#'
#' @param scores data.frame of the scores
#' @param sCond the currently active selective condition
#' @param tp the currently active time point
#' @param params the parameter sheet object
#' @param outDir the output directory
#' @return nothing, writes plot to output directory
synNonDelta <- function(scores,sCond,tp,params,outDir,srOverride){
  scores$pos <- as.integer(gsub("\\D+","",scores$aaChange))
  #positional weighted averages of synonymous variants
  #FIXME: Use residual error as weights
  syns <- as.df(with(scores[scores$type=="synonymous" & is.na(scores$filter),],tapply(1:length(pos),pos,function(idxs){
    sds <- if (!srOverride) logPhi.sd[idxs] else rep(1,length(idxs))
    joint <- join.datapoints(logPhi[idxs],sds,rep(2,length(idxs)))
    p <- unique(pos[idxs])
    c(pos=p,joint)
  })))
  #positional weighted averages of nonsense (stop) variants
  stops <- as.df(with(scores[scores$type=="nonsense" & is.na(scores$filter),],tapply(1:length(pos),pos,function(idxs){
    sds <- if (!srOverride) logPhi.sd[idxs] else rep(1,length(idxs))
    joint <- join.datapoints(logPhi[idxs],sds,rep(2,length(idxs)))
    p <- unique(pos[idxs])
    c(pos=p,joint)
  })))
  
  allpos <- as.character(1:params$template$proteinLength)
  #This is a hack to get around a bug in `[.dataframe` which pulls !approximate! matches!!
  #see: https://stackoverflow.com/questions/34233235/r-returning-partial-matching-of-row-names
  pull <- function(tbl,idxs) as.df(lapply(idxs, function(idx) {
    if (idx %in% rownames(tbl)) {
      tbl[idx,2:4]
    } else setNames(rep(NA,3),colnames(tbl)[2:4])
  }))
  joint <- cbind(pos=as.integer(allpos),syn=pull(syns,allpos),non=pull(stops,allpos))
  runningWeights <- as.df(lapply(joint$pos,function(p) {
    idxs <- which(abs(joint$pos-p) < 5)
    synsubset <- na.omit(joint[idxs,c("syn.mj","syn.sj","syn.dfj")])
    synAv <- join.datapoints(synsubset[,1],synsubset[,2],synsubset[,3])
    nonsubset <- na.omit(joint[idxs,c("non.mj","non.sj","non.dfj")])
    nonAv <- join.datapoints(nonsubset[,1],nonsubset[,2],nonsubset[,3])
    c(pos=p,synAv=synAv[["mj"]],synSD=synAv[["sj"]],nonAv=nonAv[["mj"]],nonSD=nonAv[["sj"]])
  }))
  runningWeights$delta <- with(runningWeights,synAv - nonAv)
  runningWeights$deltaSD <- with(runningWeights,sqrt(synSD^2 + nonSD^2))
  
  outfile <- paste0(outDir,sCond,"_t",tp,"_synNonDiff.pdf")
  tagger <- pdftagger(paste(params$pdftagbase,"; selection condition:",sCond),cpp=1)
  pdf(outfile,11,8.5)
  layout(rbind(1,2),heights=c(0.3,1))
  op <- par(mar=c(0,4,1,1),oma=c(12,2,12,2))
  with(runningWeights,{
    plot(NA,type="n",xlim=range(pos),ylim=c(0,1),xlab="",ylab="",axes=FALSE)
    rect(params$regions[,"Start AA"],0,params$regions[,"End AA"],0.49,col="gray80",border=NA)
    text(rowMeans(params$regions[,c("Start AA","End AA")]),0.25,params$regions[,"Region Number"])
    rect(params$tiles[,"Start AA"],0.51,params$tiles[,"End AA"],1,col="gray90",border=NA)
    text(rowMeans(params$tiles[,c("Start AA","End AA")]),0.75,params$tiles[,"Tile Number"])
    par(mar=c(5,4,0,1))
    plot(NA,type="n",xlim=range(pos),ylim=range(c(synAv,nonAv)),
         xlab="AA position",ylab=expression("running average"~log(phi))
    )
    polygon(c(pos,rev(pos)),c(synAv+synSD/2,rev(synAv-synSD/2)),col=colAlpha("chartreuse3",0.2),border=NA)
    lines(pos,synAv,col="chartreuse3")
    polygon(c(pos,rev(pos)),c(nonAv+nonSD/2,rev(nonAv-nonSD/2)),col=colAlpha("firebrick3",0.2),border=NA)
    lines(pos,nonAv,col="firebrick3")
    abline(v=params$tiles[-1,"Start AA"],col="gray90",lty="dashed")
    abline(v=params$regions[-1,"Start AA"],col="gray80",lty="dashed")
  })
  tagger$cycle()
  par(op)
  invisible(dev.off())

}


