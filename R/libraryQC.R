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
#' @param inDir input directory, defaults to subdirectory with latest timestamp ending in _mut_count.
#' @param outDir output directory, defaults to name of input directory with _QC tag attached.
#' @param paramFile input parameter file. defaults to <dataDir>/parameters.json
#' @param mc.cores the maximum number of processes to run in parallel for multi-core
#'   processing. Warning: This also multiplies the amount of RAM used!
#' @param srOverride override flag to allow for single-replicates
#' @param wmThreshold 'well-measuredness' threshold. The marginal frequency threshold 
#'   at which variants are considered well-measured. Defaults to 5e-5, which roughly 
#'   corresponds to the outdated definition in the legacy pipeline.
#' @param allCondOverride override flag to run on ALL conditions instead of just the nonselect
#' 
#' @return NULL. Results are written to file.
#' @export
libraryQC <- function(dataDir,inDir=NA,outDir=NA,paramFile=paste0(dataDir,"parameters.json"),
                      mc.cores=6,srOverride=FALSE, wmThreshold=5e-5,allCondOverride=FALSE,depthFilterOverride=FALSE) {
  
  op <- options(stringsAsFactors=FALSE)
  
  library(yogitools)
  library(hgvsParseR)
  library(pbmcapply)
  
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
  if (!is.na(outDir)){
    if (!grepl("/$",outDir)) {
      outDir <- paste0(outDir,"/")
    }
  }
  
  logInfo("Reading parameters from",normalizePath(paramFile))
  params <- withCallingHandlers(
    parseParameters(paramFile,srOverride=srOverride),
    warning=function(w)logWarn(conditionMessage(w))
  )
  
  #find counts folder
  if (is.na(inDir)) {
    latest <- latestSubDir(parentDir=dataDir,pattern="_mut_call$|mut_count$")
    inDir <- latest[["dir"]]
    timeStamp <- latest[["timeStamp"]]
    runLabel <- latest[["label"]]
  } else { #if custom input dir was provided
    #make sure it exists
    if (!dir.exists(inDir)) {
      stop("Input folder ",inDir," does not exist!")
    }
    #try to extract a timestamp and label
    lt <- extract.groups(inDir,"([^/]+_)?(\\d{4}-\\d{2}-\\d{2}-\\d{2}-\\d{2}-\\d{2})")[1,]
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
  if (!grepl("/$",inDir)) {
    inDir <- paste0(inDir,"/")
  }
  
  #if not output directory was defined
  if (is.na(outDir)) {
    #derive one from the input
    if (grepl("_mut_count/$",inDir)) {
      outDir <- sub("_mut_count/$","_QC/",inDir)
    } else {
      outDir <- sub("/$","_QC/",inDir)
    }
  } 
  #make sure it ends in "/"
  if (!grepl("/$",outDir)) {
    outDir <- paste0(outDir,"/")
  }
  
  #make sure outdir exists
  dir.create(outDir,recursive=TRUE,showWarnings=FALSE)
  
  logInfo("Using input directory",inDir,"and output directory",outDir)
  
  #create PDF tag
  tsmVersion <- sub(".9000$","",as.character(packageVersion("tileseqMave")))
  pdftag <- with(params,sprintf("tileseqPro v%s|%s (%s): %s%s",
                                tsmVersion,project,template$geneName,
                                runLabel,timeStamp
  ))
  
  #identify nonselect conditions
  nsConditions <- getNonselects(params)
  #apply override if requested
  if (allCondOverride) {
    nsConditions <- params$conditions$names
  }
  
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
  
  allCountFile <- paste0(inDir,"/allCounts.csv")
  marginalCountFile <- paste0(inDir,"/marginalCounts.csv")
  depthTableFile <-  paste0(inDir,"/sampleDepths.csv")
  covTableFile <-  paste0(inDir,"/positionalDepths.csv")
  
  if (!all(file.exists(allCountFile, marginalCountFile, depthTableFile))) {
    stop("Invalid input directory ",inDir," ! Must contain allCounts.csv, marginalCounts.csv and sampleDepths.csv !")
  }
  
  allCounts <- read.csv(allCountFile,comment.char="#")
  marginalCounts <- read.csv(marginalCountFile,comment.char="#")
  depthTable <- read.csv(depthTableFile)
  
  # COUNT ATTRIBUTES PLOT ------------------------------------------------------
  logInfo("Plotting sequencing anomalies")
  countAttributes(depthTable,inDir,outDir,pdftag)
  
  # PLOT PHIX ERROR RATES ------------------------------------------------------
  logInfo("Plotting sequencing error rates")
  plotSeqErrorRates(inDir,outDir,pdftag)
  
  if (file.exists(covTableFile)) {
    #read positional depth table and draw a diagnosis plot for it
    covTable <- read.csv(covTableFile,row.names=1)
    colnames(covTable) <- sub("^X","",colnames(covTable))
    
    logInfo("Plotting positional depths")
    drawPositionalDepth(covTable,outDir,pdftag)
    
    #adjust depth table with mean "effective" depth
    tileDepths <- adjustTableDepth(covTable, params)
    dimnames(tileDepths) <- list(rownames(covTable),params$tiles[,1])
  }
  
  #TODO: Maybe at a later point we will want to report on silent mutations too
  #but for now we remove them
  allCounts <- allCounts[which(allCounts$aaChanges != "silent"),]
  marginalCounts <- marginalCounts[which(marginalCounts$aaChange != "silent"),]
  
  logInfo("Interpreting variant descriptors. This may take some time...")
  
  allSplitChanges <- pbmclapply(
    strsplit(allCounts$codonChanges,"\\|"), 
    splitCodonChanges, mc.cores=mc.cores
  )
  marginalSplitChanges <- splitCodonChanges(marginalCounts$codonChange)
  
  #Infer tile assignments
  logInfo("Calculating tile assignments. This may take some time...")
  # tileStarts <- params$tiles[,"Start AA"]

  allTiles <- pbmclapply(allSplitChanges,inferTiles,params=params,mc.cores=mc.cores)   
  marginalTiles <- inferTiles(params, marginalSplitChanges)
  
  # CALCULATE SEVERITY OF DROPS IN SEQUENCING DEPTH ----------------------------
  if (file.exists(covTableFile)) {
    logInfo("Calculating drops in sequencing depth for each variant...")
    
    allDepthDrops <- calcDepthDrops(allCounts,allTiles,depthTable,mc.cores)
    margDepthDrops <- calcDepthDrops(marginalCounts,marginalTiles,depthTable,mc.cores)
    
    plotDepthDrops(params, outDir, pdftag, allDepthDrops)
  }
  
  #ITERATE OVER (POSSIBLY MULTIPLE) NONSELECT CONDITIONS AND ANALYZE SEPARATELY
  for (nsCond in nsConditions) {
    
    for (tp in params$timepoints$`Time point name`) {
      
      if (!(sprintf("%s.t%s.rep1.frequency",nsCond,tp) %in% colnames(marginalCounts))) {
        #if this time point does not exist for nonselect, then skip 
        logInfo("Skipping unused time point",tp,"for condition",nsCond)
        next
      }
      
      logInfo("Processing",nsCond,"; time point",tp)
      
      info = list(nsCond=nsCond, tp=tp, outDir=outDir, pdftag=pdftag)
      
      # FILTER OUT VARIANTS AT INSUFFICIENT DEPTH ------------------------------
      
      if (file.exists(covTableFile) && !depthFilterOverride) {
        logInfo("Checking for underpowered variants")
        
        #filter all variants
        filterVariants(nsCond, params, allDepthDrops, allCounts, 
                       allSplitChanges, allTiles, mc.cores)
        #filter marginal variants
        filterVariants(nsCond, params, margDepthDrops, marginalCounts, 
                       marginalSplitChanges, marginalTiles, mc.cores, 
                       allVariants = FALSE)
      }
      
      # RAW REPLICATE CORRELATIONS PER TILE ------------------------------------
      
      nsReps <- sprintf("%s.t%s.rep%s.frequency",nsCond,tp,1:params$numReplicates[[nsCond]])
      plotReplicateCorrelation(params, nsReps, info, marginalCounts, marginalTiles)
      
      # CALC MEANS AND NORMALIZE ------------------------------------------------
      
      #pull out nonselect condition and average over replicates
      dreps <- sub("frequency","effectiveDepth",nsReps)
      if (params$numReplicates[[nsCond]] > 1) {
        nsMarginalMeans <- rowMeans(marginalCounts[,nsReps],na.rm=TRUE)
        nsMarginalDepth <- rowMeans(marginalCounts[,dreps],na.rm=TRUE)
        nsAllMeans <- rowMeans(allCounts[,nsReps],na.rm=TRUE)
        nsAllDepth <- rowMeans(allCounts[,dreps],na.rm=TRUE)
      } else {
        nsMarginalMeans <- marginalCounts[,nsReps]
        nsMarginalDepth <- marginalCounts[,dreps]
        nsAllMeans <- allCounts[,nsReps]
        nsAllDepth <- allCounts[,dreps]
      }
      
      smallestFreq <- unique(sort(fin(nsMarginalMeans)))[[2]]
      pseudoCount <- 10^floor(log10(smallestFreq))
      
      #if WT controls are present, average over them as well and subtract from nonselect
      wtCond <- getWTControlFor(nsCond,params)
      # wtCond <- findWTCtrl(params,nsCond)
      if (length(wtCond) > 0) {
        #extract WT reps
        wtReps <- getWTReps(info, params, wtCond)
        
        if (params$numReplicates[[wtCond]] > 1) {
          wtMarginalMeans <- rowMeans(marginalCounts[,wtReps],na.rm=TRUE)
        } else {
          wtMarginalMeans <- marginalCounts[,wtReps]
        }
        
        #adjust pseudocount if necessary
        smallestFreq <- unique(sort(fin(wtMarginalMeans)))[[2]]
        wtPC <- 10^floor(log10(smallestFreq))
        if (wtPC < pseudoCount) pseudoCount <- wtPC
        
        #new WT ctrl plot
        plotWTCtrl(info, params, marginalTiles, wtMarginalMeans, nsMarginalMeans, pseudoCount)

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
      simplifiedMarginal <- buildSimplifiedMarginal(marginalSplitChanges, nsMarginalMeans, 
                                                    nsMarginalDepth, marginalTiles)
      
      # JACKPOT DIAGNOSIS ---------------------------------------------
      
      logInfo("Drawing Jackpot plot")
      plotJackpotDiagnosis(info, nsMarginalMeans, wmThreshold, marginalCounts)
      
      # FRAMESHIFT HOTSPOT MAP ---------------------------------------------
      
      fsIdx <- grep("fs$|del$|ins\\w+$",marginalCounts$hgvsp)
      if (length(fsIdx) != 0) {
        plotFrameshiftHotspot(info, params, fsIdx, marginalSplitChanges, 
                              nsMarginalMeans, marginalCounts)
      }
      
      # NUCL BIAS ANALYSIS -------------------------------------------
      
      logInfo("Checking nucleotide distribution")
      drawNucleotideBias(info, params, simplifiedMarginal, tile)
      
      #calculate global filters that can be combined with tile-specific filters below
      isFrameshift <- grepl("fs$",allCounts$aaChanges)
      isInFrameIndel <- sapply(allSplitChanges, function(changes) {
        any(changes[,3] != "indel" & nchar(changes[,3]) != 3)
      })
      numChanges <- sapply(allSplitChanges, nrow)
      maxChanges <- max(numChanges)
      
      # CENSUS -------------------------------------------------
      
      logInfo("Calculating variant census")
      #Census per tile
      tileCensi <- do.call(rbind,lapply(params$tiles[,1], tileCensus, 
                                        allTiles=allTiles, 
                                        maxChanges=maxChanges, 
                                        nsAllMeans=nsAllMeans, 
                                        isFrameshift=isFrameshift, 
                                        isInFrameIndel=isInFrameIndel, 
                                        numChanges=numChanges))
      rownames(tileCensi) <- params$tiles[,1]
      
      #calculate lambda and Poisson fits for each census
      tileLambdas <- do.call(rbind,lapply(rownames(tileCensi), 
                                          calculateFits, 
                                          tileCensi=tileCensi))
      rownames(tileLambdas) <- params$tiles[,1]
      
      # OVERALL CENSUS ----------------------------------------------------------
      
      logInfo("Extrapolating variant censi for each region")
      #determine which tiles are in which regions
      tilesPerRegion <- tilesInRegions(params)
      pdf(paste0(outDir,nsCond,"_t",tp,"_census.pdf"),8.5,11)
      opar <- par(mfrow=c(4,3),oma=c(2,2,2,2))
      tagger <- pdftagger(paste(pdftag,"; condition:",nsCond,"; time point:",tp),cpp=12)
      
      regionCensi <- lapply(names(tilesPerRegion), regionCensus, 
                            tilesPerRegion=tilesPerRegion, 
                            tileCensi=tileCensi, 
                            tileLambdas=tileLambdas, 
                            maxChanges=maxChanges, 
                            tagger=tagger)
      par(opar)
      invisible(dev.off())
      
      # COVERAGE MAP -------------------------------------------------------------
      
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
          freq=sum(freq[is]),
          depth=mean(depth[is],na.rm=TRUE)
        )
      })))
      
      drawCoverageMap(params, info, aaMarginal, wmThreshold, tilesPerRegion, 
                      tileDepths, depthTable, tileCensi, tileLambdas)
      
      # COMPLEXITY ANALYSIS ------------------------------------------------------
      
      logInfo("Running complexity analysis.")
      
      cplxPlot <- analyzeComplexity(allCounts, simplifiedMarginal)
      drawComplexity(info, cplxPlot)
      
      # WELL-MEASUREDNESS ANALYSIS -----------------------------------------------
      
      logInfo("Running Well-measuredness analysis")
      reachable <- reachableChanges(params)
      
      # data(trtable)
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
      #TODO: Make this dynamic based on pseudocount value (which indicates lowest magnitude)
      thresholds <- 10^seq(-7,-2,0.1)
      
      #iterate over regions and analyze separately
      pdf(paste0(outDir,nsCond,"_t",tp,"_wellmeasured.pdf"),8.5,11)
      opar <- par(mfrow=c(3,2),oma=c(2,2,2,2))
      tagger <- pdftagger(paste(pdftag,"; condition:",nsCond,"; time point:",tp),cpp=6)
      coverageCurves <- lapply(names(tilesPerRegion), function(ri) {
        
        tiles <- tilesPerRegion[[ri]]
        # rStart <- min(params$tiles[tiles,"Start AA"])
        rStart <- min(params$tili(tiles)[,"Start AA"])
        # rEnd <- max(params$tiles[tiles,"End AA"])
        rEnd <- max(params$tili(tiles)[,"End AA"])
        rLength <- rEnd-rStart+1
        
        #if there's no data, return an empty plot
        if (!any(aaMarginal$tile %in% tiles)) {
          plot.new()
          rect(0,0,1,1,col="gray80",border="gray30",lty="dotted")
          text(0.5,0.5,"no data")
          mtext(paste0("Region #",ri))
          tagger$cycle()
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
        #export coverage table
        outTblFile <- paste0(outDir,nsCond,"_t",tp,"_wm_region",ri,".csv")
        write.csv(coverageCurves,outTblFile,row.names=FALSE)
        
        #draw the plot
        plotCoverageCurve(coverageCurves, thresholds, ri, wmThreshold, tagger)
        return(coverageCurves)
        
      })
      par(opar)
      invisible(dev.off())
      
      # MUT-TYPE ANALYSIS --------------------------------------------------------
      
      logInfo("Plotting mutation type breakdown per tile")
      simplifiedMarginal$fromaa <- sapply(simplifiedMarginal$from, 
                                          function(codon) if (is.na(codon)) "" else trtable[[codon]]
      )
      mutbreakdown <- do.call(cbind,lapply(params$tiles[,1], mut_breakdown, 
                                           simplifiedMarginal=simplifiedMarginal,
                                           freq=freq, 
                                           fromaa=fromaa, 
                                           toaa=toaa, 
                                           to=to))
      colnames(mutbreakdown) <- params$tiles[,1]
      
      plotMutbreakdown(info, mutbreakdown)
      
    }
  }
  
  options(op)
  logInfo("QC analysis complete.")
  return(NULL)
}


# HELPER FUNCTIONS -------------------------------------------------------------

#' Get all nonselect conditions from parameter sheet object
#'
#' @param params parameter sheet object
#' @param allCondOverride override flag to run on ALL conditions instead of just
#' nonselect
#' 
#' @return list of nonselect conditions
#' @export
#' 
getNsConditions <- function(params, allCondOverride) {
  # identify nonselect conditions
  nsConditions <- getNonselects(params)
  
  # apply override if requested
  if (allCondOverride) {
    nsConditions <- params$conditions$names
  }
  
  # if this is a pure QC run, there are no conditions declared as nonselect, so
  # we treat *all* conditions as nonselect.
  if (length(nsConditions) == 0) {
    if (nrow(params$conditions$definitions) > 0) {
      nsConditions <- params$conditions$definitions[,3]
    } else {
      nsConditions <-  params$conditions$names
    }
  }
  return(nsConditions)
}


#' Output alignment quality plot with percentages of unmapped or mismapped reads
#' per sample    
#'
#' @param params parameter sheet object
#' @param depthTable raw depth table
#' @param inDir input directory
#' @param outDir output directory
#' @param pdftag PDF tag string
#'
#' @return NULL
#' @export
#' 
countAttributes <- function(depthTable,inDir,outDir,pdftag) {
  countHeaders <- lapply(depthTable$countfile, function(cfile) {
    #rebase
    cfile <- paste0(inDir,basename(cfile))
    lines <- readLines(cfile,20)
    values <- do.call(c,lapply(
      strsplit(lines[grep("^#",lines)],":\\s*"),
      function(xs) setNames(
        paste(xs[-1],collapse=":"),
        sub("^#","",trimws(xs[[1]]))
      )
    ))
    values <- values[names(values) != "Comment"]
    values <- as.list(values)
    isNum <- suppressWarnings(which(!is.na(as.numeric(values))))
    values[isNum] <- as.numeric(values[isNum])
    values
  })
  commonFields <- Reduce(intersect,lapply(countHeaders,names))
  countHeaders <- as.df(lapply(countHeaders,`[`,commonFields))
  rawdepth <- countHeaders$`Raw read depth`
  unmapped <- countHeaders$`Number of read pairs did not map to gene`
  offtile <- countHeaders$`Number of reads outside of the tile`
  
  toPlot <- rbind(
    unmapped=100*unmapped/rawdepth,
    offtile=100*offtile/rawdepth
  )
  
  pdffile <- paste0(outDir,"seqReads.pdf")
  pdf(pdffile,8.5,11)
  opar <- par(las=1,mar=c(5,10,1,1),oma=c(2,2,2,2),cex=.6)
  tagger <- pdftagger(paste(pdftag),cpp=1)
  barplot(
    toPlot,names.arg=depthTable$Sample.ID,
    beside=TRUE,horiz=TRUE,border=NA,
    xlab="% reads",xlim=c(0,max(toPlot)*1.2),
    col=c(1,2)
  )
  grid(NULL,NA)
  tagger$cycle()
  legend("right",c("no match","wrong tile"),fill=c(1,2))
  par(opar)
  dev.off()
}


#' Output basecall quality plot with observed error rate at given basecall PHRED 
#' quality score
#'
#' @param inDir input directory
#' @param outDir output directory
#' @param pdftag  PDF tag string
#'
#' @return NULL
#' @export
#' 
plotSeqErrorRates <- function(inDir,outDir,pdftag) {
  
  phredFiles <- list.files(inDir,pattern="calibrate_phred",full.names = TRUE)
  
  if (length(phredFiles) > 0) {
    sampleIDs <- sub("_calibrate_phred.csv$","",basename(phredFiles))
    phredTables <- setNames(lapply(phredFiles,read.csv),sampleIDs)  
    
    allPhreds <- cbind(phredTables[[1]][,1:2],do.call(cbind,lapply(phredTables,function(tbl)tbl$observed)))
    nafilter <- apply(allPhreds[,-(1:2)],1,function(x)!all(is.na(x)))
    
    # labels <- na.omit(allPhreds)[,1]
    # probs <- na.omit(allPhreds)[,-1]
    probs <- as.matrix(allPhreds[nafilter,-1][-1,])
    
    pchars <- allPhreds[nafilter,1][-1]
    qnums <- sapply(pchars,function(x) as.integer(charToRaw(x))-33)
    labels <- sprintf("Q%d (='%s')",qnums,pchars)
    
    #any zero probabilities are the result of undersampling, so we remove them
    probs[probs==0] <- NA
    
    if (nrow(probs)==0) {
      logWarn("Invalid PHRED error rate data. Skipping analysis...")
      return(NULL)
    }
    
    pdffile <- paste0(outDir,"seqErrorRates.pdf")
    pdf(pdffile,11,8.5)
    tagger <- pdftagger(paste(pdftag),cpp=1)
    op <- par(mar=c(5,4,1,1),oma=c(2,2,2,2))
    plotCols <- gray.colors(length(sampleIDs)+1)
    barplot(t(probs),
            beside=TRUE,col=plotCols,border=NA,ylab="error rate",
            xlab="Quality score",names.arg=labels,log="y",ylim=c(1e-5,1)
    )
    grid(NA,NULL)
    legend("topright",c("PHRED spec",sampleIDs),fill=plotCols,border=NA)
    tagger$cycle()
    par(op)
    dev.off()
    
  } else {
    logWarn("No PHRED calibration files found. Skipping analysis...")
  }
}


#' Output effective depth plot with effective sequencing depth across the tiles 
#' for each condition and replicate
#'
#' @param covTable dataframe with position-specific effective depth at each 
#' nucleotide position in each sample
#' @param outDir output directory
#' @param pdftag PDF tag string
#'
#' @return NULL
#' @export
#' 
drawPositionalDepth <- function(covTable,outDir,pdftag) {
  pos <- as.integer(sub("^X","",colnames(covTable)))
  
  #this may be a problem if there are more than 9 conditions...
  idxCols <- c("steelblue","firebrick","chartreuse","gold","purple","slategray",
               "royalblue","hotpink","springgreen","orange","darkorchid")
  
  conds <- sort(rownames(covTable))
  condGroups <- sub("\\.rep\\d+","",conds)
  cgTable <- table(condGroups)
  lineColors <- setNames(do.call(c,lapply(1:length(cgTable), function(j) {
    colorRampPalette(paste0(idxCols[[j]],c(1,4)))(cgTable[[j]])
  })),conds)
  
  pdf(paste0(outDir,"effectiveSeqDepth.pdf"),11,8.5)
  tagger <- pdftagger(pdftag,cpp=1)
  opar <- par(mar=c(5,4,1,1),oma=c(2,2,2,2))
  plot(NA,type="n",
       xlim=c(min(pos),1.2*max(pos)),ylim=c(0,max(covTable,na.rm = TRUE)),
       xlab="CDS nucl. position" ,ylab="Effective seq. depth"   
  )
  invisible(lapply(1:nrow(covTable),function(i) {
    lines(pos,covTable[i,],col=lineColors[rownames(covTable)[[i]]])
  }))
  legend("right",names(lineColors),col=lineColors,lty=1)
  tagger$cycle()
  par(opar)
  invisible(dev.off())
}


#' Adjust dataframe 'covTable' with mean "effective" sequencing depth
#'
#' @param covTable dataframe with position-specific effective depth at each 
#' nucleotide position in each sample
#' @param params parameter sheet object
#' 
#' @return table with adjusted mean effective depths
#' @export
#' 
adjustTableDepth <- function(covTable, params) {
  do.call(rbind,lapply(1:nrow(covTable),function(row) {
    apply(params$tiles,1,function(tile) {
      poss <- as.character(tile[["Start NC in CDS"]]:tile[["End NC in CDS"]])
      poss <- intersect(colnames(covTable),poss)
      mean(as.matrix(covTable[row,poss]),na.rm=TRUE)
    })
  }))
}


#' Extract affected positions and codon details from vector 'css'
#'
#' @param ccs vector containing codon changes
#' 
#' @return  table containing start amino acid, position, and end amino acid
#' @export
#' 
splitCodonChanges <- function(ccs) {
  #This is a faster implementation that should require less memory and runtime
  rawTable <- do.call(rbind,lapply(ccs,function(cc) {
    #convert to ASCII codes
    codes <- as.integer(charToRaw(cc))
    #codes <= 57 are numbers, above are letters.
    numIdx <- which(codes <= 57)
    numStart <- min(numIdx)
    numEnd <- max(numIdx)
    c(substr(cc,1,numStart-1),
      substr(cc,numStart,numEnd),
      substr(cc,numEnd+1,nchar(cc))
    )
  }))
  rawTable
}


#' Determine which tiles the positions from split codon changes table 'cct' 
#' belong to 
#'
#' @param params parameter sheet object
#' @param cct table with split codon changes
#' 
#' @return vector of corresponding tile numbers for each codon change
#' @export
#'
inferTiles <- function(params, cct) {
  # sapply(cct$pos,function(pos) max(which(tileStarts <= pos)))
  # params$pos2tile(cct$pos)
  params$pos2tile(as.integer(cct[,2]))
}


#' Calculate relative drops in effective sequencing depth for each variant
#'
#' @param xCounts marginal or all counts table
#' @param xTiles tile assignments for each variant in the table
#' @param depthTable raw depth table
#' @param mc.cores number of processes to run in parallel for multi-core processing
#'
#' @return a table listing the relative drop for each variant in each
#' condition-time-replicate combination
#' @export
#'
calcDepthDrops <- function(xCounts, xTiles, depthTable, mc.cores=6) {
  edcols <- grep("effectiveDepth$",colnames(xCounts))
  ednames <- colnames(xCounts)[edcols]
  edconds <- extract.groups(
    colnames(xCounts)[edcols],
    "^([A-Za-z0-9]+)\\.t([A-Za-z0-9]+)\\.rep(\\d+)\\.effectiveDepth$"
  )
  uniqueTiles <- sort(unique(depthTable$Tile.ID))
  rawDepths <- do.call(rbind,lapply(1:length(edcols),function(j) {
    cond <- edconds[j,1]
    tp <- edconds[j,2]
    repl <- edconds[j,3]
    ds <- sapply(uniqueTiles, function(tile) {
      k <- with(depthTable,which(
        Tile.ID==tile & Condition==cond & 
          Time.point==tp & Replicate==repl
      ))
      if (length(k) != 0) {
        depthTable[k,"alignedreads"]
      } else 0
    })
    setNames(ds,uniqueTiles)
  }))
  do.call(rbind,pbmclapply(1:nrow(xCounts), function(i) {
    # tile <- xTiles[[i]][[1]]
    tile <- names(which.max(table(xTiles[[i]])))
    if (is.null(tile) || !(as.numeric(tile) %in% uniqueTiles)) {
      return(rep(NA,length(edcols)))
    }
    setNames(sapply(1:length(edcols), function(j){
      ed <- xCounts[i,edcols[[j]]]
      rd <- rawDepths[j,as.character(tile)]
      sapply(1 - ed/rd,max,0)
    }),ednames)
  },mc.cores=mc.cores))
}


#' Output depth drop distribution plot showing distribution of relative 
#' drops in effective sequencing depth
#'
#' @param params parameter sheet object
#' @param info list containing details on condition (nscond), timepoint (tp), 
#' pdftag and outdir
#' @param allDepthDrops table listing the relative drop for all variants in 
#' each condition-time-replicate combination 
#' 
#' @return NULL
#' @export
#' 
plotDepthDrops <- function(params, outDir, pdftag, allDepthDrops) {
  pdffile <- paste0(outDir,"depthDrops.pdf")
  pdf(pdffile,11,8.5)
  tagger <- pdftagger(paste(pdftag),cpp=1)
  op <- par(mar=c(5,4,1,1),oma=c(2,2,2,2))
  hist(100*allDepthDrops,breaks=seq(0,100),
       col="gray",border=NA,main="",
       xlab="% drop in effective depth\nper variant"
  )
  abline(v=params$varcaller$maxDrop*100,col="red",lty="dashed")
  tagger$cycle()
  par(op)
  invisible(dev.off())
}


#' Filter variants that have greater depth drop than max allowed depth drop or 
#' do not have minimum required depth 
#'
#' @param nsCond the nonselect condition 
#' @param params parameter sheet object
#' @param xDepthDrops table listing the relative drop for variants in xCounts table
#' @param xCounts marginal or all counts table
#' @param xSplitChanges table with split codon changes for variants in xCounts
#' @param xTiles tile assignments for each variant in xCounts table
#' @param mc.cores number of processes to run in parallel for multi-core processing
#' @param allVariants flag, whether the xCounts table contains all variants
#'
#' @return NULL
#' @export
#'
filterVariants <- function(nsCond, params, xDepthDrops, xCounts, xSplitChanges, xTiles, mc.cores, allVariants=TRUE) {
  
  regex <- paste0(c(nsCond,getWTControlFor(nsCond,params)),".*effectiveDepth",collapse="|")
  aasToCover <- median(params$tiles[,"End AA"]-params$tiles[,"Start AA"])*20
  minDepth <- params$scoring$countThreshold*aasToCover
  
  relCols <- grep(regex,colnames(xDepthDrops))
  filter1 <- sapply(1:nrow(xCounts), function(i) {
    any(is.na(xDepthDrops[i,relCols])) || any(xDepthDrops[i,relCols] > params$varcaller$maxDrop)
  })
  
  relCols <- grep(regex,colnames(xCounts))
  filter2 <- do.call(c,pbmclapply(1:nrow(xCounts), function(i) {
    any(is.na(xCounts[i,relCols])) || any(xCounts[i,relCols] < minDepth)
  },mc.cores=mc.cores))
  
  filter <- filter1 | filter2
  
  if (any(filter)) {
    
    if (allVariants && sum(filter)/length(filter) > 0.1) {
      logWarn("
  #########################################
  MASSIVE SEQUENCING DATA ANOMALY DETECTED!
      ! These results may be unusable !
  #########################################"
      )
    }
    
    logWarn(sprintf(
      "Removing %d variants (%.02f%%) with insufficient effective depth from analysis",
      sum(filter),100*sum(filter)/length(filter)
    ))
    logWarn("Examples of dropped variants :\n   ",paste(head(xCounts$hgvsc[filter]),collapse="\n   "))
    
    xCounts <- xCounts[!filter,]
    xSplitChanges <- xSplitChanges[!filter]
    xTiles <- xTiles[!filter]
  }
}


#' Output replicate correlation plot between technical replicates from nonselect condition
#'
#' @param params parameter sheet object
#' @param nsReps number of nonselect replicates
#' @param info list containing details on condition (nscond), timepoint (tp), 
#' pdftag and outdir
#' @param marginalCounts marginal counts table
#' @param marginalTiles table with tile assignments for each entry in marginal 
#' counts table
#'
#' @return NULL
#' @export
#' 
plotReplicateCorrelation <- function(params, nsReps, info, marginalCounts, marginalTiles) {
  nreps <- params$numReplicates[[info$nsCond]]
  
  if (nreps == 2) {
    
    logInfo("Drawing per-tile replicate correlation plot")
    pdf(paste0(info$outDir,info$nsCond,"_t",info$tp,"_tileRepCorr.pdf"),8.5,11)
    tagger <- pdftagger(paste(info$pdftag,"; condition:",info$nsCond,"timepoint:",info$tp),cpp=12)
    opar <- par(mfrow=c(4,3),oma=c(2,2,2,2))
    for (tilei in 1:nrow(params$tiles)) {
      tile <- params$tiles[tilei,"Tile Number"]
      subset <- marginalCounts[which(marginalTiles == tile),]
      if (nrow(subset) > 0) {
        r <- cor(fin(log10(subset[,nsReps[1:2]]+1e-7)))[1,2]
        rtxt <- sprintf("R = %.02f",r)
        plot(
          subset[,nsReps[[1]]]+1e-7,subset[,nsReps[[2]]]+1e-7,
          log="xy",pch=20,col=colAlpha(1,0.2),
          main=sprintf("Tile %d\nR = %.02f",tile,r),
          xlab=paste(info$nsCond,"rep.1"),ylab=paste(info$nsCond,"rep.2")
        )
      } else {
        plot.new()
        rect(0,0,1,1,col="gray80",border="gray30",lty="dotted")
        text(0.5,0.5,"no data")
        mtext(paste("Tile",tile),side=3)
      }
      tagger$cycle()
    }
    par(opar)
    invisible(dev.off())
    
  } else if (nreps > 2) {
    
    panel.cor <- function(x, y,...){
      usr <- par("usr")
      par(usr=c(0,1,0,1))
      r <- cor(fin(cbind(x,y)))[1,2]
      txt <- sprintf("R = %.02f",r)
      cex.cor <- r*0.8/strwidth(txt)
      text(.5,.5,txt,cex=cex.cor)
      points(mean(usr[1:2]), mean(usr[3:4]))
      par(usr)
    }
    
    pdf(paste0(info$outDir,info$nsCond,"_t",info$tp,"_tileRepCorr.pdf"),8.5,11)
    tagger <- pdftagger(paste(info$pdftag,"; condition:",info$nsCond,"timepoint:",info$tp),cpp=1)
    for (tilei in 1:nrow(params$tiles)) {
      tile <- params$tiles[tilei,"Tile Number"]
      subset <- marginalCounts[which(marginalTiles == tile),]
      if (nrow(subset) > 0) {
        pairs(
          log10(subset[,nsReps]+1e-7), labels=sprintf("repl. %d",1:nreps),
          pch=20,col=colAlpha(1,0.2),main=paste("Tile",tile),
          lower.panel=panel.cor,oma=c(12,4,10,4)
          # xlab=paste(info$nsCond,"rep.1"),ylab=paste(info$nsCond,"rep.2")
        )
      } else {
        opar <- par(oma=c(2,2,2,2))
        plot.new()
        rect(0,0,1,1,col="gray80",border="gray30",lty="dotted")
        text(0.5,0.5,"no data")
        mtext(paste("Tile",tile),side=3)
        par(opar)
      }
      opar <- par(oma=c(2,2,2,2))
      tagger$cycle()
      par(opar)
    }
    par(opar)
    invisible(dev.off())
    
  }
}


#' Extract appropriate WT replicates for timepoint given in 'info'
#'  
#' @param info list containing details on condition (nscond), timepoint (tp), 
#' pdftag and outdir
#' @param params parameter sheet object
#' @param wtCond 
#'
#' @return vector of strings with corresponding WT replicates 
#' @export
#'

getWTReps <- function(info, params, wtCond) {
  #Find the appropriate time point for the matched WT control
  wtTps <- getTimepointsFor(wtCond,params)
  wtTp <- if (info$tp %in% wtTps) info$tp else wtTps[[1]]
  #extract WT reps
  wtReps <- sprintf("%s.t%s.rep%s.frequency",wtCond,wtTp,1:params$numReplicates[[wtCond]])
  
  return(wtReps)
}


#' Draw WT frequency plot comparing distribution of marginal frequencies in 
#' nonselect condition against WT condition 
#'
#' @param info list containing details on condition (nscond), timepoint (tp), 
#' pdftag and outdir
#' @param params parameter sheet object
#' @param marginalTiles table with tile assignments for each entry in marginal 
#' counts table
#' @param wtMarginalMeans table with mean marginal frequency across technical 
#' replicates in WT condition
#' @param nsMarginalMeans table with mean marginal frequency across technical 
#' replicates in nonselect condition
#' @param pseudoCount a small number added to each frequency to allow for log-scale visualization
#'
#' @return NULL
#' @export
#' 
plotWTCtrl <- function (info, params, marginalTiles, wtMarginalMeans, nsMarginalMeans, pseudoCount) {
  logInfo("Drawing WT control level plot")
  pdf(paste0(info$outDir,info$nsCond,"_t",info$tp,"_WTlevels.pdf"),8.5,11)
  tagger <- pdftagger(paste(info$pdftag,"; condition:",info$nsCond,"timepoint:",info$tp),cpp=6)
  opar <- par(mfrow=c(3,2),oma=c(2,2,2,2))
  
  for (tile in params$tiles[,1]) {
    rows <- which(marginalTiles == tile)
    if (length(rows) == 0) {
      plot.new()
      rect(0,0,1,1,col="gray80",border="gray30",lty="dotted")
      text(0.5,0.5,"no data")
      mtext(paste0("Tile #",tile),side=4)
    } else {
      breaks <- seq(log10(pseudoCount),0,.1)
      wtHist <- hist(log10(wtMarginalMeans[rows]+pseudoCount),breaks=breaks,plot=FALSE)
      nsHist <- hist(log10(nsMarginalMeans[rows]+pseudoCount),breaks=breaks,plot=FALSE)
      maxDens <- max(c(wtHist$density,nsHist$density),na.rm=TRUE)
      plot(NA,type="n",xlim=c(log10(pseudoCount),0),ylim=c(-maxDens,maxDens),axes=FALSE,
           ylab="wildtype density : nonselect density",xlab="mean marginal frequency",
           main=sprintf("Tile #%i",tile)
      )
      axis(1,at=seq(log10(pseudoCount),0),c(0,10^((log10(pseudoCount)+1):-1),1))
      axis(2,at=-round(maxDens):round(maxDens),c(round(maxDens):0,1:round(maxDens)))
      with(nsHist,rect(breaks[-length(breaks)],0,breaks[-1],density,col="steelblue3",border=NA))
      with(wtHist,rect(breaks[-length(breaks)],0,breaks[-1],-density,col="firebrick3",border=NA))
      grid(NULL,NULL)
      abline(h=0)
    }
    tagger$cycle()
  }
  
  par(opar)
  invisible(dev.off())
}

#' Build a simplified version of the marginal table with columns for variant, mean marginal  
#' frequency, marginal depth and tile, with 0-depth variants filtered out
#'
#' @param marginalSplitChanges table with split codon changes from marginal counts table
#' @param nsMarginalMeans table with mean marginal frequency across nonselect technical replicates 
#' @param nsMarginalDepth table with mean marginal depth across nonselect technical replicates
#' @param marginalTiles table with tile assignments for each entry in marginal 
#' counts table
#'
#' @return table with variant, mean marginal frequency, marginal depth and tile 
#' @export
#' 
buildSimplifiedMarginal <- function(marginalSplitChanges, nsMarginalMeans, nsMarginalDepth, marginalTiles) {
  colnames(marginalSplitChanges) <- c("from","pos","to")
  simplifiedMarginal <- cbind(
    as.data.frame(marginalSplitChanges),
    freq=nsMarginalMeans,depth=nsMarginalDepth,tile=marginalTiles
  )
  simplifiedMarginal$pos <- as.integer(simplifiedMarginal$pos)
  
  #filter out 0-depth variants here?
  if (any(is.na(simplifiedMarginal$freq)) || any(simplifiedMarginal$depth < 1)) {
    culprits <- which(simplifiedMarginal$depth < 1 | is.na(simplifiedMarginal$freq))
    logWarn(sprintf("Discarding %d variants with 0 depth!",length(culprits)))
    simplifiedMarginal <- simplifiedMarginal[-culprits,]
  }
  
  return(simplifiedMarginal)
}


#' Output jackpot plot with marginal frequencies of all variants in decreasing order
#'
#' @param nsMarginalMeans table with mean marginal frequency across nonselect technical replicates
#' @param wmThreshold well-measuredness threshold
#' @param marginalCounts table of marginal counts
#'
#' @return NULL
#' @export
#' 
plotJackpotDiagnosis <- function(info, nsMarginalMeans, wmThreshold, marginalCounts) {
  
  pdf(paste0(info$outDir,info$nsCond,"_t",info$tp,"_jackpot.pdf"),11,8.5)
  tagger <- pdftagger(paste(info$pdftag,"; condition:",info$nsCond,"timepoint:",info$tp),cpp=1)
  opar <- par(oma=c(2,2,2,2),mar=c(7,6,3,3)+.1)
  
  decOrder <- order(nsMarginalMeans,decreasing=TRUE)
  cm <- yogitools::colmap(
    valStops = c(0,wmThreshold/10,wmThreshold,wmThreshold*100,1e-2),
    colStops = c("firebrick3","gold","darkolivegreen3","gold","firebrick3")
  )
  plotcol <- cm(nsMarginalMeans[decOrder])
  plot(
    seq(0,1,length.out=length(decOrder)),
    nsMarginalMeans[decOrder],pch=20,col=plotcol,
    ylab="marginal frequency",xlab="fraction of variants"
  )
  top10 <- decOrder[1:10]
  steps <- seq(0.1,0.9,length.out=5)
  arrows(steps,nsMarginalMeans[top10],(1:10)/length(decOrder),nsMarginalMeans[top10],length=0.05,col="gray")
  text(steps,nsMarginalMeans[top10],marginalCounts$hgvsp[top10],pos=4,cex=.8)
  
  tagger$cycle()
  par(opar)
  invisible(dev.off())
}


#' Output frameshift hotspot plot with distribution of frameshift variants across 
#' template sequence 
#'
#' @param info list containing details on condition (nscond), timepoint (tp), 
#' pdftag and outdir
#' @param params parameter sheet object
#' @param fsIdx frameshift indices from marginal counts table
#' @param marginalSplitChanges table with split codon changes from marginal counts table
#' @param nsMarginalMeans table with mean marginal frequency across nonselect technical replicates
#' @param marginalCounts table of marginal counts
#' 
#' @return NULL
#' @export
#' 
plotFrameshiftHotspot <- function(info, params, fsIdx , marginalSplitChanges, nsMarginalMeans, marginalCounts) {
  
  idxPlusCoord <- do.call(rbind,tapply(fsIdx, marginalSplitChanges[fsIdx,2], function(idx) {
    x <- as.numeric(marginalSplitChanges[idx,2]) + seq_along(idx)/(length(idx)+1)
    cbind(idx,x)
  }))
  
  #re-order according to x-coordinate
  fsIdx <- idxPlusCoord[,1]
  xCoords <- idxPlusCoord[,2]
  yCoords <- nsMarginalMeans[fsIdx]
  
  top10 <- head(fsIdx[order(yCoords,decreasing=TRUE)],10)
  top10X <- xCoords[sapply(top10,function(i)which(fsIdx==i))]
  top10Y <- nsMarginalMeans[top10]
  top10Labels <- marginalCounts$hgvsc[top10]
  
  cm <- colmap(valStops = c(0,5e-4,5e-3),colStops = c("black","gold","firebrick3"))
  plotcols <- sapply(cm(yCoords),colAlpha,0.5)
  
  # logInfo("Drawing frameshift map")
  pdf(paste0(info$outDir,info$nsCond,"_t",info$tp,"_fsMap.pdf"),11,8.5)
  tagger <- pdftagger(paste(info$pdftag,"; condition:",info$nsCond,"timepoint:",info$tp),cpp=1)
  opar <- par(oma=c(2,2,2,2),mar=c(5,5,0,1)+.1)
  layout(rbind(1,2),heights=c(0.1,0.9))
  # plotcols <- sapply(1:2,colAlpha,0.5)
  # funs <- list(mean=mean,median=median)
  
  xRange <- c(min(params$regions[,"Start AA"]),max(params$regions[,"End AA"]))
  
  par(mar=c(0,5.1,0,1.1))
  plot(NA,type="n",xlim=xRange, ylim=c(0,1),xlab="",ylab="",axes=FALSE)
  rect(params$regions[,"Start AA"],0,params$regions[,"End AA"],0.49,col="gray80",border=NA)
  text(rowMeans(params$regions[,c("Start AA","End AA"),drop=FALSE]),0.25,params$regions[,"Region Number"])
  rect(params$tiles[,"Start AA"],0.51,params$tiles[,"End AA"],1,col="gray90",border=NA)
  text(rowMeans(params$tiles[,c("Start AA","End AA"),drop=FALSE]),0.75,params$tiles[,"Tile Number"])
  
  par(mar=c(5,5,0,1)+.1)
  plot(NA,type="n",log="y",
       xlim=xRange,
       ylim=c(1e-7,1e-1),
       xlab="AA position",ylab="frameshift marginal freq."
  )
  segments(xCoords,1e-7,xCoords,yCoords+1e-7,col=plotcols)
  grid(NA,NULL)
  text(top10X,top10Y,top10Labels,cex=0.5,pos=3)
  
  tagger$cycle()
  par(opar)
  invisible(dev.off())
}


#' Helper function to perform nucleotide bias analysis and draw the corresponding plot
#'
#' @param simplifiedMarginal simplified marginal frequency table
#' @param tile the tile number being analyzed
#' @param draw flag, whether to draw nucleotide bias analysis plot
#'
#' @return rates at which each nucelotide is substituted for another across the tile
#' 
nucleotideBiasAnalysis <- function(simplifiedMarginal,tile,draw=TRUE) {
  
  counters <- replicate(2,array(0, dim=c(4,4,3),
                                dimnames=list(toChars("ACGT"), toChars("ACGT"), 1:3)
  ), simplify=FALSE)
  names(counters) <- c("single","multi")
  
  bigCounters <- c(single=0,multi=0)
  
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
    freq <- simplifiedMarginal[i,"freq"]
    bigCounters[type] <- bigCounters[type]+freq
    for (j in 1:3) {
      from <- substr(simplifiedMarginal[i,"from"],j,j)
      to <- substr(simplifiedMarginal[i,"to"],j,j)
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
  
  bigRates <- bigCounters/sum(bigCounters)
  
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
    bigLabels <- sprintf("%s (%.01f%%)",c("SNV","MNV"),100*bigRates)
    text(c(7,23),4.66,bigLabels,cex=0.9)
    axis(2,at=(4:1)-.5,toChars("ACGT"))
    axis(1,at=c(1:4,6:9,11:14,17:20,22:25,27:30)-.5,rep(toChars("ACGT"),6))
    mtext(paste0("Tile #",tile),side=4)
    par(op)
  }
  
  return(rates)
}


#' Draw Nucleotide Bias Analysis plots
#'
#' @param info list containing details on condition (nscond), timepoint (tp), 
#' pdftag and outdir
#' @param params parameter sheet object 
#' @param simplifiedMarginal simplified marginal frequency table
#' @param tile the tile number being analyzed
#' 
#' @return NULL
#' @export
#' 
drawNucleotideBias <- function(info, params, simplifiedMarginal,tile){
  pdf(paste0(info$outDir,info$nsCond,"_t",info$tp,"_nucleotide_bias.pdf"),8.5,11)
  opar <- par(mfrow=c(6,1),oma=c(2,2,2,2))
  #set up pdf tagger
  tagger <- pdftagger(paste(info$pdftag,"; condition:",info$nsCond,"timepoint:",info$tp),cpp=6)
  nuclRates <- lapply(params$tiles[,"Tile Number"], function(tile) {
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
    #add one pdf tag per page
    tagger$cycle()
  })
  par(opar)
  invisible(dev.off())
}


#' Perform census on tile to determine what mutations are present
#'
#' @param tile the tile number being analyzed
#' @param allTiles table containing all tiles
#' @param maxChanges maximum number of mutations allowed on x-axis of plot
#' @param nsAllMeans table with mean frequency across nonselect technical replicates
#' @param isFrameshift vector containing flags of whether a variant is a frameshift
#' @param isInFrameIndel vector containing flags of whether a variant is an indel
#' @param numChanges vector containing number of nucleotide changes in a given mutational context
#' 
#' @return table containing population breakdown of variant types within tile
#' @export
#' 
tileCensus <- function(tile, allTiles, maxChanges, nsAllMeans, isFrameshift, isInFrameIndel, numChanges) {
  #filter for entries in given tile
  isInTile <- sapply(allTiles,function(tiles) !any(is.na(tiles)) && tile %in% tiles)
  if (!any(isInTile)) {
    return(c(fs=NA,indel=NA,WT=NA,setNames(rep(NA,maxChanges),1:maxChanges)))
  }
  # WT frequency
  wtFreq <- 1 - sum(nsAllMeans[isInTile],na.rm=TRUE)
  #frameshift freq
  fsFreq <- sum(nsAllMeans[isInTile & isFrameshift],na.rm=TRUE)
  #indel frequency
  indelFreq <- sum(nsAllMeans[isInTile & !isFrameshift & isInFrameIndel],na.rm=TRUE)
  #substitution mutation census
  isSubst <- which(isInTile & !isFrameshift & !isInFrameIndel)
  census <- tapply(nsAllMeans[isSubst], numChanges[isSubst],sum,na.rm=TRUE)
  census <- c(
    fs=fsFreq,
    indel=indelFreq,
    WT=wtFreq,
    sapply(as.character(1:maxChanges),function(n) if (n %in% names(census)) census[[n]] else 0)
  )
  return(census)
}


#' Calculate lambda and Poisson fits for each census
#'
#' @param tile tile number being analyzed
#' @param tileCensi table containing population breakdown of variant types within each tile
#' 
#' @return vector of lambda and Poisson fits for census
#' @export
#' 
calculateFits <- function(tile, tileCensi) {
  freqs <- tileCensi[tile,-c(1,2)]
  if (any(is.na(freqs))) {
    return(c(lambda=NA,rmsd=NA))
  }
  ns <- 1:length(freqs)-1
  lambda <- sum(ns*(freqs/sum(freqs)))
  rmsd <- sqrt(mean((freqs-dpois(ns,lambda))^2))
  return(c(lambda=lambda,rmsd=rmsd))
}


#' Perform census on overall region to determine what mutations are present
#'
#' @param ri region index; region being analyzed 
#' @param tilesPerRegion a list of vectors containing the tile IDs for each region 'ri'
#' @param tileCensi table containing population breakdown of variant types within each tile
#' @param tileLambdas average number of mutations in each tile
#' @param maxChanges maximum number of mutations allowed on x-axis of plot
#' @param tagger list that contains functions for adding tags to page margins
#' 
#' @return table containing population breakdown of variant types within region
#' @export
#' 
regionCensus <- function(ri, tilesPerRegion, tileCensi, tileLambdas, maxChanges, tagger){
  
  relevantTiles <- as.character(tilesPerRegion[[ri]])
  #if there is no data for any of those tiles
  if (all(is.na(tileCensi[relevantTiles,"WT"]))) {
    plot.new()
    rect(0,0,1,1,col="gray80",border="gray30",lty="dotted")
    text(0.5,0.5,"no data")
    mtext(paste0("Extrapolation for Region #",ri))
    tagger$cycle()
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
  tagger$cycle()
  return(overallCensus)
}


#' Helper function to plot a census dataset
#'
#' @param census table of population breakdown of variants within tile/region
#' @param lamda average number of mutations in tile/region
#' @param d raw depth
#' @param e effective depth 
#' @param main str, title of plot
#' 
#' @return NULL
#'
plotCensus <- function(census,lambda,d=NULL,e=NULL,main="") {
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
  topmid <- 0.2*u[[3]]+0.8*u[[4]]
  uppermid <- 0.5*u[[3]]+0.5*u[[4]]
  lowermid <- 0.8*u[[3]]+0.2*u[[4]]
  text(midx,lowermid,bquote(lambda == .(round(lambda,digits=3))))
  if (!is.null(d)) {
    # text(midx,topmid,paste("#reads =",d))
    text(midx,topmid,bquote(d["raw"] == .(d)))
  }
  if (!is.null(e)) {
    text(midx,uppermid,bquote(d["eff"] %~~% .(e)))
  }
  par(op)
  invisible(return(NULL))
}


#' Output heatmap of marginal frequencies for all possible amino acid changes 
#' across all tiles
#'
#' @param params parameter sheet object
#' @param info list containing details on condition (nscond), timepoint (tp), 
#' pdftag and outdir
#' @param aaMarginal marginal frequency table at amino acid-level 
#' @param wmThreshold 'well-measuredness' threshold 
#' @param tilesPerRegion a list of vectors containing the tile IDs for each region 'ri'
#' @param tileDepths table with adjusted mean effective sequencing depth for each tile
#' @param depthTable raw depth table
#' @param tileCensi table containing population breakdown of variant types within each tile
#' @param tileLambdas average number of mutations in each tile
#'
#' @return NULL
#' @export
#' 
drawCoverageMap <- function(params, info, aaMarginal, wmThreshold, tilesPerRegion, tileDepths, depthTable, tileCensi, tileLambdas) {
  #plan PDF page layout with 3 rows, covering 4 tiles each.
  #each row being subdivided into a top track and bottom track
  #the top track will contain individual histograms, while the bottom track
  #will contain a joint heatmap for the 4 tiles.
  tpr <- 4 #tiles per row
  rpp <- 3 #rows per page
  ntiles <- nrow(params$tiles)
  lastTile <- max(params$tiles[,"Tile Number"])
  #build layout plan, which sequencing tiles go in which row/column
  tileLayout <- lapply(seq(1,ntiles,tpr),function(i) params$tiles[i:min(ntiles,i+tpr-1),"Tile Number"])
  nrows <- length(tileLayout)
  #break down the layout by PDF page
  pageLayout <- lapply(seq(1,nrows,rpp), function(i) tileLayout[i:min(nrows,i+rpp-1)])
  
  #now do the actual plotting.
  pdf(paste0(info$outDir,info$nsCond,"_t",info$tp,"_coverage.pdf"),8.5,11)
  opar <- par(oma=c(2,2,2,2))
  tagger <- pdftagger(paste(info$pdftag,"; condition:",info$nsCond,"; time point:",info$tp),cpp=1)
  invisible(lapply(pageLayout, function(tileSets) {
    
    #build a matrix that indicates the grid positions of plots in order of drawing
    plotIndex <- do.call(rbind,lapply(1:rpp, buildPlotMatrix, tpr=tpr, tileSets=tileSets, 
                                      lastTile=lastTile))
    layout(plotIndex)
    
    #For each row on the current page...
    lapply(tileSets, function(tiles) {
      #plot coverage map
      # startPos <- min(params$tiles[tiles,"Start AA"])
      startPos <- min(params$tili(tiles)[,"Start AA"])
      # endPos <- max(params$tiles[tiles,"End AA"])
      endPos <- max(params$tili(tiles)[,"End AA"])
      # seps <- params$tiles[tiles,"Start AA"][-1]
      seps <- params$tili(tiles)[,"Start AA"][-1]
      coverageSubmap(startPos,endPos,aaMarginal,seps,thresholds=c(wmThreshold/10,wmThreshold*10))
      #and plot the corresponding censi
      lapply(tiles, function(tile) {
        tileRegion <- names(tilesPerRegion)[sapply(tilesPerRegion,function(ts)tile %in% ts)]
        if (exists("tileDepths")) {
          reps <- params$numReplicates[[info$nsCond]]
          edepths <- as.matrix(tileDepths[
            sprintf("%s.t%s.rep%d",info$nsCond,info$tp,1:reps), as.character(tile)
          ])
        } else {
          edepths <- with(depthTable,depth[
            Condition==info$nsCond & Time.point==info$tp & Tile.ID==tile
          ])
        }
        depths <- with(depthTable,alignedreads[
          Condition==info$nsCond & Time.point==info$tp & Tile.ID==tile
        ])
        plotCensus(
          tileCensi[as.character(tile),],
          lambda=tileLambdas[as.character(tile),"lambda"],
          d=if(length(depths)>0) min(depths,na.rm=TRUE) else NULL,
          e=if(length(edepths)>0) min(edepths,na.rm=TRUE) else NULL,
          main=paste0("Tile #",tile)
        )
        mtext(paste0("(Region #",tileRegion,")"),line=0.5,cex=.7)
      })
    })
    tagger$cycle()
  }))
  drawCoverageLegend(wmThreshold,aaMarginal)
  par(opar)
  invisible(dev.off())
}


#' Helper function to build layout for coverage map
#'
#' @param ri region index
#' @param tpr number of tiles per row
#' @param tileSets 
#' @param lastTile number of the last tile in region
#' 
#' @return table that indicates the grid positions of plots for region in order of drawing
#' 
buildPlotMatrix <- function(ri, tpr, tileSets, lastTile){
  bottom <- (ri-1)*5+1
  topRow <- bottom+1:tpr
  bottomRow <- rep(bottom,tpr)
  currTiles <- min(tileSets[[min(ri,length(tileSets))]])-1+1:tpr
  bottomRow[currTiles > lastTile] <- -1
  #position map on bottom of the row, and census plots for the corresponding tiles above
  rbind(topRow, bottomRow)
}


#' Helper function to draw a subsection of the coverage map
#'
#' @param startPos start position for subsection of coverage map
#' @param endPos end position for subsection of coverage map
#' @param aaMarginal marginal frequency table at amino acid-level
#' @param seps tile separator lines
#' @param thresholds heatmap colour transition thresholds
#' 
#' @return NULL
#' 
coverageSubmap <- function(startPos,endPos,aaMarginal,seps=NULL,thresholds=c(1e-6,5e-5)) {
  mapwidth <- endPos-startPos+2
  aas <- toChars("AVLIMFYWRHKDESTNQCGP*")
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


#' Draw legend for coverage maps
#'
#' @param wmThreshold well-measuredness threshold
#' @param aaMarginal marginal frequency table at amino acid-level
#' 
#' @return NULL
#' 
drawCoverageLegend <- function(wmThreshold,aaMarginal) {
  cmap <- yogitools::colmap(log10(c(wmThreshold/10,wmThreshold*10,1)),c("white","orange","firebrick3"))
  op <- par(mar=c(5,4,0,0)+.1)
  plot(NA,type="n",xlim=c(1e-7,1),ylim=c(0,5),log="x",axes=FALSE,ylab="",xlab="marginal frequency")
  axis(1)
  stops <- 10^seq(-7,0,length.out=100)
  histo <- hist(aaMarginal$freq,breaks=c(0,stops),plot=FALSE)
  bars <- histo$counts/max(histo$counts)
  rect(stops[-100],0,stops[-1],1,col=cmap(log10(stops[-100])),border=NA)
  rect(stops[-100],1,stops[-1],1+4*bars[-100],col="gray",border=NA)
  abline(v=wmThreshold,lty="dashed")
  par(op)
}


#' Output coverage curves showing the progression of how many variants are observed at 
#' increasing marginal frequency thresholds
#'
#' @param coverageCurves dataframe with cumulative distribution of variants observed 
#' at different marginal frequency thresholds
#' @param thresholds colour transition thresholds
#' @param ri region number
#' @param wmThreshold well-measuredness threshold
#' @param tagger list that contains functions for adding tags to page margins
#' 
#' @return NULL
#' @export
#' 
plotCoverageCurve <- function(coverageCurves, thresholds, ri, wmThreshold, tagger) {
  plotcolors <- c(
    `codon changes`="black",`SNVs`="gray50",
    `AA changes`="firebrick3",`SNV-reachable AA changes`="firebrick2"
  )
  linetypes <- rep(c("solid","dashed"),2)
  plot(NA,type="n",
       log="x",xlim=range(thresholds),ylim=c(0,1),
       xlab="Marginal read frequency cutoff",
       ylab="Fraction of possible variants measured",
       main=paste0("Region #",ri)
  )
  for (i in 1:4) {
    lines(coverageCurves[,c(1,i+1)],col=plotcolors[[i]],lty=linetypes[[i]],lwd=2)
  }
  grid()
  abline(v=wmThreshold,lty="dotted",col="gray50")
  text(wmThreshold,1,sprintf("legacy cutoff (%.0e)",wmThreshold),col="gray50",pos=4)
  legend("bottomleft",names(plotcolors),col=plotcolors,lty=linetypes,lwd=2)
  tagger$cycle()
}


#' Determine complexity of all varinats in allCounts; look at marginal frequency 
#' against how many unique contexts the variant is seen in
#'
#' @param allCounts all counts table
#' @param simplifiedMarginal simplified marginal frequency table 
#' 
#' @return table with "complexity" vs marginal frequency
#' @export
#' 
analyzeComplexity <- function(allCounts, simplifiedMarginal) {
  
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
  
  return(cplxPlot)
}


#' Output complexity plot with the marginal frequency of variants plotted 
#' against the number of unique genetic contexts in which they were observed
#'
#' @param info list containing details on condition (nscond), timepoint (tp), 
#' pdftag and outdir
#' @param cplxPlot table with "complexity" vs marginal frequency
#' 
#' @return NULL
#' @export
#' 
drawComplexity <- function(info, cplxPlot) {
  #draw the plot for each tile
  pdf(paste0(info$outDir,info$nsCond,"_t",info$tp,"_complexity.pdf"),8.5,11)
  opar <- par(mfrow=c(3,2),oma=c(2,2,2,2))
  tagger <- pdftagger(paste(info$pdftag,"; condition:",info$nsCond,
                            "; time point:",info$tp),cpp=6)
  invisible(tapply(1:nrow(cplxPlot),cplxPlot[,"tile"],function(is) {
    if (length(is) > 1) {
      tile <- unique(cplxPlot[is,3])
      if (!all(cplxPlot[is,2] == 0,na.rm=TRUE)) {
        plot(cplxPlot[is,1:2],log="xy",
             xlab="#unique contexts in tile",ylab="marginal frequency",
             main=paste0("Tile #",tile)
        )
        tagger$cycle()
      }
    }
  }))
  par(opar)
  invisible(dev.off())
}


#' Return the sums of marginal frequencies for different types of mutations in 
#' each tile
#'
#' @param tile tile number being analyzed
#' @param simplifiedMarginal simplified marginal frequency table
#' @param freq vector of frequencies for each mutation in simplifiedMarginal
#' @param fromaa starting amino acid
#' @param toaa ending amino acid
#' @param to length of codon change ****
#'
#' @return list with the sums of marginal frequencies for different mutations types
#' @export
#' 
mut_breakdown <- function(tile, simplifiedMarginal, freq, fromaa, toaa, to) {
  if (!any(simplifiedMarginal$tile == tile)) {
    return(setNames(rep(NA,5),c("fs","indel","stop","syn","mis")))
  }
  marginalSubset <- simplifiedMarginal[which(simplifiedMarginal$tile == tile),]
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
}


#' Output variant-type breakdown for each tile
#'
#' @param info list containing condition details (nscond), timepoint (tp), 
#' pdftag and outdir
#' @param mutbreakdown table containing the sums of marginal frequencies 
#' for different mutation types
#' 
#' @return NULL
#' @export
#' 
plotMutbreakdown <- function(info, mutbreakdown) {
  plotcols <- c(
    frameshift="firebrick3",indel="orange",nonsense="gold",
    synonymous="steelblue3",missense="darkolivegreen3"
  )
  pdf(paste0(info$outDir,info$nsCond,"_t",info$tp,"_mutationtypes.pdf"),8.5,11)
  opar <- par(oma=c(2,2,2,2))
  tagger <- pdftagger(paste(info$pdftag,"; condition:",info$nsCond,"; time point:",info$tp),cpp=2)
  barplot(mutbreakdown,col=plotcols,border=NA,horiz=TRUE,
          xlab="sum of marginal frequencies",ylab="Tile",main="Mutation types"
  )
  legend("right",names(plotcols),fill=plotcols)
  tagger$cycle()
  par(opar)
  invisible(dev.off())
}


#' Helper object for adding tags to pdf page margins
#'
#' @param tag str, description tag to be added
#' @param cpp number of cycles per page
#' 
#' @return list that contains functions for adding tags to page margins
#' @export
#' 
pdftagger <- function(tag, cpp=1) {
  .cyclesPerPage <- cpp
  .cycle <- 1
  cycle <- function() {
    #add one pdf tag per page
    if (.cyclesPerPage==1 || (.cycle%%.cyclesPerPage==1)) {
      mtext(tag,outer=TRUE,line=0,cex=0.5)
    }
    .cycle <<- .cycle+1
  }
  reInit <- function(cpp=1) {
    .cyclesPerPage <<- cpp
    .cycle <- 1
  }
  list(cycle=cycle,reInit=reInit)
}

