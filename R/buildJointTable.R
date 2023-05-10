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
#' @param outDir output directory, defaults to input directory.
#' @param paramFile input parameter file. defaults to <dataDir>/parameters.json
#' @param mc.cores the number of CPU cores to use in parallel
#' @param srOverride the single-replicate override flag
#' @param covOverride override the requirement for coverage files. For backwards-
#'   compatibility with older versions prior to v0.7
#' @return NULL. Results are written to file.
#' @export
buildJointTable <- function(dataDir,inDir=NA,outDir=NA,
                            paramFile=paste0(dataDir,"parameters.json"),
                            mc.cores=6,srOverride=FALSE,covOverride=FALSE) {

  op <- options(stringsAsFactors=FALSE)

  library(yogitools)
  library(hash)
  library(hgvsParseR)
  library(parallel)
  library(pbmcapply)

  #make sure data and out dir exist and ends with a "/"
  if (!grepl("/$",dataDir)) {
    dataDir <- paste0(dataDir,"/")
  }
  if (!dir.exists(dataDir)) {
    #we don't use the logger here, assuming that whichever script wraps our function
    #catches the exception and writes to the logger (or uses the log error handler)
    stop("Workspace folder ",dataDir," does not exist!")
  }
  
  
  #Read parameters
  if (!canRead(paramFile)) {
    stop("Unable to read parameter file!")
  }
  logInfo("Reading parameters from",normalizePath(paramFile))
  params <- withCallingHandlers(
    parseParameters(paramFile,srOverride=srOverride),
    warning=function(w)logWarn(conditionMessage(w))
  )
  
  #configure input and output folders
  if (is.na(inDir)) {
    latest <- latestSubDir(parentDir=dataDir,pattern="_mut_call$|mut_count$")
    inDir <- latest[["dir"]]
  } else {
    if (!dir.exists(inDir)) {
      stop("Input folder ",inDir," does not exist!")
    }
  }
  if (!grepl("/$",inDir)) {
    inDir <- paste0(inDir,"/")
  }
  if (is.na(outDir)) {
    outDir <- inDir
  } else if (!dir.exists(outDir)) {
    dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
  }

  if (!grepl("/$",outDir)) {
    outDir <- paste0(outDir,"/")
  }
  
  logInfo("Using input directory",inDir,"and output directory",outDir)
  
  
  ######################
  # BUILD SAMPLE TABLE #
  ######################
  
  #prepare a table of all samples
  logInfo("Accounting for all input count files")
  sampleTable <- params$samples
  
  #find count table files and assign each to its respective sample
  countfiles <- list.files(inDir,pattern="counts_sample_.+\\.csv$",full.names=TRUE)
  if (length(countfiles)==0) {
    error("No count files found in ",inDir,"! (Are they named correctly?)")
  }
  filesamples <- extract.groups(countfiles,"counts_sample_([^/]+)\\.csv$")[,1]
  sampleTable$countfile <- sapply(as.character(sampleTable$`Sample ID`), function(sid) {
    i <- which(filesamples == sid)
    if (length(i) == 0 || is.null(i)) {
      NA
    } else {
      countfiles[[i]]
    }
  })

  #if any samples/files are unaccounted for, throw an error!
  if (any(is.na(sampleTable$countfile))) {
    missingFiles <- with(sampleTable,`Sample ID`[which(is.na(countfile))])
    stop("Missing count files for the following samples: ",paste(missingFiles,collapse=", "))
  }

  if (!covOverride) {
    #find coverage files and assign each to its respective sample
    covfiles <- list.files(inDir,pattern="coverage_.+\\.csv$",full.names=TRUE)
    if (length(covfiles)==0) {
      error("No coverage files found in ",inDir,"! (Are they named correctly?)")
    }
    filesamples <- extract.groups(covfiles,"coverage_([^/]+)\\.csv$")[,1]
    sampleTable$coveragefile <- sapply(as.character(sampleTable$`Sample ID`), function(sid) {
      i <- which(filesamples == sid)
      if (length(i) == 0 || is.null(i)) {
        NA
      } else {
        covfiles[[i]]
      }
    })
    
    #if any samples/files are unaccounted for, throw an error!
    if (any(is.na(sampleTable$coveragefile))) {
      missingFiles <- with(sampleTable,`Sample ID`[which(is.na(coveragefile))])
      stop("Missing coverage data files for the following samples: ",paste(missingFiles,collapse=", "))
    }
  }

  #####################################
  # Read and prepare all input tables #
  #####################################

  logInfo("Reading count data")
  allCounts <- lapply(sampleTable$countfile,function(cfile) {
    tryCatch(
      parseCountFile(cfile),
      error=function(e) stop("Error reading file ",cfile," : ",conditionMessage(e))
    )
  })
  
  if (!covOverride) {
    logInfo("Reading depth data")
    #parse coverage files to calculate the number of rejected variant calls per position
    allRejects <- lapply(1:nrow(sampleTable), function(i) {
      parseCoverageFile(sampleTable[i,"coveragefile"],params,sampleTable[i,"Tile ID"])
    })
    #the number of reads minus the number of rejections gives us the positional depth
    positionalDepth <- lapply(1:nrow(sampleTable), function(i) {
      as.integer(attr(allCounts[[i]],"depth")) - allRejects[[i]]
    })
    #re-organize the positional depth by condition timepoint and replicate
    condRID <- with(sampleTable,sprintf("%s.t%s.rep%d",Condition,`Time point`,Replicate))
    allPos <- sort(as.integer(Reduce(union,lapply(positionalDepth,names))))
    positionalDepth <- do.call(rbind,tapply(1:nrow(sampleTable),condRID,function(is) {
      pds <- do.call(c,positionalDepth[is])
      if (any(duplicated(names(pds)))) {
        stop("Overlapping tiles in coverage data!")
      }
      # pds[order(as.integer(names(pds)))]
      setNames(pds[as.character(allPos)],allPos)
    }))
    if (any(is.na(positionalDepth))) {
      positionalDepth[which(is.na(positionalDepth))] <- 0
    }
    
      
  }

  #tabulate raw sequencing depths for each sample
  sampleTable$alignedreads <- sapply(allCounts,function(counts) as.integer(attr(counts,"depth"))) 
  sampleTable$wtreads <- sapply(allCounts,function(counts) as.integer(attr(counts,"wtpairs"))) 
  sampleTable$mutreads <- sapply(allCounts,function(counts) as.integer(attr(counts,"mutpairs"))) 
  sampleTable$depth <- sampleTable$wtreads + sampleTable$mutreads
  #and save the sample table for future reference
  logInfo("Exporting sequencing depth information.")
  outfile <- paste0(outDir,"/sampleDepths.csv")
  write.csv(sampleTable,outfile,row.names=FALSE)
  
  if (!covOverride) {
    outfile <- paste0(outDir,"/positionalDepths.csv")
    write.csv(as.data.frame(positionalDepth),outfile)
  }

  #Double-check that we have coverage for all samples
  if (any(is.na(sampleTable$depth)) || any(sampleTable$depth <= 0)) {
    uncovered <- with(sampleTable,c(which(is.na(depth)),which(depth <= 0)))
    noreads <- with(sampleTable,c(which(is.na(alignedreads)),which(alignedreads <= 0)))
    uncovered <- setdiff(uncovered,noreads)
    if (length(noreads > 0)) {
      stop("No aligned reads found for sample(s): ",
        paste(sampleTable[noreads,"Sample ID"],collapse=", "),
        "\nPlease check your tileseqMut output for failed alignments."
      )
    }
    if (length(uncovered) > 0) {
      stop("No variants were called for sample(s): ",
        sampleTable[uncovered,"Sample ID"],
        "\nPlease check your parameter sheet and your tileseqMut output."
      )
    }
  }
  
  #helper function to calculate the combined effective depth
  #for co-occurring variants
  combineDepths <- function(depths,rawDepth) {
    if (length(depths)==1) {
      return(depths)
    } else {
      prod(depths[[1]],depths[-1]/rawDepth)
    }
  }
  
  #helper function to extract positions from HGVS strings
  #only positions relevant to effective depth calculation!!
  extractHGVSPositions <- function(hgvss) {
    #break multi-variants into components and iterate
    pbmclapply(strsplit(hgvss,";"), function(chunksList) {
      #process each chunk
      do.call(c,lapply(chunksList, function(chunk) {
        #check if there is a positional range or just a single position present
        if (grepl("_",chunk)) {
          #for a range, extract the start and stop
          startStop <- as.integer(yogitools::extract.groups(chunk,"(\\d+)_(\\d+)")[1,])
          #for delins use the start + length of insertion
          if (grepl("delins",chunk)) {
            l <- nchar(gsub(".+delins|\\]","",chunk))
            startStop[[1]]:(startStop[[1]]+(l-1))
          } else if (grepl("del",chunk)) {
            #for deletions, only the start counts
            startStop[[1]]
          } else {
            #for insertions, only the "end" indicates the actual position of the base
            startStop[[2]]
          }
        } else {
          #for a single position, just extract it
          as.integer(gsub("\\D+","",chunk))
        }
      }))
    },mc.cores=mc.cores)
  }
  

  #find the union of all variants
  allVars <- Reduce(union,lapply(allCounts,function(x)x$HGVS))
  logInfo(sprintf(
    "Translating %d unique variant calls to protein level. This may take some time...", 
    length(allVars)
  ))
  #translate them to amino acid level
  builder <- new.hgvs.builder.p(aacode=3)
  cbuilder <- new.hgvs.builder.c()
  # transTable <- as.df(pbmclapply(allVars, translateHGVS, params, builder, cbuilder, mc.cores=mc.cores))
  cdsSeq <- params$template$cdsSeq
  transTable <- as.df(pbmclapply(allVars, function(mut) {
    tryCatch({
      hgvsParseR::translateHGVS(mut, cdsSeq, builder, cbuilder)
    },error=function(e){
      c(
        hgvsc=mut,hgvsp=as.character(e),codonChanges=NA,codonHGVS=NA,
        aaChanges=NA,aaChangeHGVS=NA
      )
    })
  },mc.cores=mc.cores))

  #save result for faster debugging
  outfile <- paste0(outDir,"/transTable.rData")
  save(transTable,file=outfile)
  
  #Check for errors
  if (any(is.na(transTable[,3]))) {
    for (i in which(is.na(transTable[,3]))) {
      logWarn("Translation for variant ",allVars[[i]],"failed: ",transTable[i,2])
    }
    stop("Translations for ",sum(is.na(transTable[,3]))," variants failed! See log for details.")
  }
  
  # save(transTable,file="transTable.rda")

  ################
  # Merge tables #
  ################

  #helper function to sort unique numbers by abundance
  mostUnique <- function(xs) {
    tbl <- table(xs)
    as.numeric(names(tbl)[order(tbl,decreasing=TRUE)])
  }
  
  #helper function to find the most applicable tile assignments for set of positions
  ncpos2tile <- function(posGroups) {
    #and find tile assignments for all variants
    bestTiles <- pbmclapply(posGroups, function(poss) {
      rows <- sapply(poss,function(pos) {
        i <- which(
          params$tiles[,"Start NC in CDS"] <= pos & 
            params$tiles[,"End NC in CDS"] >= pos
        )
        if (length(i)==0) NA else i
      }) 
      mostUnique(params$tiles[rows,"Tile Number"])
    },mc.cores=mc.cores)
    
    #check for problematic cases and complain if necessary
    if (any(sapply(bestTiles,length)>1)) {
      culprits <- which(sapply(bestTiles,length)>1)
      for (culprit in culprits) {
        logWarn("The following variant positions cross tile-boundaries: ",
                names(posGroups)[[culprit]]
                # ,"\n(",sapply(posGroups[[culprit]],paste,collapse=","),")"
        )
      }
    }
    #reduce to just the most likely tile per variant
    sapply(bestTiles, function(ts) {
      if (length(ts) > 0) {
        ts[[1]]
      } else {
        NA
      }
    })
  }

  nuc2tile <- function(pos) {
    rows <- sapply(pos,function(pos) {
      i <- which(
        params$tiles[,"Start NC in CDS"] <= pos & 
        params$tiles[,"End NC in CDS"] >= pos
      )
      if (length(i)==0) NA else i
    }) 
    params$tiles[rows,"Tile Number"]
  }
  
  #Build a lookup table for the cleaned HGVS strings
  logInfo("Merging equivalent variants...")
  cleanHGVS <- hash::hash(allVars,transTable[,"hgvsc"])
  transTableClean <- as.df(tapply(1:nrow(transTable),transTable$hgvsc,function(is) {
    transTable[is[[1]],]
  }))

  # save(transTableClean,file="transTableClean.rda")


  #after cleanup we don't need the original table anymore, so we can remove it
  #to save RAM
  rm(transTable)
  
  logInfo("Indexing variant positions...")
  #extract positions for all variants
  allPositions <- extractHGVSPositions(transTableClean[,"hgvsc"])
  names(allPositions) <- transTableClean[,"hgvsc"]
  
  logInfo("Indexing tile assignments...")
  #and their most applicable tiles
  allTiles <- ncpos2tile(allPositions)
  
  #index sample table by condition-timepoint-replicate
  sampleTable$condRID <- sapply(1:nrow(sampleTable),function(i)with(sampleTable[i,],
     sprintf("%s.t%s.rep%d",Condition,`Time point`,Replicate)
  ))
  
  logInfo("Merging tables and calculating frequencies...")
  
  #join the tables
  condTables <- tapply(1:nrow(sampleTable),sampleTable$condRID,function(isamples) {
    
    out <- data.frame(
      count=rep(0,nrow(transTableClean)),
      frequency=rep(0,nrow(transTableClean)),
      effectiveDepth=rep(0,nrow(transTableClean)),
      row.names=transTableClean[,"hgvsc"]
    )
    
    for (isample in isamples) {
      currTile <- sampleTable[isample,"Tile ID"]
      counts <- allCounts[[isample]]
      rawDepth <- as.integer(attr(counts,"depth"))
      crid <- sampleTable[isample,"condRID"]
      
      applVars <- which(allTiles==currTile)
      if (covOverride) {
        #in case of old data without positional coverage we use the old method
        depth <- as.integer(attr(counts,"wtpairs")) + as.integer(attr(counts,"mutpairs"))
        out[applVars,"effectiveDepth"] <- depth
      } else {
        #look up the positional depths for this condition
        localDepth <- positionalDepth[crid,]
        #pull up the list of positions for these variants
        posList <- allPositions[applVars]
        #for multi-mutants, we must combine the depths of the applicable positions
        combinedDepth <- sapply(posList, function(poss) {
          valid <- which(nuc2tile(poss)==currTile)
          if (length(valid) > 0) {
            ds <- localDepth[as.character(poss[valid])]
            combineDepths(ds,rawDepth)
          } else {
            logWarn("Invalid variant outside of tile!")
            NA
          }
        })
        out[applVars,"effectiveDepth"] <- combinedDepth
      }
      
      #condense counts by equivalent HGVS
      if (nrow(counts) > 0) {
        cCounts <- tapply(counts[,"count"],hash::values(cleanHGVS,counts$HGVS),sum,na.rm=TRUE)
        rows <- names(cCounts)
        # out[rows,"count"] <- out[rows,"count"] + counts[,"count"]
        out[rows,"count"] <- out[rows,"count"] + cCounts
      }
    }
    
    out$frequency <- out$count/out$effectiveDepth
    
    return(out)
    
  })
  
  #add back HGVS strings and translations in the final column binding step
  jointTable <- cbind(
    # HGVS=allVars,
    # HGVS_pro=values(trIdx,keys=allVars),
    transTableClean,
    do.call(cbind,condTables)
  )
  #after joining we don't need the individual tables anymore, so we can remove them
  #to save RAM
  rm(condTables)
  

  ########################
  # WRITE OUTPUT TO FILE #
  ########################
  logInfo("Writing results to file.")
  outfile <- paste0(outDir,"/allCounts.csv")
  cat("# COMBINATORY VARIANT COUNTS #",
      "\n# project name:", params$project,
      "\n# gene name:",params$template$geneName,
      "\n# tileseqMave version:",as.character(packageVersion("tileseqMave")),
      "\n# parameter sheet:",normalizePath(paramFile),"\n",
      file=outfile
  )
  suppressWarnings(
    write.table(jointTable,outfile,sep=",",append=TRUE,row.names=FALSE,qmethod="double")
  )
  
  ##################################
  # CALCULATE MARGINAL FREQUENCIES #
  ##################################

  logInfo("Calculating marginal frequencies...")

  #split codon change HGVSs by codon
  codonChangeHGVSs <- strsplit(gsub("c\\.|c\\.\\[|\\]","",transTableClean$codonHGVS),";")
  #and attach the proper prefix
  codonChangeHGVSs <- lapply(codonChangeHGVSs, function(x) paste0("c.",x))
  #do the same for AA change HGVSs
  aaChangeHGVSs <- strsplit(gsub("p\\.|p\\.\\[|\\]","",transTableClean$aaChangeHGVS),";")
  aaChangeHGVSs <- lapply(aaChangeHGVSs, function(x) paste0("p.",x))
  #now split the native codon changes and AA changes
  codonChangeStrs <- strsplit(transTableClean$codonChanges,"\\|")
  aaChangeStrs <- strsplit(transTableClean$aaChanges,"\\|")
  
  #helper function to merge two HGVS terms in-cis
  cisMerge <- function(hs) {
    if (length(hs) < 2) return(hs)
    paste0("c.[",paste(substr(hs,3,nchar(hs)),collapse=";"),"]")
  }
  
  #build codon change index
  ccIdx <- hash()
  ccStrIdx <- hash()
  aaIdx <- hash()
  aaStrIdx <- hash()
  for (i in 1:length(codonChangeHGVSs)) {
    # js <- min(length(codonChangeHGVSs[[i]]), length(aaChangeHGVSs[[i]]))
    #js is the number of aa changes in this variant
    js <- length(aaChangeHGVSs[[i]])
    #in case of frameshifts, all other mutations get ignored, so we need to figure out 
    #which one triggered the frameshift(s). In that case, we have fewer aa changes than codon changes listed
    if (length(codonChangeHGVSs[[i]]) != js) {
      #extract positions from aaChanges and  codon changes and find the match
      aapos <- as.integer(extract.groups(aaChangeStrs[[i]],"(\\d+)")[,1])
      ncpos <- as.integer(extract.groups(codonChangeHGVSs[[i]],"(\\d+)")[,1])
      ccpos <- floor((ncpos-1)/3+1)
      #ks are the codon change HGVSs that match the positions for each j
      ks <- lapply(aapos,function(ap) {
        ds <- abs(ccpos-ap)
        matches <- which(ds==min(ds))
        matches[order(ncpos[matches])]
      })
    } else {
      #ks are the codon change HGVSs that match the positions for each j
      ks <- as.list(1:js)
    }
    #now iterate over aa changes and index them
    for (j in 1:js) {
      cc <- cisMerge(codonChangeHGVSs[[i]][ks[[j]]])
      if (has.key(cc,ccIdx)) {
        ccIdx[[cc]] <- c(ccIdx[[cc]],i)
      } else {
        ccIdx[[cc]] <- i
        aaIdx[[cc]] <- aaChangeHGVSs[[i]][[j]]
        ccStrIdx[[cc]] <- codonChangeStrs[[i]][[j]]
        aaStrIdx[[cc]] <- aaChangeStrs[[i]][[j]]
      }
    }
  }

  #build marginals using index
  depthCols <- grep("effectiveDepth",colnames(jointTable))
  marginalCCs <- keys(ccIdx)
  marginalCounts <- as.df(pbmclapply(marginalCCs, function(cc) {
    # jointTable[ccIdx[[cc]],-(1:6)]
    c(
      list(hgvsc=cc,hgvsp=aaIdx[[cc]]),
      codonChange=ccStrIdx[[cc]],aaChange=aaStrIdx[[cc]],
      # colSums(jointTable[ccIdx[[cc]],-(1:6)],na.rm=TRUE)
      setNames(lapply(7:ncol(jointTable),function(coli) {
        if (coli %in% depthCols) {#depths are averaged
          mean(jointTable[ccIdx[[cc]],coli],na.rm=TRUE)
        } else {#counts and frequencies are summed
          sum(jointTable[ccIdx[[cc]],coli],na.rm=TRUE)
        }
      }),colnames(jointTable)[-(1:6)])
    )
  },mc.cores=mc.cores))
  
  #if positional coverage data is available, we need to correct the frequencies
  if (!covOverride) {
    logInfo("Indexing marginal variant positions...")
    margPositions <- extractHGVSPositions(marginalCCs)
    names(margPositions) <- marginalCCs
    logInfo("Indexing marginal variant tile assignments...")
    margTiles <- ncpos2tile(margPositions)
    
    logInfo("Correcting marginal frequencies for position-specific depth...")
    for (rowi in 1:nrow(marginalCounts)) {
      tili <- margTiles[[rowi]]
      #pull up the list of positions for these variants
      poss <- margPositions[[rowi]]
      valid <- which(nuc2tile(poss)==tili)
      poss <- poss[valid]
      #process across all condition-timepoint-replicate groups
      if (!is.na(tili) && length(poss) > 0 && !any(is.na(poss))) {
        for (crid in unique(sampleTable$condRID)) {
          #pull up the relevant raw-depth
          rawDepth <- with(sampleTable,alignedreads[condRID==crid & `Tile ID`==tili])
          #pull up the relevant positional depths 
          relevPos <- intersect(colnames(positionalDepth),as.character(poss))
          ds <- positionalDepth[crid,relevPos]
          if (length(ds) > 0) {#skip invalid cases
            #combine them using the formula
            effDepth <- combineDepths(ds,rawDepth)
            #and overwrite the old results
            cnt <- marginalCounts[rowi,paste0(crid,".count")]
            marginalCounts[rowi,paste0(crid,".effectiveDepth")] <- effDepth
            marginalCounts[rowi,paste0(crid,".frequency")] <- cnt/effDepth
          }
        }
      }
    }
    
  }

  logInfo("Writing results to file.")
  outfile <- paste0(outDir,"/marginalCounts.csv")
  cat("# MARGINAL VARIANT COUNTS #",
      "\n# project name:", params$project,
      "\n# gene name:",params$template$geneName,
      "\n# tileseqMave version:",as.character(packageVersion("tileseqMave")),
      "\n# parameter sheet:",normalizePath(paramFile),"\n",
      file=outfile
  )
  suppressWarnings(
    write.table(marginalCounts,outfile,sep=",",append=TRUE,row.names=FALSE,qmethod="double")
  )

  logInfo("Done.")

  options(op)
  return(NULL)
}
