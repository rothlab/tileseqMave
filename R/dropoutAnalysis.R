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

#' run variant dropout troubleshooter
#' 
#' @param dataDir working data directory
#' @param countDir input directory for counts, defaults to subdirectory with latest timestamp ending in _mut_count.
#' @param outDir output directory, defaults to name of input directory with _QC tag attached.
#' @param paramFile input parameter file. defaults to <dataDir>/parameters.json
#' @param srOverride single replicate override (use with caution)
#' @return NULL. Results are written to file.
#' @export
dropoutAnalysis <- function(dataDir,countDir=NA, #scoreDir=NA, 
                        outDir=NA, 
                        paramFile=paste0(dataDir,"parameters.json"),
                        srOverride=FALSE,tile=NA) {

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
  
  # if (is.na(tile)) {
  #   #if no tile is defined, default to first tile
  #   tile <- params$tiles[1,1]
  # }
  
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
  
  # #if no score directory was provided
  # if (is.na(scoreDir)) {
  #   #try to guess its name
  #   if (grepl("_mut_count/$",countDir)) {
  #     scoreDir <- sub("_mut_count/$","_scores/",countDir)
  #   } else {
  #     scoreDir <- sub("/$","_scores/",countDir)
  #   }
  #   if (!dir.exists(scoreDir)) {
  #     stop("No matching score directory found for ",countDir,"!\n",
  #          "(Expecting ",scoreDir,")\n",
  #          "Use --scores option to define explicitly."
  #     )
  #   }
  # } else {
  #   if (!dir.exists(scoreDir)) {
  #     stop("Score directory ",scoreDir," does not exist!")
  #   }
  # }
  # #make sure it ends in "/"
  # if (!grepl("/$",scoreDir)) {
  #   scoreDir <- paste0(scoreDir,"/")
  # }
  
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
          # "score directory", scoreDir,
          "and output directory",outDir
  )
  #create PDF tag
  # pdftag <- with(params,sprintf("%s (%s): %s%s",project,template$geneName,runLabel,timeStamp))
  # params$pdftagbase <- pdftag
  
  
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
  
  #filter down to selected tile
  if (!is.na(tile) && is.numeric(tile)) {
    positions <- as.integer(extract.groups(marginalCounts$codonChange,"(\\d+)")[,1])
    tiles <- params$pos2tile(positions)
    marginalCounts <- marginalCounts[which(tiles==tile),]
  }
  
  #helper function to compare dropout sets
  compareSets <- function(set1,set2,title,n=nrow(marginalCounts)) {
    l1 <- length(set1)
    l2 <- length(set2)
    common <- length(intersect(set1,set2))
    nullDistr <- replicate(1000,{
      length(intersect(sample(n,l1),sample(n,l2)))
    })
    mNull <- mean(nullDistr)
    q95 <- quantile(nullDistr,c(0.025,0.975))
    sprintf("%s: %.0f from %.0f/%.0f (expected %.01f [%.01f;%.01f])",title,common,l1,l2,mNull,q95[[1]],q95[[2]])
  }
  
  #iterate over conditions
  for (sCond in getSelects(params)) {
    
    #iterate over timepoints
    for (tp in params$timepoints$`Time point name`) {
      
      logInfo("Processing",sCond,"time",tp)
      
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
      
      smallestSelect <- unique(sort(na.omit(marginalCounts[,repMatrix["select",1]])))[[2]]
      pseudoCount <- 10^floor(log10(smallestSelect))
      
      # plotCol <- sapply(1+(
      #   log10(marginalCounts[,repMatrix["select",1]]+pseudoCount) < log10(marginalCounts[,repMatrix["nonselect",1]])/4-4.5
      # ),colAlpha,0.5)
      # 
      quickCounts <- log10(marginalCounts[,as.vector(repMatrix[1:2,])]+pseudoCount)
      maxCount <- apply(quickCounts,1,max)
      
      plotCol <- sapply(1+(
        quickCounts[,repMatrix["nonselect",1]] < (maxCount-1)
      ),colAlpha,0.2)
      
      outfile <- paste0(outDir,sCond,"_t",tp,"_dropoutAll.pdf")
      pdf(outfile,10,10)
      opar <- par(mfrow=c(2,2),mar=c(5,4,1,1))
      plot(
        jitter(maxCount,amount=0.2),
        jitter(quickCounts[,repMatrix["nonselect",1]],amount=0.2),
        pch=20,col=plotCol,
        xlab="max condition",ylab="nonselect rep 1"
      )
      plot(
        jitter(maxCount,amount=0.2),
        jitter(quickCounts[,repMatrix["nonselect",2]],amount=0.2),
        pch=20,col=plotCol,
        xlab="max condition",ylab="nonselect rep 2"
      )
      plot(
        jitter(maxCount,amount=0.2),
        jitter(quickCounts[,repMatrix["select",1]],amount=0.2),
        pch=20,col=plotCol,
        xlab="max condition",ylab="select rep 1"
      )
      plot(
        jitter(maxCount,amount=0.2),
        jitter(quickCounts[,repMatrix["select",2]],amount=0.2),
        pch=20,col=plotCol,
        xlab="max condition",ylab="select rep 2"
      )
      par(opar)
      invisible(dev.off())
      # xs <- 10^seq(-7,-1)
      # lines(xs,xs*10)
      
      droppedS1 <- which(quickCounts[,repMatrix["select",1]] < maxCount-1)
      droppedNS1 <- which(quickCounts[,repMatrix["nonselect",1]] < maxCount-1)
      droppedS2 <- which(quickCounts[,repMatrix["select",2]] < maxCount-1)
      droppedNS2 <- which(quickCounts[,repMatrix["nonselect",2]] < maxCount-1)
      
      logInfo("DROP FROM ALL:")
      n <- sum(maxCount > -6)
      logInfo(compareSets(droppedS1,droppedNS1,"max->s1:max->ns1",n=n))
      logInfo(compareSets(droppedNS1,droppedNS2,"max->ns1:max->ns2",n=n))
      logInfo(compareSets(droppedS1,droppedS2,"max->s1:max->s2",n=n))
      logInfo(compareSets(droppedS1,droppedNS2,"max->s1:max->ns2",n=n))
      
      
      
      
      
      
      outfile <- paste0(outDir,sCond,"_t",tp,"_dropouts.pdf")
      pdf(outfile,10,10)
      plotCol <- sapply(1+(
        log10(marginalCounts[,repMatrix["nonselect",2]]+pseudoCount) < log10(marginalCounts[,repMatrix["nonselect",1]]+pseudoCount)-1
      ),colAlpha,0.2)
      opar <- par(mfrow=c(2,2),mar=c(5,4,1,1))
      plot(
        jitter(log10(marginalCounts[,repMatrix["nonselect",1]]+pseudoCount),amount=0.2),
        jitter(log10(marginalCounts[,repMatrix["nonselect",2]]+pseudoCount),amount=0.2),
        pch=20,col=plotCol,
        xlab="nonselect rep 1",ylab="nonselect rep 2"
      )
      plot(
        jitter(log10(marginalCounts[,repMatrix["nonselect",1]]+pseudoCount),amount=0.2),
        jitter(log10(marginalCounts[,repMatrix["select",1]]+pseudoCount),amount=0.2),
        pch=20,col=plotCol,
        xlab="nonselect rep 1",ylab="select rep 1"
      )
      plot(
        jitter(log10(marginalCounts[,repMatrix["select",1]]+pseudoCount),amount=0.2),
        jitter(log10(marginalCounts[,repMatrix["select",2]]+pseudoCount),amount=0.2),
        pch=20,col=plotCol,
        xlab="select rep 1",ylab="select rep 2"
      )
      plot(
        jitter(log10(marginalCounts[,repMatrix["nonselect",2]]+pseudoCount),amount=0.2),
        jitter(log10(marginalCounts[,repMatrix["select",1]]+pseudoCount),amount=0.2),
        pch=20,col=plotCol,
        xlab="nonselect rep 2",ylab="select rep 1"
      )
      par(opar)
      invisible(dev.off())
      # xs <- 10^seq(-7,-1)
      # lines(xs,xs*10)
      
      logInfo("DROP RELATIVE:")
      
      dropSel <- which(marginalCounts[,repMatrix["select",1]] < marginalCounts[,repMatrix["nonselect",1]]/10)
      riseSel <- which(marginalCounts[,repMatrix["select",1]] > marginalCounts[,repMatrix["nonselect",1]]*10)
      
      dropRep2 <- which(marginalCounts[,repMatrix["nonselect",2]] < marginalCounts[,repMatrix["nonselect",1]]/10)
      riseRep2 <- which(marginalCounts[,repMatrix["nonselect",2]] > marginalCounts[,repMatrix["nonselect",1]]*10)
      
      dropSel2 <- which(marginalCounts[,repMatrix["select",1]] < marginalCounts[,repMatrix["nonselect",2]]/10)
      riseSel2 <- which(marginalCounts[,repMatrix["select",1]] > marginalCounts[,repMatrix["nonselect",2]]*10)

      logInfo(compareSets(dropRep2,dropSel,"ns2<-ns1->s1"))
      logInfo(compareSets(riseRep2,dropSel2,"ns1<-ns2->s1"))
      logInfo(compareSets(riseSel,riseSel2,"ns1<-s1->ns2"))
      

      
      dropSelMirr <- which(marginalCounts[,repMatrix["select",2]] < marginalCounts[,repMatrix["nonselect",2]]/10)
      riseSelMirr <- which(marginalCounts[,repMatrix["select",2]] > marginalCounts[,repMatrix["nonselect",2]]*10)

      dropRepMirr <- which(marginalCounts[,repMatrix["select",2]] < marginalCounts[,repMatrix["select",1]]/10)
      riseRepMirr <- which(marginalCounts[,repMatrix["select",2]] > marginalCounts[,repMatrix["select",1]]*10)

      logInfo(compareSets(dropRep2,dropRepMirr,"ns1->ns2:s1->s2"))
      logInfo(compareSets(riseRep2,riseRepMirr,"ns1<-ns2:s1<-s2"))
      logInfo(compareSets(dropSel,dropSelMirr,"ns1->s1:ns2->s2"))
      logInfo(compareSets(riseSel,riseSelMirr,"ns1->s1:ns2->s2"))
      
    }
  }
  
  options(op)
}