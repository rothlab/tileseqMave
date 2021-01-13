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

#' perform scaling of logPhi into final scores
#' 
#' @param dataDir working data directory
#' @param scoreDir input directory for scores, defaults to subdirectory with latest timestamp ending in _scores.
#' @param outDir output directory, defaults to name of input directory with _QC tag attached.
#' @param paramFile input parameter file. defaults to <dataDir>/parameters.json
#' @param srOverride single replicate override
#' @param bnOverride bottleneck filter override
#' @return NULL. Results are written to file.
#' @export
scaleScores <- function(dataDir, scoreDir=NA, outDir=NA, 
                        paramFile=paste0(dataDir,"parameters.json"),
                        srOverride=FALSE, bnOverride=FALSE,autoPivot=FALSE) {
  
  
  op <- options(stringsAsFactors=FALSE)
  
  library(yogitools)
  # library(hgvsParseR)
  # library(pbmcapply)
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
  
  if (!autoPivot) {
    checkPivots(params)
  }
  
  #find counts folder
  if (is.na(scoreDir)) {
    latest <- latestSubDir(parentDir=dataDir,pattern="_scores$")
    scoreDir <- latest[["dir"]]
    timeStamp <- latest[["timeStamp"]]
    runLabel <- latest[["label"]]
  } else { #if custom input dir was provided
    #make sure it exists
    if (!dir.exists(scoreDir)) {
      stop("Input score folder ",scoreDir," does not exist!")
    }
    #try to extract a timestamp and label
    lt <- extract.groups(scoreDir,"([^/]+_)?(\\d{4}-\\d{2}-\\d{2}-\\d{2}-\\d{2}-\\d{2})")[1,]
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
  if (!grepl("/$",scoreDir)) {
    scoreDir <- paste0(scoreDir,"/")
  }
  
 
  #if not output directory was defined
  if (is.na(outDir)) {
    #use existing score directory as output directory
    outDir <- scoreDir
  } 
  #make sure it ends in "/"
  if (!grepl("/$",outDir)) {
    outDir <- paste0(outDir,"/")
  }
  #make sure outdir exists
  dir.create(outDir,recursive=TRUE,showWarnings=FALSE)
  
  logInfo("Using score directory", scoreDir,
          "and output directory",outDir
  )
  
  #iterate over conditions
  for (sCond in getSelects(params)) {
    
    #iterate over timepoints
    for (tp in params$timepoints$`Time point name`) {
      
      logInfo("Processing condition",sCond, "; time",tp)
      
      #load enrichment table for this condition
      inFile <- paste0(scoreDir,"/",sCond,"_t",tp,"_enrichment.csv")
      if (!file.exists(inFile)) {
        logWarn("No score file found! Skipping...")
        next
      }
      enrichments <- read.csv(inFile,comment.char="#")
      rownames(enrichments) <- enrichments$hgvsc
      enrichments$position <- as.integer(extract.groups(enrichments$codonChange,"(\\d+)"))
      enrichments$region <- params$pos2reg(enrichments$position)
      enrichments$tile <- params$pos2tile(enrichments$position)
      
      regions <- params$regions[,"Region Number"]
      scoreTable <- do.call(rbind,lapply(regions, function(region) {
        
        logInfo("Processing region",region)
        
        msc <-  enrichments[which(enrichments$region == region),]
        
        # Scaling to synonymous and nonsense medians -----------------------------
        logInfo("Scaling...")
        #check if overrides were provided for the distribution modes
        normOvr <- getPivotOverrides(params,sCond,tp,region)
        #and apply
        msc <- cbind(msc,enactScale(msc,msc$aaChange,params$scoring$sdThreshold,normOvr))
        
        # flooring and stderr calculation only where more than one replicate exists ----
        if (params$numReplicates[[sCond]] > 1) {
          #Flooring and SD adjustment
          msc <- cbind(msc,flooring(msc,params))
          
          #add stderr
          msc$se.floored <- msc$sd.floored/sqrt(params$numReplicates[[sCond]])
        }
        
        return(msc)
        
      }))
      
      #delete unnecessary columns
      scoreTable$position <- NULL
      scoreTable$region <- NULL
      scoreTable$tile <- NULL
      
      # export to file --------------------------------------------------------
      timeStamp <- if (exists("timeStamp")) timeStamp else "?"
      paramHeader <- paste0(
        "# project name: ", params$project,"\n",
        "# gene name: ",params$template$geneName,"\n",
        "# tileseqMave version: ",packageVersion("tileseqMave"),"\n",
        "# variant call date: ",timeStamp,"\n",
        "# enrichment processing date: ",Sys.time(),"\n",
        "# condition: ",sCond,"\n",
        "# time point: ",tp,"\n",
        "# parameter sheet: ",normalizePath(paramFile),"\n",
        "# count threshold: ",params$scoring$countThreshold,"\n",
        "# wt filter quantile: ",params$scoring$wtQuantile,"\n",
        "# pseudo-replicates: ",params$scoring$pseudo.n,"\n",
        "# syn/non sd threshold: ",params$scoring$sdThreshold,"\n",
        "# cv deviation threshold: ",params$scoring$cvDeviation,"\n",
        "# assay direction: ",params$assay[["selection"]],"\n"
      )
      
      
      outFile <- paste0(outDir,sCond,"_t",tp,"_complete.csv")
      cat("# COMPLETE UNFILTERED SCORES #\n",paramHeader,file=outFile,sep="")
      suppressWarnings(
        write.table(scoreTable,outFile,sep=",",append=TRUE,row.names=FALSE,qmethod="double")
      )
      
      #configure filter level for export
      exportFilter <- if (bnOverride) {
        sapply(scoreTable$filter,function(f) is.na(f) || grepl(f=="bottleneck"))
      } else {
        is.na(scoreTable$filter)
      }
      
      #simplified (MaveDB-compatible) format
      simpleCols <- if (srOverride) {
        #when single-replicate override is enabled, no flooring is available, 
        #so we default to unfloored values
        c("hgvsc","hgvsp","score","score.sd","score.sd")
      } else {
        c("hgvsc","hgvsp","score.floored","sd.floored","se.floored")
      }
      simpleTable <- scoreTable[
        exportFilter,
        simpleCols
      ]
      colnames(simpleTable) <- c("hgvs_nt","hgvs_pro","score","sd","se")
      
      #export to file
      logInfo("Writing simplified table to file.")
      outFile <- paste0(outDir,sCond,"_t",tp,"_simple.csv")
      cat("# NUCLEOTIDE-LEVEL SCORES #\n",paramHeader,file=outFile,sep="")
      suppressWarnings(
        write.table(simpleTable,outFile,sep=",",append=TRUE,row.names=FALSE,qmethod="double")
      )
      
      #collapse by amino acid change
      logInfo("Collapsing amino acid changes...")
      flooredAA <- collapseByAA(scoreTable,params,sCond,"score.floored","sd.floored",srOverride,bnOverride)
      unflooredAA <- collapseByAA(scoreTable,params,sCond,"score","score.sd",srOverride,bnOverride)
      
      logInfo("Writing AA-centric table to file.")
      outFile <- paste0(outDir,sCond,"_t",tp,"_simple_aa_floored.csv")
      cat("# FLOORED AMINO ACID-LEVEL SCORES #\n",paramHeader,file=outFile,sep="")
      suppressWarnings(
        write.table(flooredAA,outFile,sep=",",append=TRUE,row.names=FALSE,qmethod="double")
      )
      
      outFile <- paste0(outDir,sCond,"_t",tp,"_simple_aa.csv")
      cat("# UNFLOORED AMINO ACID-LEVEL SCORES #\n",paramHeader,file=outFile,sep="")
      suppressWarnings(
        write.table(unflooredAA,outFile,sep=",",append=TRUE,row.names=FALSE,qmethod="double")
      )
      
      
    }
  }
  
}

#' Infer a function that models expected error relative to data point magnitutes
#'
#' @param values underlying datapoints
#' @param errors the error associated with each datapoint
#'
#' @return a function that maps data point magnitude to expected error
#' @export
expectedError <- function(values,errors) {
  #calculate running median error
  runningMedian <- as.df(lapply(sort(values), function(mid) {
    grp <- which(abs(values-mid) < 0.1)
    c(value=mid,medErr=median(errors[grp],na.rm=TRUE))
  }))
  #fit a smooth polynomial line to the median trail
  # x <- runningMedian$value
  # z <- lm(y~.,data=fin(data.frame(y=log10(runningMedian$medErr),x1=x,x2=x^2,x3=x^3)))
  
  #fit a mirrored logarithm function to the median trail
  theta0 <- c(.5,.1,1,0)
  modelFun <- function(x,theta=theta0) {
    mirror <- theta[[1]]
    i <- theta[[2]]
    m1 <- theta[[3]]
    a1 <- theta[[4]]
    m1*log10(abs(x-mirror)+i)+a1
  }
  # curve(modelFun)
  objFun <- function(theta) {
    pred <- modelFun(runningMedian$value,theta)
    out <- mean(abs(pred-log10(runningMedian$medErr)),na.rm=TRUE)
    if (is.na(out)) 1e12 else out
  }
  opt <- optimization::optim_nm(fun = objFun,start = theta0)
  
  #build a function to derive the expected error based on the fit above
  expect.error <- function(xs) {
    # 10^predict.lm(z,newdata=data.frame(x1=xs,x2=xs^2,x3=xs^3))
    10^modelFun(xs,theta=opt$par)
  }
  # plot(values,errors,pch=".",log="y")
  # lines(runningMedian,col=2)
  # curve(expect.error,add=TRUE)
  return(expect.error)
}

#' Calculate the residual error
#' 
#' Models expected error relative to data point magnitudes and then uses the model
#' to calculate the residual error, i.e. how much greater or lesser the error is
#' compared to expectation.
#'
#' @param values vector of underlying datapoints
#' @param errors the error associated with each datapoint
#'
#' @return the residual error vector
#' @export
residualError <- function(values, errors) {
  
  #build a function to derive the expected error based on the fit above
  expect.error <- expectedError(values, errors)
  
  #calculate residual error
  residual.error <- errors-expect.error(values)
  #remove excessive negative outliers
  q01 <- quantile(residual.error,0.01,na.rm=TRUE)
  residual.error <- sapply(residual.error,max,q01)
  
  return(residual.error)
}


#' Helper function to extract scaling pivot point override values from the parameter sheet
#' if they exist
#' @param params the parameter sheet
#' @param sCond the current selective condition ID
#' @param tp the time point ID
#' @param region the region ID
#' @return a vector with syn and stop overrides, or \code{NA} if they weren't defined.
#' @export
getPivotOverrides <- function(params,sCond,tp,region) {
  pullValue <- function(type) {
    with(params$pivots,{
      row <- which(Condition == sCond & `Time point` == tp & Region == region & Type == type)
      if (length(row)==0) NA else Value[row]
    })
  }
  if (length(params$pivots) == 0 || nrow(params$pivots) == 0) {
    return(c(syn=NA,non=NA))
  } else {
    return(c(syn=pullValue("synonymous"),non=pullValue("nonsense")))
  }
}

checkPivots <- function(params) {
  helpMsg <- "==>Please either define pivots or use the autoPivot flag.<=="
  if (is.null(params$pivots)) {
    stop("No scale pivots defined in parameter sheet!\n",helpMsg)
  }
  errors <- do.call(c,lapply(getSelects(params), function(sCond) {
    do.call(c,lapply(getTimepointsFor(sCond,params),function(tp) {
      missing <- sapply(params$regions$`Region Number`, function(regi) {
        any(is.na(getPivotOverrides(params,sCond,tp,regi)))
      })
      if (any(missing)) {
        sprintf("No scale pivots defined for region %i in %s time-point %s",
                regi,sCond,tp
        )
      } else NULL
    }))
  }))
  if (length(errors) > 0){
    stop(paste(c(errors,helpMsg),collapse="\n"))
  }
}

#' Calculate scores by scaling bce values to syn/nonsense medians
#' @param msc data.frame containing the enrichment ratios (phi) and log(phi)
#' @param aac the vector of corresponding amino acid changes
#' @param sdThreshold stdev threshold for finding the syn/stop means
#' @return data.frame with variant type, score and stdev
#' @export
enactScale <- function(msc,aac,sdThreshold,overrides=c(syn=NA,non=NA)) {
  
  # #determine mutation types
  # fromAA <- substr(aac,1,1)
  # toAA <- substr(aac,nchar(aac),nchar(aac))
  # msc$type <- as.vector(mapply(function(from,to) {
  #   if (from==to) "synonymous" else if (to=="*") "nonsense" else "missense"
  # },fromAA,toAA))
  # 
  
  #apply filter
  mscFiltered <- msc[is.na(msc$filter),]
  if (!all(is.na(mscFiltered$bce.sd))) {
    resErr <- residualError(mscFiltered$bce,mscFiltered$bce.sd)
    mscFiltered$resErr <- resErr+min(resErr,na.rm=TRUE)
    numNsSurvive <- with(mscFiltered,sum(is.na(filter) & resErr < sdThreshold & type == "nonsense",na.rm=TRUE))
    if (numNsSurvive >= 10) {
      with(mscFiltered,mscFiltered[which(is.na(filter) & resErr < sdThreshold),])
    } else {
      nonsenseSDs <- with(mscFiltered,resErr[is.na(filter) & type=="nonsense"])
      # q10Threshold <- quantile(with(msc,bce.sd[is.na(filter) & type=="nonsense"]),0.1)
      r10Threshold <- sort(nonsenseSDs)[[min(10,length(nonsenseSDs))]]
      logWarn(sprintf("sdThreshold %.03f is too restrictive! Using sd < %.03f instead.",sdThreshold,r10Threshold))
      mscFiltered <- with(mscFiltered,mscFiltered[which(is.na(filter) & resErr < r10Threshold),])
    }
  } else {
    logWarn("No SD-filter applied to mode finder, as no error estimates are available.")
  }
  
  #calculate medians
  synonymousMedian <- with(mscFiltered,median(
    bce[which(type == "synonymous" & bce > 0)]
    ,na.rm=TRUE))
  nonsenseMedian <- with(mscFiltered,median(
    bce[which(type == "nonsense" & bce < 1)]
    ,na.rm=TRUE))
  
  logInfo(sprintf("Auto-detected nonsense median: %.03f",nonsenseMedian))
  logInfo(sprintf("Auto-detected synonymous median: %.03f",synonymousMedian))
  
  #FIXME: If manual pivot points are provided, this check is moot!
  if (is.na(nonsenseMedian) || is.na(synonymousMedian)) {
    logWarn(
      "No nonsense or no synonymous variants of sufficient quality were found!\n",
      "!!Scaling is impossible! Most downstream analyses will not work!!"
    )
    return(data.frame(type=msc$type,score=NA,score.sd=NA))
  }
  
  #apply any potential overrides
  if (!is.na(overrides[["syn"]])) {
    logInfo(sprintf(
      "Applying manual synonymous pivot: %.03f -> %.03f",
      synonymousMedian, overrides[["syn"]]
    ))
    synonymousMedian <- overrides[["syn"]]
  }
  if (!is.na(overrides[["non"]])) {
    logInfo(sprintf(
      "Applying manual nonsense pivot: %.03f -> %.03f",
      nonsenseMedian, overrides[["non"]]
    ))
    nonsenseMedian <- overrides[["non"]]
  }
  
  #safety check, that the medians are in the right orientation
  if (synonymousMedian > nonsenseMedian) {
    #use medians to normalize
    score <- (msc$bce - nonsenseMedian) / (synonymousMedian - nonsenseMedian)
    score.sd <- msc$bce.sd / (synonymousMedian - nonsenseMedian)
  } else {
    #otherwise, we CANNOT assign a correct score!
    logWarn(paste(
      "Synonymous median fell below nonsense median! This means:",
      " * Scores cannot be calculated!",
      " * Most downstream analyses will not work!",
      sep="\n"
    ))
    score <- score.sd <- rep(NA,nrow(msc))
  }
  
  return(data.frame(score=score,score.sd=score.sd))
}

#' Apply flooring and SD adjustment to negative scores
#' 
#' This assumes that no cell can be "deader than dead", and that negative scores are just due to 
#' measurement error. From this follows that the more negative the score, the more certain we 
#' can be that the cell is dead. We can therefore adjust the SD to reflect this, by setting 
#' the score to zero and shifting the SD such that the implied probability of being above 0 
#' remains the same.
#' 
#' @param msc the score matrix
#' @param params the global parameter object
#' @return a matrix containing the floored scores and SDs
flooring <- function(msc,params) {
  #depending on the selection type of the assay, we're either flooring at the top or bottom.
  switch(params$assay[["selection"]],
         Positive={
           #the target null-like score towards which we will shift these values
           targetScore <- 0
           #the quantile for which we want to keep the p-value fixed
           quantile <- 1
           #the row numbers containing the cases to be fixed
           toFix <- which(msc$score < targetScore)
         },
         Negative={
           #if we're dealing with an inverse assay, we have to apply a ceiling instead of flooring
           #the target functional (but dead) score towards which we will shift these values
           targetScore <- 1
           #the quantile for which we want to keep the p-value fixed
           quantile <- 0
           #the row numbers containing the cases to be fixed
           toFix <- which(msc$score > targetScore)
         }
  )
  score.floored <- msc$score
  sd.floored <- msc$score.sd
  score.floored[toFix] <- targetScore
  #the equivalent sds of a normal distribution with the target mean based on the above area
  sd.floored[toFix] <- with(msc[toFix,], score.sd*(quantile-targetScore)/(quantile-score))
  
  return(cbind(score.floored=score.floored,sd.floored=sd.floored))
}

#' join multiple datapoints weighted by stdev
#' @param ms data means
#' @param sds data stdevs
#' @param dfs degrees of freedom for each data point
#' @param ws datapoint weights. By default these get auto-calculated from sds
#' @return a vector containing the joint mean and joint stdev
join.datapoints <- function(ms,sds,dfs,ws=(1/sds)/sum(1/sds)) {
  #weighted mean
  mj <- sum(ws*ms)
  #weighted joint variance
  vj <- sum(ws*(sds^2+ms^2)) -mj^2
  #return output
  return(c(mj=mj,sj=sqrt(vj),dfj=sum(dfs)))
}


#' Collapse score table by amino acid changes
#' 
#' @param scoreTable the full score
#' @param params the parameter object
#' @param sCond the current selective condition
#' @param scoreCol the name of the column containing the scores
#' @param sdCol the name of the column containing the stdev
#' @param sdOverride the sdOverride flag
#' @return a \code{data.frame} containing the collapsed score table
#' @export
collapseByAA <- function(scoreTable,params,sCond,scoreCol="score",sdCol="score.sd",srOverride=FALSE,bnOverride=FALSE) {
  
  #configure filter level for export
  exportFilter <- if (bnOverride) {
    sapply(scoreTable$filter,function(f) is.na(f) || f=="bottleneck")
  } else {
    is.na(scoreTable$filter)
  }
  
  filteredTable <- scoreTable[exportFilter,]
  
  if (!srOverride) {
    resErr <- residualError(filteredTable[,scoreCol],filteredTable[,sdCol])
    nerr <- resErr - min(resErr,na.rm=TRUE) + min(abs(resErr),na.rm=TRUE)
    aaTable <- as.df(tapply(1:nrow(filteredTable),filteredTable$hgvsp, function(is) {
      joint <- join.datapoints(
        ms=filteredTable[is,scoreCol],
        sds=filteredTable[is,sdCol],
        dfs=rep(params$numReplicates[[sCond]],length(is)),
        ws=(1/nerr[is])/sum(1/nerr[is])
      )
      list(
        hgvs_pro=unique(filteredTable[is,"hgvsp"]),
        score=joint[["mj"]],
        sd=joint[["sj"]],
        df=joint[["dfj"]],
        se=joint[["sj"]]/sqrt(joint[["dfj"]])
      )
    }))
  } else {
    aaTable <- as.df(tapply(1:nrow(filteredTable),filteredTable$hgvsp, function(is) {
      list(
        hgvs_pro=unique(filteredTable[is,"hgvsp"]),
        score=mean(fin(filteredTable[is,scoreCol])),
        sd=NA,
        df=length(fin(filteredTable[is,scoreCol]))
      )
    }))
  }
  return(aaTable)
}


