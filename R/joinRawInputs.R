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

#' join raw input datasets
#' 
#' Joins multiple input datasets into a single one, treating each as additional replicates
#' 
#' @param countfiles A vector of filenames for the input count files
#' @param outfile The path of the output file. Defaults to NULL, in which case the
#'      output is written to a file in the same directory as the firsts input file, named
#'      similarly to said input file with an added "_merged" tag.
#' @param logger A yogilogger object, may be omitted, in which case messages are written
#'      directly to stdout.
#' @return the merged data as a \code{data.frame}.
#' @export
#' @examples
# ' countfiles <- sprintf("workspace/input/CBS/rawData_B6-%d_Q5_%s.txt",
# '   c(0,0,1,1),
# '   rep(c("new","old"),2)
# ' )
# ' joinRawInputs(countfiles)
#'
joinRawInputs <- function(countfiles,outfile=NULL,logger=NULL) {

  options(stringsAsFactors=FALSE)

  library(yogitools)
  library(hash)
  library(pbmcapply)

  ################################
  # Check validity of parameters #
  # And set up helper functions  #
  ################################

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
      logger$fatal(...)
    } else {
      do.call(cat,c("Error:",list(...),"\n"))
    }
  }

  canRead <- function(filename) file.access(filename,mode=4) == 0

  #check that all input files exist
  if (!all(canRead(countfiles))) {
    missing <- countfiles[which(!file.exists(countfiles))]
    errMsg <- paste("Input file(s)",missing,"do(es) not exist!")
    logErr(errMsg)
    stop(errMsg)
  }

  #check that there at least two inputs
  if (length(countfiles) < 2) {
    errMsg <- "Less than two input files. Nothing to merge!"
    logErr(errMsg)
    stop(errMsg)
  }

  #if no outputfile is given, derive from first input file
  if (is.null(outfile)) {
    outfile <- sub("\\.(\\w+)$","_merged.txt",countfiles[[1]])
  } 

  #if output directory doesn't exit, create it (with warning).
  outdir <- sub("[^/]+$","",outfile)
  if (!dir.exists(outdir)) {
    logWarn(outdir,"does not exist. Creating now.")
    dir.create(outdir,recursive=TRUE)
  }


  ##############################################
  # Read input data and perform pre-processing #
  ##############################################

  logInfo("Reading input files.")

  #read input files into list of dataframes
  countdata <- lapply(countfiles,read.delim)

  #check that all columns are where they are supposed to be
  hasGoodFormat <- sapply(countdata, function(rawCounts) {
    identical(colnames(rawCounts)[1:6],c(
      "wt_aa","pos","mut_aa","wt_codon","mut_codon","annotation"
    )) &&
    all(grepl("^\\w+\\d+$",colnames(rawCounts)[-(1:6)]))
  })
  if (!all(hasGoodFormat)) {
    problemFiles <- paste(countfiles[[which(!hasGoodFormat)]],collapse=", ")
    errMsg <- paste(problemFiles, "contain(s) invalid column names!")
    logErr(errMsg)
    stop(errMsg)
  }

  # Organize conditions and replicates
  #-----------------------------------

  logInfo("Organizing conditions and replicates")

  #extract condition names and replicates into a list of matrices,
  #tabulating column names by condition and replicate name
  conditionMatrices <- lapply(countdata, function(rawCounts) {
    conditions <- colnames(rawCounts)[-(1:6)]
    condStruc <- extract.groups(conditions,"^(\\w+)(\\d+)$")
    condNames <- unique(condStruc[,1])
    repNames <- unique(condStruc[,2])
    condMatrix <- do.call(rbind,lapply(condNames,paste0,repNames))
    dimnames(condMatrix) <- list(condNames,repNames)
    return(condMatrix)
  })

  #check that the conditions are the same in all inputs
  hasSameConditions <- all(sapply(conditionMatrices[-1],function(condMat) {
    identical(sort(rownames(conditionMatrices[[1]])),sort(rownames(condMat)))
  }))
  if (!hasSameConditions) {
    errMsg <- "The input files do not have the same conditions!"
    logErr(errMsg)
    stop(errMsg)
  }


  #make sure the conditions are in the same order for processing
  condNames <- rownames(conditionMatrices[[1]])
  conditionMatrices <- lapply(conditionMatrices, function(condMat) condMat[condNames,])
  
  #plan output columns
  totalReplicates <- sum(sapply(conditionMatrices,ncol))
  columnPlan <- zbind(
    do.call(cbind,conditionMatrices),
    do.call(cbind,lapply(1:length(conditionMatrices),function(i) {
      apply(conditionMatrices[[1]],1:2,function(x) i)
    }))
  )
  dimnames(columnPlan) <- list(condNames, 1:totalReplicates, c("column","table"))

  # Organize variant sets
  #----------------------

  logInfo("Indexing variants")

  #construct unique variant IDs for each table
  vIds <- lapply(countdata,function(rawCounts) {
    with(rawCounts,paste0(wt_codon,pos,mut_codon))
  })
  #list of all unique variant IDs
  uniqueVids <- Reduce(union,vIds)

  #build a look-up table for the row-number corresponding to each unique ID
  rowIdxs <- do.call(cbind,lapply(vIds, function(vidList) {
    vidIdx <- hash::hash(vidList, 1:length(vidList))
    sapply(uniqueVids, function(uvid) {
      if (hash::has.key(uvid,vidIdx)) {
        vidIdx[[uvid]]
      } else {
        NA
      }
    }) 
  }))


  ##############################
  # Implement the merging plan #
  ##############################

  logInfo("Merging datasets")

  jointTable <- do.call(rbind,pbmclapply(1:length(uniqueVids), function(i) {
    #table j i the first table in which this variant appears
    j <- which(!is.na(rowIdxs[i,]))[[1]]
    metadata <- countdata[[j]][rowIdxs[i,j],1:6]
    jointData <- do.call(c,lapply(condNames, function(condName) {
      out <- mapply(function(tableIdx,colName) {
        countdata[[tableIdx]][rowIdxs[i,tableIdx],colName]
      },
      tableIdx=as.integer(columnPlan[condName,,"table"]),
      colName=columnPlan[condName,,"column"]
      )
      setNames(out,paste0(condName,1:totalReplicates))
    }))
    cbind(metadata,as.list(jointData))
  },mc.cores=8))

  ##############################
  # Write result to filesystem #
  ##############################

  logInfo("Writing results to file")

  write.table(jointTable,outfile,quote=FALSE,sep="\t",row.names=FALSE)

  logInfo("Done!")

  return(jointTable)

}
