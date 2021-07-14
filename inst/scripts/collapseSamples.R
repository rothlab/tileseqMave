#!/usr/bin/env Rscript

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

#####################################################
# Collapse samples to combine sequencing runs or conditions
#####################################################

options(
  stringsAsFactors=FALSE,
  ignore.interactive=TRUE,
  future.cmdargs=c()
)

#load libraries
library(tileseqMave)
library(argparser)
library(yogilog)
library(yogitools)

#process command line arguments
p <- arg_parser(
  "Collapses count for arbitrary samples in a tileseq run.",
  name="collapseSamples.R"
)
p <- add_argument(p, "instructions", help="instructions csv file. Required columns: input.file1; input.file2; output.file")
p <- add_argument(p, "--output", help="output directory. ",default="collapsed/")
p <- add_argument(p, "--logfile", help="log file.",default="collapseSamples.log")
p <- add_argument(p, "--noCoverage", help="do not use coverage data", flag=TRUE)
# p <- add_argument(p, "--cores", default=6L, help="number of CPU cores to use in parallel for multi-threading")
p <- add_argument(p, "--silent", help="Turn off message printing to stdout",flag=TRUE)
args <- parse_args(p)

outdir <- args$output

if (!grepl("/$",outdir)) {
  outdir <- paste0(outdir,"/")
}

#set up logger and shunt it into the error handler
logger <- new.logger(args$logfile,stdout=!args$silent)
registerLogger(logger)
registerLogErrorHandler(logger)
logVersion()

if (!canRead(args$instructions)) {
  stop("Unable to read instruction file ",args$instructions)
}
instructions <- read.csv(args$instructions)

dir.create(outdir,recursive=TRUE,showWarnings=FALSE)

logInfo("Validating input")
#validate input instructions
nInputs <- ncol(instructions)-1
if (nInputs < 2) {
  stop("Must have at least 2 input columns!")
}
if (colnames(instructions)[[nInputs+1]] != "output.file") {
  stop("Last column in instructions file must be 'output.file'!")
}
inCols <- paste0("input.file",1:nInputs)
if (any(colnames(instructions)[1:nInputs] != inCols)) {
  stop("Input columns must be called 'input.file1', 'input.file2', etc.!")
}

invisible(lapply(inCols, function(inCol) {
  filenames <- na.omit(instructions[,inCol])
  if (!all(yogitools::canRead(filenames))) {
    culprits <- filenames[which(!yogitools::canRead(filenames))]
    stop("File(s) ",paste(culprits,collapse=", ")," cannot be read!")
  }
  valid <- grepl("^counts_sample_.+\\.csv$",basename(filenames))
  if (!all(valid)) {
    stop("All input files must follow naming convention 'counts_sample_<ID>.csv'")
  }
}))

if (!args$noCoverage) {
  logInfo("Finding corresponding coverage data")
  #find the corresponding coverage files
  covInstructions <- data.frame(lapply(inCols, function(inCol) {
    cfiles <- instructions[,inCol]
    sapply(cfiles, function(cfile) {
      if (is.na(cfile)) {
        return(NA)
      }
      cdir <- dirname(cfile)
      cname <- basename(cfile)
      id <- gsub("^counts_sample_|\\.csv$","",cname)
      covFile <- paste0(cdir,"/coverage_",id,".csv")
      if (!yogitools::canRead(covFile)) {
        stop("Cannot read coverage file ",covFile)
      }
      covFile
    })
  }))
  colnames(covInstructions) <- inCols
}

#helper function to write joint counts to output file
writeCounts <- function(jointCounts,countOutFile) {
  header <- sprintf("#Sample:%s
#Tile:%d
#Condition:%s
#Replicate:%d
#Timepoint:%s
#Raw read depth:%d
#Number of read pairs without mutations:%d
#Number of read pairs did not map to gene:%d
#Number of reads outside of the tile:%d
#Final read-depth:%d
#Total read pairs with mutations:%d
#Comment: Total read pairs with mutations = Read pairs with mutations that passed the posterior threshold
#Comment: Final read-depth = raw read depth - reads didn't map to gene - reads mapped outside of the tile",
    attr(jointCounts,"sample"),
    as.integer(attr(jointCounts,"tile")),
    attr(jointCounts,"condition"),
    as.integer(attr(jointCounts,"replicate")),
    attr(jointCounts,"timepoint"),
    as.integer(attr(jointCounts,"raw.depth")),
    as.integer(attr(jointCounts,"wtpairs")),
    as.integer(attr(jointCounts,"unmapped")),
    as.integer(attr(jointCounts,"mismapped")),
    as.integer(attr(jointCounts,"depth")),
    as.integer(attr(jointCounts,"mutpairs"))
  )
  con <- file(countOutFile,open="w")
  cat(header,file=con)
  write.csv(jointCounts,con,row.names=FALSE,quote=FALSE)
  close(con)
}

#process instructions
invisible(lapply(1:nrow(instructions), function(i) {
  #pull up relevant files
  countFiles <- instructions[i,inCols]
  countOutFile <- paste0(outdir,instructions[i,nInputs+1])
  outId <- gsub("^counts_sample_|\\.csv$","",countOutFile)
  
  logInfo("Processing ",outId)
  
  #parse counts
  counts <- lapply(countFiles, parseCountFile)
  #exract header information
  countHeaders <- as.df(lapply(counts,function(cs) attributes(cs)[-(1:3)]))
  
  if (length(unique(countHeaders$tile)) > 1) {
    stop("Proposed merger is across different tiles!")
  }
  
  #construct joint header
  outHeader <- list(
    sample=outId,
    tile=unique(countHeaders$tile),
    condition=paste(unique(countHeaders$condition),collapse=""),
    replicate=min(as.integer(countHeaders$replicate)),
    timepoint=paste(unique(countHeaders$timepoint),collapse=""),
    raw.depth=sum(as.numeric(countHeaders$raw.depth)),
    depth=sum(as.numeric(countHeaders$depth)),
    wtpairs=sum(as.numeric(countHeaders$wtpairs)),
    mutpairs=sum(as.numeric(countHeaders$mutpairs)),
    unmapped=sum(as.numeric(countHeaders$unmapped)),
    mismapped=sum(as.numeric(countHeaders$mismapped))
  )
  
  #extract all unique variants
  allHGVS <- Reduce(union,lapply(counts,`[[`,"HGVS"))
  
  #merge count tables
  allCounts <- do.call(cbind,lapply(counts,function(ctab){
    rownames(ctab) <- ctab$HGVS
    cs <- ctab[allHGVS,"count"]
    cs[is.na(cs)] <- 0
    cs
  }))
  jointCounts <- data.frame(
    HGVS=allHGVS,
    count=rowSums(allCounts, na.rm=TRUE)
  )
  attributes(jointCounts) <- c(attributes(jointCounts),outHeader)
  
  #write to file
  writeCounts(jointCounts,countOutFile)
  
  
  if (!args$noCoverage) {
    coverFiles <- covInstructions[i,inCols]
    
    coverOutFile <- if (grepl("counts_sample_",countOutFile)) {
      sub("counts_sample_","coverage_",countOutFile)
    } else {
      paste0(dirname(countOutFile),"/coverage_",basename(countOutFile))
    }
  
    coverages <- lapply(coverFiles, read.csv)
    
    allPos <- sort(Reduce(union,lapply(coverages,`[[`,"pos")))
    
    covMatrix <- do.call(zbind,lapply(coverages, function(covi) {
      rownames(covi) <- covi$pos
      #FIXME: this might trigger the approx-match bug!!
      as.matrix(covi[as.character(allPos),-1])
    }))
    covJoint <- apply(covMatrix,1:2,sum)
    covOut <- data.frame(pos=allPos,covJoint)
    
    #write to file
    write.csv(covOut,coverOutFile,row.names=FALSE,quote=FALSE)
  
  }
  
}))

logInfo("Done!")
