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

##############################################
#This script calculates calibration probabilities for PHRED scores
##############################################

options(
  stringsAsFactors=FALSE,
  ignore.interactive=TRUE,
  future.cmdargs=c()
)

library(yogitools)
library(tileseqMave)
library(argparser)
library(bitops)
library(yogilog)
library(parallel)

p <- argparser::arg_parser(
  "Calibrate PHRED scores in a WTctrl SAM file to their empirical error probabilities",
  name="calibratePhred.R"
)

p <- add_argument(p, "sam", help="the SAM file to analyze. Must be a WT ctrl alignment.")
p <- add_argument(p, "--output", help="output file. Defaults to name of the <sam>_phredCalibration.csv")
p <- add_argument(p, "--parameters", help="parameter file containing the reference sequence. Defaults to parameters.json in the current directory.")
p <- add_argument(p, "--fastaref", help="fasta file containing the reference sequence. This can be given as an alternative to the parameter sheet.")
p <- add_argument(p, "--logfile", help="log file. Defaults to 'calibratePhred.log' in the same directory")
p <- add_argument(p, "--srOverride", help="Manual override to allow singleton replicates. USE WITH EXTREME CAUTION!",flag=TRUE)
p <- add_argument(p, "--maxReads", help="Maximum number of reads to process. (To finish faster)",default=1e5)
p <- add_argument(p, "--silent", help="Turn off message printing to stdout",flag=TRUE)
p <- add_argument(p, "--cores", default=6L, help="Number of CPU cores to use")
args <- parse_args(p)


paramFile <- if (is.na(args$parameters)) "parameters.json" else args$parameters
fastaFile <- args$fastaref
logfile <- if (is.na(args$logfile)) "calibratePhred.log" else args$logfile

#set up logger and shunt it into the error handler
logger <- new.logger(logfile,stdout=!args$silent)
registerLogger(logger)
registerLogErrorHandler(logger)
logVersion()

sam.file <- args$sam

outfile <- if (is.na(args$output)) {
  sub(".sam$","_phredCalibration.csv",basename(sam.file))
} else {
  args$output
}

if (!canRead(sam.file)) {
  stop("Cannot read SAM file ",sam.file,"!")
}
if (is.na(fastaFile) || !canRead(fastaFile)) {
  logWarn("Unable to read reference FASTA file, defaulting to parameter sheet...")
  if (!canRead(paramFile)) {
    stop("Cannot read parameter file ",paramFile,"!")
  }
  #parse parameter file
  params <- tileseqMave::parseParameters(paramFile,srOverride=args$srOverride)
  wtseq <- params$template$seq
} else {
  fastaLines <- readLines(fastaFile)
  if (sum(grepl("^>",fastaLines)) != 1) {
    stop("FASTA file must contain exactly one single entry!")
  }
  wtseq <- paste(trimws(fastaLines[-1]),collapse="")
  if (!grepl("^[ACGT]+$",wtseq)) {
    stop("Reference sequence is not a valid DNA sequence! (Must only contain A,C,G,T)")
  }
}

#set readCutoff
maxReads <- as.numeric(args$maxReads)
if (is.na(maxReads) || maxReads < 1) {
  maxReads <- Inf
}

#helper function to convert Phred scores to specified error probabilities
phredToProb <- function(phred) {
  if (phred == " ") return(NA)
  phredNum <- as.integer(charToRaw(phred))-33
  10^(-phredNum/10)
}
#helper vectors of valid Phred scores and their probabilities
phredChars <- toChars(" !\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ")
phredProbs <- sapply(phredChars,phredToProb)
#SAM specification flags and their bit masks
flagMasks <- c(
  multiSegment=0x1, allSegmentsOK=0x2, segmentUnmapped=0x4,
  nextSegmentUnmapped=0x8, revComp=0x10, nextRevComp=0x20, 
  firstSegment=0x40, lastSegment=0x80, secondary=0x100, 
  failQC=0x200, duplicate=0x400, supplementary=0x800
)

#helper function to parse the SAM tab-delimited format
#read.delim doesn't work, because lines can have differing
#numbers of columns!
parseSpecialTab <- function(lines) {
  splits <- strsplit(lines,"\\t")
  maxCols <- max(sapply(splits,length))
  df <- as.df(lapply(splits, function(elem) {
    if (length(elem) < maxCols) {
      c(elem,rep(NA,maxCols-length(elem)))
    } else elem
  }))
  if (maxCols > 12) {
    colnames(df) <- c(
      "cname","flag","rname","pos","mapq","cigar","mrnm","mpos",
      "isize","seq","qual","tags", 13:maxCols
    )
  } else {
    colnames(df) <- c(
      "cname","flag","rname","pos","mapq","cigar","mrnm","mpos",
      "isize","seq","qual","tags"
    )
  }
  df$flag <- as.integer(df$flag)
  df$pos <- as.integer(df$pos)
  df$mapq <- as.integer(df$mapq)
  df$mpos <- as.integer(df$mpos)
  df$isize <- as.integer(df$isize)
  df
}

#read sam file
logInfo("Processing SAM file")

# TODO: it may be useful to sample random lines e.g. 
# shuf -n $maxReads <file>
# also, potentially enriching for matches using samtools?

#set up multithreading
nthreads <- args$cores

#scope out the file size
#the sam files are roughly 430-465 bytes per line, but we will purposefully
#overestimate this as 500bytes per line to get a conservative number
totalLines <- round(file.info(sam.file)$size/500)
#if the approximate number is close to our cutoff, we want a more precise number
if (totalLines < 1.2*maxReads) {
  totalLines <- as.integer(strsplit(
    system2("wc",c("-l",sam.file),stdout=TRUE),
    "\\s+"
  )[[1]][[1]])
} 

#if the file is smaller than our maximum line number, it becomes our new maximum
maxReads <- min(maxReads,totalLines)
#subdivide the task into chunks (the number of lines to process by each thread)
chunkSize <- floor(maxReads/nthreads)

#start multiple threads
result <- mclapply(1:nthreads,function(threadi) {
  
  startLine <- chunkSize*(threadi-1)+1
  endLine <- startLine+chunkSize-1
  
  #set up global count matrix to be filled during chunk processing
  counts <- matrix(0,nrow=length(phredChars),ncol=2,dimnames=list(phredChars,c("err","tot")))
  
  #keep track of how many lines were processed or errors were encountered
  events <- yogitools::new.counter()
  
  # con <- file(sam.file,open="r")
  # pipe(sprintf("tail +%d %s|head -%d",startLine,sam.file,chunkSize))
  con <- pipe(sprintf("sed -n %d,%dp %s",startLine,endLine,sam.file),open="r")

  #read file in 1000 line chunks
  while(length(lines <- readLines(con,1000))>0) {
    
    #parse lines
    sam <- parseSpecialTab(lines)
    
    #parse CIGAR strings
    cigar <- yogitools::global.extract.groups(sam$cigar,"(\\d+)([SHNMDIP]{1})")
    
    #parse unmapped flag
    unmapped <- bitops::bitAnd(sam$flag,flagMasks[["segmentUnmapped"]])>0
    
    #if nothing maps, there's nothing worth counting here.
    if (all(unmapped)) {
      next
    }
    
    #inner loop to process individual lines in the chunk
    for (i in 1:nrow(sam)) {
      
      #if the segment is unmapped, ignore it.
      if (is.na(unmapped[[i]]) || unmapped[[i]]) {
        next
      }
      
      #the read sequence
      rseq <- sam[i,"seq"]
      #the phred quality track
      qual <- sam[i,"qual"]
      
      #a positional cursor keeping track of where we are in the template sequence
      tcursor <- sam[i,"pos"]
      #a position cursor for where we are in the read sequence
      rcursor <- 1
      
      #string buffers for the expanded template, read, and quality alignments
      tbuffer <- ""
      rbuffer <- ""
      qbuffer <- ""
      
      #load the parsed cigar string
      cig <- cigar[[i]]
      
      #variable for keeping track of error states in the inner loop
      cigarError <- FALSE
      
      #process all positions in this line
      for (j in 1:nrow(cig)) {
        #k is the length of the current cigar element (i.e. match length etc)
        k <- as.integer(cig[j,1])
        #op is the cigar operation (i.e. match / deletion / insertion etc)
        op <- cig[j,2]
        #build the alignment in the string buffers based on the operations
        switch(op,
               M={#match/mismatch
                 tbuffer <- paste0(tbuffer,substr(wtseq,tcursor,tcursor+k-1))
                 rbuffer <- paste0(rbuffer,substr(rseq,rcursor,rcursor+k-1))
                 qbuffer <- paste0(qbuffer,substr(qual,rcursor,rcursor+k-1))
                 tcursor <- tcursor+k
                 rcursor <- rcursor+k
               },
               D={#deletion
                 spaces <- paste(rep(" ",k),collapse="")
                 tbuffer <- paste0(tbuffer,substr(wtseq,tcursor,tcursor+k-1))
                 rbuffer <- paste0(rbuffer,spaces)
                 qbuffer <- paste0(qbuffer,spaces)
                 tcursor <- tcursor+k
               },
               I={#insertion
                 spaces <- paste(rep(" ",k),collapse="")
                 tbuffer <- paste0(tbuffer,spaces)
                 rbuffer <- paste0(rbuffer,substr(rseq,rcursor,rcursor+k-1))
                 qbuffer <- paste0(qbuffer,substr(qual,rcursor,rcursor+k-1))
                 rcursor <- rcursor+k
               },
               S={#soft clipping
                 rcursor <- rcursor+k
               },
               H={#hard clipping
                 rcursor <- rcursor+k
               },
               {#other
                 # logWarn("Unsupproted cigar character: ",op," in row ",chunk+i)
                 cigarError <- TRUE
               }
        )
      }
      
      if (cigarError) {
        events$inc("cigarError")
        next
      }
      
      matches <- data.frame(t=toChars(tbuffer),r=toChars(rbuffer),q=toChars(qbuffer))
      matches$err <- with(matches,t != r)
      if (sum(matches$err) > 10) {
        events$inc("suspicious")
      }
      icounts <- do.call(rbind,with(matches,tapply(err,q,function(e)c(sum(e),length(e)))))
      counts[rownames(icounts),] <- counts[rownames(icounts),] + icounts
      
    }
    
    events$add("processed",length(lines))
  }
  close(con)
  
  return(list(counts=counts,events=events))
  
},mc.cores=nthreads)

#add all the count matrices together
counts <- Reduce(`+`,lapply(result,`[[`,1))

#join all the event logs together
events <- new.counter()
for (eventLog in sapply(result,function(res)res$events$export())) {
  events$import.add(eventLog)
}
events <- events$ls()

logInfo(sprintf("Successfully processed %d reads.",events$processed))

if ("suspicious" %in% names(events)) {
  logWarn(sprintf("%d reads had over 10 errors each!",events$suspicious))
}

totRate <- sum(counts[,1])/sum(counts[,2])
highRate <- sum(counts[-(1:30),1])/sum(counts[-(1:30),2])
logInfo(sprintf(
  "Total error rate: %.02f%%; high quality error rate: %.02f%%",
  100*totRate,100*highRate
))

logInfo("Writing results to file ",outfile)
result <- data.frame(
  specification=phredProbs,observed=counts[,1]/counts[,2],
  errors=counts[,1],total=counts[,2]
)
rownames(result)[[1]] <- "del"

write.csv(result,outfile)

logInfo("Done!")

