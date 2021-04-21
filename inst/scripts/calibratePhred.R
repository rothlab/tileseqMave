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

p <- argparser::arg_parser(
  "Calibrate PHRED scores in a WTctrl SAM file to their empirical error probabilities",
  name="calibratePhred.R"
)

p <- add_argument(p, "sam", help="the SAM file to analyze. Must be a WT ctrl alignment.")
p <- add_argument(p, "--output", help="output file. Defaults to name of the <sam>_phredCalibration.csv")
p <- add_argument(p, "--parameters", help="parameter file. Defaults to parameters.json in the current directory.")
p <- add_argument(p, "--logfile", help="log file. Defaults to 'calibratePhred.log' in the same directory")
p <- add_argument(p, "--srOverride", help="Manual override to allow singleton replicates. USE WITH EXTREME CAUTION!",flag=TRUE)
p <- add_argument(p, "--maxReads", help="Maximum number of reads to process. (To finish faster)",default=1e5)
p <- add_argument(p, "--silent", help="Turn off message printing to stdout",flag=TRUE)

args <- parse_args(p)


paramFile <- if (is.na(args$parameters)) "parameters.json" else args$parameters
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
if (!canRead(paramFile)) {
  stop("Cannot read parameter file ",paramFile,"!")
}
#parse parameter file
params <- tileseqMave::parseParameters(paramFile,srOverride=args$srOverride)
wtseq <- params$template$seq
# wtseq <- strsplit(scan("MTHFR_seq.txt",what="character",sep="\n",quiet=TRUE)[[1]]," ")[[1]][[2]]

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

#read sam file
logger$info("Processing SAM file")

#set up global count matrix to be filled during chunk processing
counts <- matrix(0,nrow=length(phredChars),ncol=2,dimnames=list(phredChars,c("err","tot")))


#TODO: may be better to pipe 'shuf -n $maxReads <file>"
con <- file(sam.file,open="r")
#read file in 1000 line chunks
chunk <- 0
while(length(lines <- readLines(con,1000))>0) {
  
  #read chunk as into table via textConnection
  tcon <- textConnection(lines)
  sam <- read.delim(tcon,header=FALSE,stringsAsFactors=FALSE)
  close(tcon)
  
  colnames(sam) <- c(
    "cname","flag","rname","pos","mapq","cigar","mrnm","mpos",
    "isize","seq","qual","tags", 13:ncol(sam)
  )
  
  #parse CIGAR strings
  cigar <- yogitools::global.extract.groups(sam$cigar,"(\\d+)([SHNMDIP]{1})")
  
  #parse FLAGS
  flags <- do.call(rbind,lapply(sam$flag,function(x)bitops::bitAnd(x,flagMasks)>0))
  colnames(flags) <- names(flagMasks)
  flags <- to.df(flags)
  
  #inner loop to process individual lines in the chunk
  for (i in 1:nrow(sam)) {
    
    #if the segment is unmapped, ignore it.
    if (is.na(flags[i,"segmentUnmapped"]) || flags[i,"segmentUnmapped"]) {
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
               logWarn("Unsupproted cigar character: ",op," in row ",chunk+i)
               cigarError <- TRUE
             }
      )
    }
    
    if (cigarError) {
      next
    }
    
    matches <- data.frame(t=toChars(tbuffer),r=toChars(rbuffer),q=toChars(qbuffer))
    matches$err <- with(matches,t != r)
    if (sum(matches$err) > 10) {
      logWarn("Suspicious error rate in row ",chunk+i)
    }
    icounts <- do.call(rbind,with(matches,tapply(err,q,function(e)c(sum(e),length(e)))))
    counts[rownames(icounts),] <- counts[rownames(icounts),] + icounts
    
    # setTxtProgressBar(pb, i)
    
  }
  
  # close(pb)
  chunk <- chunk+length(lines)
  logInfo("Processed",chunk,"lines...")
  
  if (chunk > maxReads) {
    logInfo("Exceeded maxReads parameter. Terminating...")
    break
  }
}
close(con)

logger$info("Writing results to file ",outfile)
result <- data.frame(
  specification=phredProbs,observed=counts[,1]/counts[,2],
  errors=counts[,1],total=counts[,2]
)
rownames(result)[[1]] <- "del"

write.csv(result,outfile)

totRate <- sum(counts[,1])/sum(counts[,2])
highRate <- sum(counts[-(1:10),1])/sum(counts[-(1:10),2])
logInfo(sprintf(
  "Total error rate: %.02f%%; high quality error rate: %.02f%%",
  100*totRate,100*highRate
))

# fit <- coefficients(lm(log10(result[!is.na(result[,2]),2])~which(!is.na(result[,2]))))

# idx <- which(!is.na(result[-1,2]))
# 
# pdffile <- sub(".sam$","_qualFreqs.pdf",basename(sam.file))
# pdf(pdffile,11,6)
# op <- par(mfrow=c(2,1),cex=0.9)
# barcols <- c("gray","firebrick3")
# xs <- barplot(t(as.matrix(result[-1,-3])),beside=TRUE,col=barcols,ylim=c(0,1.2),
#               main="Sequencing error rates",xlab="PHRED",ylab="Error probabability"
# )
# text(colMeans(xs)[idx],result[idx+1,2],
#      paste0("n=",trimws(format(result[idx+1,3],big.mark=","))),
#      pos=3,cex=0.7
# )
# grid(NA,NULL)
# legend("right",c("Specificiation","Empirical"),fill=barcols)
# barplot(t(as.matrix(result[-1,-3])),beside=TRUE,col=barcols,log="y",
#         xlab="PHRED",ylab="Error probabability"
# )
# grid(NA,NULL)
# par(op)
# dev.off()

