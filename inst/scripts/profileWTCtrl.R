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
#This is a helper script to identify biases in WT ctrl mutations
##############################################

options(stringsAsFactors=FALSE)

#load libraries
library(tileseqMave)
library(argparser)
library(yogilog)
library(yogitools)

#process command line arguments
p <- arg_parser(
  "Analyze WT ctrl mutations",
  name="profileWTCtrl.R"
)
p <- add_argument(p, "--workspace", help="workspace data directory. Defaults to current working directory")
p <- add_argument(p, "--input", help="input directory containing count data. Defaults to subdirectory _mut_count")
p <- add_argument(p, "--parameters", help="Parameter sheet. Defaults to parameters.json in workspace directory")
p <- add_argument(p, "--outfile", help="output file prefix. Defaults to wtProfile in the workspace directory")
p <- add_argument(p, "--logfile", help="log file. Defaults to profileWTCtrl.log in the same directory.")
p <- add_argument(p, "--srOverride", help="Manual override to allow singleton replicates. USE WITH EXTREME CAUTION!",flag=TRUE)
args <- parse_args(p)


dataDir <- args$workspace
if (is.na(dataDir)) {
  dataDir <- getwd()
}
if (!grepl("/$",dataDir)) {
  dataDir <- paste0(dataDir,"/")
}
if (!dir.exists(dataDir)) {
  stop("Workspace folder does not exist!")
}

logfile <- if (is.na(args$logfile)) paste0(dataDir,"profileWTCtrl.log") else args$logfile
#set up logger and shunt it into the error handler
logger <- new.logger(logfile)
registerLogger(logger)
registerLogErrorHandler(logger)
logVersion()


paramfile <- if (is.na(args$parameters)) paste0(dataDir,"parameters.json") else args$parameters
params <- parseParameters(paramfile, srOverride = args$srOverride)

input <- args$input
if (is.na(input)) {
  latest <- latestSubDir(parentDir=dataDir,pattern="_mut_count$")
  input <- latest[["dir"]]
}
if (!dir.exists(input)) {
  stop("Input ",input," does not exist!")
}
infile <- list.files(input, full.names = TRUE, pattern = "marginalCounts.csv$")
if (!file.exists(infile)) {
  stop("No marginalCounts.csv file found in input directory ",input,
       "Have you run joinCounts.R yet?"
  )
}

outfile <- args$outfile
if (is.na(outfile)) {
  outfile <- paste(dataDir,"wtProfile")
}

def <- params$conditions$definitions
wtCtrls <- unique(def[which(def[,"Relationship"]=="is_wt_control_for"),"Condition 1"])

margs <- read.csv(infile,comment.char="#")
#filter out indels etc
margs <- margs[grepl("^[ACGT]{3}\\d+[ACGT]{3}$",margs$codonChange),]
ncChanges <- as.df(lapply(margs$codonChange, function(cc) {
  n <- nchar(cc)
  list(
    f1=substr(cc,1,1),
    f1=substr(cc,2,2),
    f1=substr(cc,3,3),
    pos=as.integer(substr(cc,4,n-3)),
    t3=substr(cc,n-2,n-2),
    t3=substr(cc,n-1,n-1),
    t3=substr(cc,n,n)
  )
}))
numChanges <- apply(ncChanges[,-4],1, function(ncc) {
  sum(ncc[1:3]!=ncc[4:6])
})

formatVector <- function(x) {
  paste0(
    paste(names(x),collapse="\t"),"\n",
    paste(sprintf("%.02e",x),collapse="\t")
  )
}

tp <- params$timepoints[1,1]

invisible(lapply(wtCtrls,function(wtCtrl) {
  
  reps <- 1:params$numReplicates[[wtCtrl]]
  cnames <- sprintf("%s.t%s.rep%d.frequency",wtCtrl,tp,reps)
  mfreq <- rowMeans(margs[,cnames])
  
  snvVmnvFreq <- tapply(mfreq,numChanges,mean)
  
  ncTable <- do.call(rbind,lapply(1:nrow(ncChanges), function(i) {
    ncPos <- (ncChanges[i,"pos"]-1)*3 + 1:3
    from <- ncChanges[i,1:3]
    to <- ncChanges[i,5:7]
    idx <- which(from!=to)
    hgvs <- sprintf("%d%s>%s",ncPos,from,to)[idx]
    data.frame(hgvs=hgvs,freq=mfreq[[i]])
  }))
  totFreqs <- with(ncTable,tapply(freq,hgvs,sum))
  totIDs <-names(totFreqs)
  totFreqs <- data.frame(
    id=totIDs,
    from=substr(totIDs,nchar(totIDs)-2,nchar(totIDs)-2),
    to=substr(totIDs,nchar(totIDs),nchar(totIDs)),
    pos=as.integer(substr(totIDs,1,nchar(totIDs)-3)),
    totFreq=totFreqs
  )
  perNCFreq <- with(totFreqs,tapply(totFreq,paste0(from,">",to),mean))
  perBaseFreq <- tapply(perNCFreq,substr(names(perNCFreq),1,1),sum)  
  
  con <- file(paste0(outfile,"_",wtCtrl,".tsv"),open="w")
  cat(
    "#Rates of SVNs and MNVs per codon",
    formatVector(setNames(snvVmnvFreq,paste0(1:3,"-nt"))),
    "\n#Mutation rate per original base",
    formatVector(perBaseFreq),
    "\n#Rate per nucleotide change",
    formatVector(perNCFreq),
    "\n#Individual positional mutation rates",
    file=con,
    sep="\n"
  )
  write.table(totFreqs,con,sep="\t",row.names=FALSE)
  close(con)
  
  pdf(paste0(outfile,"_",wtCtrl,".pdf"),20,4)
  cmap <- yogitools::colmap(c(0,max(totFreqs$totFreq)),c("white","orange"))
  ncs <- toChars("ACGT")
  x <- totFreqs$pos
  y <- sapply(totFreqs$to,function(to)5-which(ncs==to))
  sqCol <- cmap(totFreqs$totFreq)
  opar <- par(mar=c(5,4,1,1))
  plot(NA,type="n",
       xlim=c(.5,max(x)+.5),ylim=c(.5,4.5),
       axes=FALSE,xlab="base position",ylab="base change"
  )
  axis(1)
  axis(2,4:1,ncs)
  rect(.5,.5,max(x)+.5,max(y)+.5,col="gray",border=NA)
  rect(x-.5,y-.5,x+.5,y+.5,col=sqCol,border=NA)
  par(opar)
  dev.off()
  
}))
cat("Done!\n")


