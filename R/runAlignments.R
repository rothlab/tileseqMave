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

runAlignments <- function(dataDir, inDir, outDir=NA, 
  paramFile=paste0(dataDir,"parameters.json"), srOverride=FALSE,
  usePhiX=FALSE) {

  op <- options(stringsAsFactors=FALSE)

  library(yogitools)
  library(yogiseq)

  #Check if required external software is installed
  retVal <- system2("bowtie2","--version")
  if (retVal != 0) {
    stop("bowtie2 does not appear to be installed!")
  }

  if (!grepl("/$",dataDir)) {
    dataDir <- paste0(dataDir,"/")
  }
  if (!dir.exists(dataDir)) {
    #we don't use the logger here, assuming that whichever script wraps our function
    #catches the exception and writes to the logger (or uses the log error handler)
    stop("Data folder does not exist!")
  }
  if (!grepl("/$",inDir)) {
    inDir <- paste0(inDir,"/")
  }
  if (!dir.exists(inDir)) {
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
  
  #if no output directory was defined
  if (is.na(outDir)) {
    outDir <- paste0(dataDir,"alignments/")
  } 
  #make sure it ends in "/"
  if (!grepl("/$",outDir)) {
    outDir <- paste0(outDir,"/")
  }
  #make sure outdir exists
  dir.create(outDir,recursive=TRUE,showWarnings=FALSE)

  logInfo("Using input directory",inDir,"and output directory",outDir)

  # VALIDATE PRESENCE OF FASTQ FILES -------------------

  allFQ <- list.files(inDir,full.names=TRUE,pattern="\\.fastq\\.gz$")
  findFQ <- function(sid) {
    r1Hits <- grep(sprintf("^%s_.*R1.*\\.fastq\\.gz",sid),basename(allFQ))
    r2Hits <- grep(sprintf("^%s_.*R2.*\\.fastq\\.gz",sid),basename(allFQ))
    if (length(r1Hits) != 1 || length(r2Hits) != 1) {
      stop("No unique FASTQ file for sample",sid) 
    }
    c(r1=allFQ[[r1Hits]],r2=allFQ[[r2Hits]])
  }
  sampleFASTQs <- as.df(lapply(params$samples[,"Sample ID"], findFQ))
  if (usePhiX) {
    phiXFASTQ <- findFQ("Undetermined")
  }

  #build bowtie2 library
  dir.create(paste0(outDir,"db/"))
  dbPrefix <- paste0(outDir,"db/reference")
  buildRefLib <- function(seqs,names,dbPrefix) {
    dbFASTA <- paste0(dbPrefix,".fasta")
    con <- file(dbFASTA,open="w")
    yogiseq::writeFASTA(con,with(params$template,setNames(seq,geneName)))
    close(con)
    retVal <- system2("bowtie2-build",c(dbFASTA,dbPrefix))
    if (retVal != 0) {
      stop("Failed to generate bowtie2 library!")
    }
  }

  #submit alignment jobs
  jobids <- sapply(1:nrow(sampleFASTQs),function(i) {
    cmd <- sprintf(
      "bowtie2 --rdg 12,1 --rfg 12,1 --local --fr -x %s -1 %s -2 %s|samtools view -b -o %s%s.bam -",
      sampleFASTQs[i,1], sampleFASTQs[i,2], outDir, params$samples[i,"Sample ID"]
    )
    jobid <- submitJob(cmd)
    return(jobid)
  })
  if (usePhiX) {
    cmd <- sprintf(
      "bowtie2 --rdg 12,1 --rfg 12,1 --local --fr -x %s -1 %s -2 %s|samtools view -b -o %s%s.bam -",
      sampleFASTQs[i,1], sampleFASTQs[i,2], outDir, params$samples[i,"Sample ID"]
    )
  }

  waitForJobs(jobids)
  
  options(op)

}