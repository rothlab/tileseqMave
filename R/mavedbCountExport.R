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

#' Convert count files from legacy pipeline to MaveDB format
#'
#' @param countfile the path to the "rawData.txt" file produced by the legacy pipeline.
#' @param outdir path to desired output directory
#' @return nothing. output is written to various files in the output directory
#' @export
mavedbCountExport <-  function(countfile,outdir,logger=NULL) {

  library(hgvsParseR)
  # library(yogilog)
  library(yogitools)

  options(stringsAsFactors=FALSE)

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

  ##############
  # Read and validate input data
  ##############

  #countfile <- "/home/jweile/projects/ccbr2hgvs/HMGCR_S_resultfile/rawData.txt"
  canRead <- function(filename) file.access(filename,mode=4) == 0
  stopifnot(
    canRead(countfile)
  )

  rawCounts <- read.delim(countfile)

  stopifnot(
    c(
      "wt_codon","pos","mut_codon","wt_aa","mut_aa",
      "nonselect1","nonselect2","select1","select2",
      "controlNS1","controlNS2","controlS1","controlS2"
    ) %in% colnames(rawCounts)
  )

  #make sure outdir ends with a "/"
  if (!grepl("/$",outdir)) {
    outdir <- paste0(outdir,"/")
  }
  #and if it doesn't exist, create it
  if (!dir.exists(outdir)) {
    dir.create(outdir,recursive=TRUE)
  }

  ##################
  # Build HGVS variant descriptor strings
  ##################

  #for nucleotide level descriptors
  cbuilder <- new.hgvs.builder.c()
  #for protein-level descriptors
  pbuilder <- new.hgvs.builder.p(aacode=3)

  #for nucleotide level descriptors
  hgvsc <- apply(rawCounts,1,function(row) {
    pos <- as.numeric(row["pos"])
    #codon start indices
    cstart <- pos*3-2
    wt <- row["wt_codon"]
    mut <- row["mut_codon"]
    #calculate differences between codons
    diffs <- sapply(1:3,function(i)substr(wt,i,i)!=substr(mut,i,i))
    ndiff <- sum(diffs)
    if (ndiff == 1) { #one difference is a SNV
      offset <- which(diffs)
      wtbase <- substr(wt,offset,offset)
      mutbase <- substr(mut,offset,offset)
      snvpos <- cstart+offset-1
      return(cbuilder$substitution(snvpos,wtbase,mutbase))
    } else if (ndiff > 1) { #multiple differences is a delIns
      return(cbuilder$delins(cstart,cstart+2,mut))
    } else {
      stop("mutation must differ from wt!")
    }
  })

  hgvsp <- apply(rawCounts,1,function(row) {
    pos <- as.numeric(row["pos"])
    wt <- row["wt_aa"]
    mut <- row["mut_aa"]
    if (mut=="_") mut <- "*" #correct stop character
    if (wt == mut) {
      return(pbuilder$synonymous(pos,wt))
    } else {
      return(pbuilder$substitution(pos,wt,mut))
    }
  })

  logInfo(sprintf(
    "Parsed data for %d variants covering %d amino acid changes",
    length(hgvsc),length(unique(hgvsp))
  ))

  logInfo("Writing output to file.")

  outTable <- cbind(hgvs_nt=hgvsc,hgvs_pro=hgvsp,rawCounts[,-(1:6)])

  outfile <- paste0(outdir,"mavedb_counts_perNt.csv")
  write.csv(outTable,outfile,row.names=FALSE)

}

