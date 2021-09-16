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
#This is a command line tool for locally running mavevis
##############################################

options(stringsAsFactors=FALSE)

library(mavevis)
library(tileseqMave)
library(hgvsParseR)
library(argparser)
library(yogilog)
library(yogitools)


#process command line arguments
p <- arg_parser(
  "Runs mavevis locally on a MaveDB-formatted file to produce a genophenogram.",
  name="mavevisLocal.R"
)
p <- add_argument(p, "--workspace", help="workspace data directory. Defaults to current working directory")
p <- add_argument(p, "--input", help="either input file or input directory containing the count data. Defaults to subdirectory with latest timestamp ending in _scores")
p <- add_argument(p, "--output", help="either output pdf file or output directory. Default adjusts automatically depending on input type")
p <- add_argument(p, "--parameters", help="parameter file. Defaults to parameters.json in the data directory.")
p <- add_argument(p, "--pdb", help="PDB structures. Semicolon-separated list of #-separated pairings between PDB IDs and chain IDs.")
p <- add_argument(p, "--squish", help="output pdf file. Defaults to the name of the input file with pdf extension.",flag=TRUE)
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

paramfile <- if (is.na(args$parameters)) paste0(dataDir,"parameters.json") else args$parameters
params <- parseParameters(paramfile, srOverride = args$srOverride)
uniprot <- params$template$uniprot

input <- args$input
if (is.na(input)) {
  latest <- latestSubDir(parentDir=dataDir,pattern="_scores$")
  input <- latest[["dir"]]
}

if (dir.exists(input)) {
 infiles <- list.files(input, full.names = TRUE, pattern = "_simple_aa.csv$")
 outdir <- args$output
 if (is.na(outdir)) {
   outdir <- input
 }
 dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
} else if (file.exists(input)) {
  if (!grepl("\\.csv$",input)) {
    stop("The given input file is not a CSV file!")
  }
  infiles <- input
  outdir <- NA
} else {
  stop("Input ",input," does not exist!")
}


#get the PDB argument, parse and validate it.
# pdbArg <- getArg("pdb",default=NULL)
if (!is.na(args$pdb)) {
  if (!grepl("#",args$pdb)){
    stop("PDB argument must indicate PDB ID and chain ID separated by a hash character.")
  }
  pdb <- strsplit(args$pdb,";")[[1]]
  pdbIds <- sapply(strsplit(pdb,"#"),`[[`,1)
  pdbChains <- sapply(strsplit(pdb,"#"),`[[`,2)

  if (!all(grepl("^[0-9][A-Za-z0-9]{3}$",pdbIds))) {
    stop("One or more of the given PDB IDs is invalid!")
  }
  if (!all(grepl("^[A-Z]{1}$",pdbChains))) {
    stop("One or more of the PDB chain identifiers is invalid!")
  }
}

#iterate over input files
for (infile in infiles) {
  
  #if this is a single file (not from a directory)
  if (is.na(outdir)) {
    if (is.na(args$output)) {
      pdffile <-  sub("\\.csv$",".pdf",infile)
    } else {
      pdffile <- args$output
      if (!grepl("pdf$",pdffile)) {
        pdffile <- paste0(pdffile,".pdf")
      }
    }
  } else { #if this is one of many files in directory
    pdffile <- paste0(outdir,"/",sub("csv$","pdf",sub(".+/","",infile)))
  }

  cat("Reading and parsing input data...")
  indata <- read.csv(infile,comment.char="#")
  if (!all(c("hgvs_pro","score") %in% colnames(indata))) {
    stop("Input file must be in MaveDB format!")
  }
  mutdata <- parseHGVS(indata$hgvs_pro,aacode=1)
  mutdata[which(mutdata$variant=="*"),"type"] <- "nonsense"
  data <- cbind(mutdata,indata[,-1])
  
  cat("done!\n")

  if (nrow(data)==0 || all(is.na(data$score))) {
    warning("No (valid) data! Skipping...")
    next
  }
  
  #derive WT sequence
  cat("Deriving WT sequence...\n")
  ancestrals <- with(data,tapply(ancestral,start,unique))
  wt.aa <- ancestrals[as.character(1:max(data$start))]
  wt.aa[[1]] <- "M"
  
  if (any(is.na(wt.aa))) {
    warning("Unable to fully derive WT sequence! Defaulting to Uniprot sequence.")
    wt.aa <- yogitools::toChars(mavevis::getUniprotSeq(uniprot))
  }
  
  td <- new.trackdrawer(length(wt.aa),nox=TRUE)
  
  if (!is.null(uniprot)) {
    cons <- calc.conservation(uniprot)
    td$add.constrack(cons)
  }
  
  if (!is.na(args$pdb)) {
  
    strucfeats <- mapply(calc.strucfeats,pdbIds,pdbChains,SIMPLIFY=FALSE)
  
    #consolidate secondary structure information from all structures
    sscols <- lapply(strucfeats,`[`,,"secstruc")
    fillup.max <- max(sapply(sscols,length))
    sscols <- lapply(sscols,function(xs) c(xs,rep(NA,fillup.max-length(xs))))
    if (length(sscols) > 1) {
      ss.consensus <- apply(do.call(cbind,sscols),1,function(xs) if (!all(is.na(xs))) names(which.max(table(xs))) else NA)
    } else {
      ss.consensus <- sscols[[1]]
    }
  
    #consolidate SASA information from all structures
    accols <- lapply(strucfeats,`[`,,"all.rel")
    fillup.max <- max(sapply(sscols,length))
    accols <- lapply(accols,function(xs) c(xs,rep(NA,fillup.max-length(xs))))
    acc.consensus <- apply(do.call(cbind,accols),1,median,na.rm=TRUE)
  
    td$add.ss.track(ss.consensus)
    td$add.track(acc.consensus,"Rel. ASA","steelblue3")
    for (sf in strucfeats) {
      burial.columns <- which(grepl("rel.burial",colnames(sf)))
      if (length(burial.columns) > 0) {
        for (col in burial.columns) {
          prot <- sub("rel.burial.","",colnames(sf)[[col]])
          td$add.track(sf[,col],prot,"orange",maxVal=1)
        }
      }
    }
  }
  
  cat("Drawing genophenogram...\n")
  
  #build genophenogram
  # img.width <- length(wt.aa) * 0.06 + 2.5
  if (args$squish) {
    img.width <- 12
  } else {
    img.width <- length(wt.aa) * 0.13 + 4
  }
  img.height <- 4.5 + 0.13 * if(is.null(td)) 0 else td$num.tracks()
  
  pdf(pdffile,width=img.width,height=img.height)
  genophenogram(
    wt.aa,
    data$start,
    data$variant,
    data$score,
    1,0,
    error=data$se,
    grayBack=TRUE,
    img.width=img.width,
    tracks=td
  )
  invisible(dev.off())
  cat("done\n")

}

cat("\nScript completed successfully!\n")



