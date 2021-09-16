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
# This tool exports counts data to MaveDB format
#####################################################

options(
  stringsAsFactors=FALSE,
  ignore.interactive=TRUE
)

#load libraries
library(tileseqMave)
library(argparser)
library(yogilog)

#process command line arguments
p <- arg_parser(
  "Exports count data in MaveDB-compatible format.",
  name="exportCounts.R"
)
p <- add_argument(p, "--workspace", help="workspace data directory. Defaults to current working directory")
p <- add_argument(p, "--input", help="either input file or input directory containing the count data. Defaults to subdirectory with latest timestamp ending in _mut_count")
p <- add_argument(p, "--output", help="either output csv file or output directory. Default adjusts automatically depending on input type")
p <- add_argument(p, "--parameters", help="parameter file. Defaults to parameters.json in the data directory.")
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

input <- args$input
if (is.na(input)) {
  latest <- latestSubDir(parentDir=dataDir,pattern="_mut_count$")
  input <- latest[["dir"]]
}

if (dir.exists(input)) {
 infiles <- list.files(input, full.names = TRUE, pattern = "Counts.csv$")
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

#iterate over input files
for (infile in infiles) {
   
  #if this is a single file (not from a directory)
  if (is.na(outdir)) {
    if (is.na(args$output)) {
      outfile <-  sub("\\.csv$","_mavedb.csv",infile)
    } else {
      outfile <- args$output
      if (!grepl("\\.csv$",outfile)) {
        outfile <- paste0(outfile,".csv")
      }
    }
  } else { #if this is one of many files in directory
    outfile <- paste0(outdir,"/",sub("\\.csv$","_mavedb.csv",basename(infile)))
  }

  cat("Processing ",infile,"\n")

  #read the input file
  data <- read.csv(infile,comment.char="#")
  cn <- colnames(data)
  cn[cn=="hgvsc"] <- "hgvs_nt"
  cn[cn=="hgvsp"] <- "hgvs_pro"
  colnames(data) <- cn
  toDelete <- which(cn %in% c("codonChanges","codonHGVS","aaChanges","aaChangeHGVS"))
  data <- data[,-toDelete]

  write.csv(data,outfile,row.names=FALSE)

}

cat("\nScript completed successfully!\n")


