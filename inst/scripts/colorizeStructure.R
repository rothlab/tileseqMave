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
# Pymol structure colorization script
#####################################################

options(
  stringsAsFactors=FALSE,
  ignore.interactive=TRUE
)

#load libraries
library(tileseqMave)
library(yogitools)
library(argparser)

# Process command line arguments------------------------------------------------

p <- arg_parser(
  "Create a PyMol script for colorizing structures according to a MAVE file",
  name="colorizeStructure.R"
)
p <- add_argument(p, "--workspace", help="tileseqMave workspace data directory. Defaults to current working directory")
p <- add_argument(p, "--input", help="either input file or input directory containing the count data. Defaults to subdirectory with latest timestamp ending in _scores")
p <- add_argument(p, "--output", help="either output pml file or output directory. Default adjusts automatically depending on input type")
p <- add_argument(p, "--chain", default="A", help="PDB chain identifier")
p <- add_argument(p, "--summary", default="median", help="Summary statistic to be used per AA position. Options: median, mean, min, max")
p <- add_argument(p, "--offset", default=0, help="PDB sequence offset. (For cases where the sequence numbering is shifted)")
args <- parse_args(p)

# Validation of parameters -----------------------------------------------------

#summary statistic function
summary <- match.arg(args$summary,choices=c("median","mean","min","max"))
fun <- switch(summary,median=median,mean=mean,min=min,max=max,default=median)

#PDB chain identifier
chain <- args$chain
if (!is.character(chain) || !grepl("^[A-Z]{1}$",chain)) {
  stop("'chain' must be single upper case letter!")
}

#Offset parameter
offset <- args$offset
if (!is.numeric(offset)) {
  stop("'offset' must be integer number!")
}

#Finding the input data and figuring out whether we're working with a single file
#or with a whole directory (and whether we are meant to find it automatically)
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

#create color map
cm <- yogitools::colmap()


# iterate over input files------------------------------------------------------
for (infile in infiles) {
  
  #if this is a single file (not from a directory)
  if (is.na(outdir)) {
    if (is.na(args$output)) {
      outfile <-  sub("\\.csv$","_colorize.pml",infile)
    } else {
      outfile <- args$output
      if (!grepl("pdf$",outfile)) {
        outfile <- paste0(outfile,".pdf")
      }
    }
  } else { #if this is one of many files in directory
    outfile <- paste0(outdir,"/",sub("\\.csv$","_colorize.pml",sub(".+/","",infile)))
  }
  
  # Core logic -----------------------------------------------------------------
  
  cat("Colorizing",infile,"\n")
  
  #read input
  indata <- read.csv(infile,comment.char="#")
  #filter out syn/stop
  indata <- indata[!grepl("Ter|=$",indata$hgvs_pro),]
  
  #extract positions
  pos <- as.integer(gsub("\\D","",indata$hgvs_pro))
  #calculate summary statistics per residue
  values <- tapply(indata$score,pos,fun,na.rm=TRUE)
  #corresponding position labels
  valpos <- as.integer(names(values))
  
  #generate script lines
  scriptLines <- sprintf("color %s , chain %s & resi %d",
    sub("#","0x",cm(values)), chain, valpos+offset
  )
  
  #write output
  cat("Writing output to",outfile,"\n")
  writeLines(scriptLines,con=outfile)
  
}

cat("Done!\n")



