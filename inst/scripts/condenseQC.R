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
# QC report condenser tool
#####################################################

options(
	stringsAsFactors=FALSE,
	ignore.interactive=TRUE
)

#load libraries
library(tileseqMave)
library(argparser)
# library(yogilog)

#process command line arguments
p <- arg_parser(
	"Condenser for PDF QC reports",
	name="condenseQC.R"
)
p <- add_argument(p, "--workspace", help="workspace data directory. Defaults to current working directory")
p <- add_argument(p, "--input", help="input directory containing the QC documents. Defaults to subdirectory with latest timestamp ending in _QC")
p <- add_argument(p, "--parameters", help="parameter file. Defaults to parameters.json in the data directory.")
p <- add_argument(p, "--srOverride", help="Manual override to allow singleton replicates. USE WITH EXTREME CAUTION!",flag=TRUE)
p <- add_argument(p, "--retainSingles", help="Retain individual files.",flag=TRUE)
p <- add_argument(p, "--allConditions", help="Indicate that libraryQC was run on all conditions (not just nonselect).",flag=TRUE)
args <- parse_args(p)

if (is.na(args$workspace)) {
  dataDir <- getwd()
} else {
  dataDir <- args$workspace
}
if (!grepl("/$",dataDir)) {
	dataDir <- paste0(dataDir,"/")
}
if (!dir.exists(dataDir)) {
	#logger cannot initialize without dataDirectory, so just a simple exception here.
	stop("Data folder does not exist!")
}
paramfile <- if (is.na(args$parameters)) paste0(dataDir,"parameters.json") else args$parameters

params <- parseParameters(paramfile,srOverride=args$srOverride)

qcDir <- args$input
if (is.na(qcDir)) {
  latest <- latestSubDir(parentDir=dataDir,pattern="_QC$")
  qcDir <- latest[["dir"]]
} else {
  if(!dir.exists(qcDir)) {
    stop("Input directory ",qcDir," does not exist!")
  }
}

nonSels <- getNonselects(params)
#if there are no nonselect conditions, then treat everything as nonselect
#(this is so we can process pure library QC runs)
if (length(nonSels) == 0 || args$allConditions) {
  nonSels <- params$conditions$names
}

cat("Condensing variant caller QC plots...")
#condense general sequencing QC plots
reportFiles <- sprintf("%s/%s.pdf",qcDir,c(
  "seqErrorRates","seqReads","effectiveSeqDepth","depthDrops"
))
#discard those that don't exist or are unreadable
reportFiles <- reportFiles[which(canRead(reportFiles))]
#if more than one file is applicable, merge them
if (length(reportFiles) > 0) {
  condensedFile <- sprintf("%s/varcallQC.pdf",qcDir)
  gsArgs <- c(
    "-dNOPAUSE",
    "-sDEVICE=pdfwrite",
    paste0("-sOUTPUTFILE='",condensedFile,"'"),
    "-dAutoRotatePages=/None",
    "-dBATCH",
    paste0("'",reportFiles,"'")
  )
  retVal <- system2("gs",gsArgs,stdout=FALSE)
  if (retVal == 0 && !args$retainSingles) {
    invisible(file.remove(reportFiles))
  }
}
fcount <- length(reportFiles)
if (fcount > 0) {
  cat(sprintf("%d files.\n",fcount))
} else {
  cat("no data!\n")
}


fcount <- 0
cat("Condensing library QC plots...")
for (nsCond in nonSels) {
  for (tp in params$timepoints$`Time point name`) {
  	#list report files in order
  	reportFiles <- sprintf("%s/%s_t%s_%s.pdf",qcDir,nsCond,tp,c(
  		"coverage","census","wellmeasured",
  		"tileRepCorr","WTlevels","jackpot",
  		"complexity","mutationtypes","nucleotide_bias","fsMap"
  	))
  	#discard those that don't exist or are unreadable
  	reportFiles <- reportFiles[which(canRead(reportFiles))]
  	#if more than one file is applicable, merge them
  	if (length(reportFiles) > 0) {
  		condensedFile <- sprintf("%s/%s_t%s_libraryQC.pdf",qcDir,nsCond,tp)
  		gsArgs <- c(
  			"-dNOPAUSE",
  			"-sDEVICE=pdfwrite",
  			paste0("-sOUTPUTFILE='",condensedFile,"'"),
  			"-dAutoRotatePages=/None",
  			"-dBATCH",
  			paste0("'",reportFiles,"'")
  		)
  		retVal <- system2("gs",gsArgs,stdout=FALSE)
  		if (retVal == 0 && !args$retainSingles) {
  			file.remove(reportFiles)
  		}
  	}
  	fcount <- fcount+length(reportFiles)
  }
}
if (fcount > 0) {
  cat(sprintf("%d files.\n",fcount))
} else {
  cat("no data!\n")
}


fcount <- 0
cat("Condensing selection QC plots...")
for (sCond in getSelects(params)) {
	for (tp in params$timepoints[,1]) {
		#list report files in order
		reportFiles <- sprintf("%s/%s_t%s_%s.pdf",qcDir,sCond,tp,c(
			"ns_replicates","phi_replicates","replicates","logPhiBias","codonCorr","filtering",
			"filtering2","synNonDiff","logPhiDistribution","errorModel","errorProfile"
		))
		#discard those that don't exist or are unreadable
		reportFiles <- reportFiles[which(canRead(reportFiles))]
		#if more than one file is applicable, merge them
		if (length(reportFiles) > 0) {
			condensedFile <- sprintf("%s/%s_t%s_selectionQC.pdf",qcDir,sCond,tp)
			gsArgs <- c(
				"-dNOPAUSE",
				"-sDEVICE=pdfwrite",
				paste0("-sOUTPUTFILE='",condensedFile,"'"),
				"-dAutoRotatePages=/None",
				"-dBATCH",
				paste0("'",reportFiles,"'")
			)
			retVal <- system2("gs",gsArgs,stdout=FALSE)
			if (retVal == 0 && !args$retainSingles) {
				file.remove(reportFiles)
			}
		}
	  fcount <- fcount+length(reportFiles)
	}
}
if (fcount > 0) {
  cat(sprintf("%d files.\n",fcount))
} else {
  cat("no data!\n")
}

cat("Done!\n")

