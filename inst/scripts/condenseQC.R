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
p <- add_argument(p, "dataDir", help="workspace data directory")
p <- add_argument(p, "--parameters", help="parameter file. Defaults to parameters.json in the data directory.")
p <- add_argument(p, "--srOverride", help="Manual override to allow singleton replicates. USE WITH EXTREME CAUTION!",flag=TRUE)
p <- add_argument(p, "--retainSingles", help="Retain individual files.",flag=TRUE)
args <- parse_args(p)

dataDir <- args$dataDir
if (!grepl("/$",dataDir)) {
	dataDir <- paste0(dataDir,"/")
}
if (!dir.exists(dataDir)) {
	#logger cannot initialize without dataDirectory, so just a simple exception here.
	stop("Data folder does not exist!")
}
paramfile <- if (is.na(args$parameters)) paste0(dataDir,"parameters.json") else args$parameters

params <- parseParameters(paramfile,srOverride=args$srOverride)

latest <- latestSubDir(parentDir=dataDir,pattern="_QC$")
qcDir <- latest[["dir"]]

for (nsCond in getNonselects(params)) {
	#list report files in order
	reportFiles <- sprintf("%s/%s_%s.pdf",qcDir,nsCond,c(
		"coverage","census","wellmeasured","complexity","mutationtypes","nucleotide_bias"
	))
	#discard those that don't exist or are unreadable
	reportFiles <- reportFiles[which(canRead(reportFiles))]
	#if more than one file is applicable, merge them
	if (length(reportFiles) > 0) {
		condensedFile <- sprintf("%s/%s_libraryQC.pdf",qcDir,nsCond)
		gsArgs <- c(
			"-dNOPAUSE",
			"-sDEVICE=pdfwrite",
			paste0("-sOUTPUTFILE=",condensedFile),
			"-dBATCH",
			reportFiles
		)
		retVal <- system2("gs",gsArgs)
		if (retVal == 0 && !args$retainSingles) {
			file.remove(reportFiles)
		}
	}

}

for (sCond in getSelects(params)) {
	for (tp in params$timepoints[,1]) {
		#list report files in order
		reportFiles <- sprintf("%s/%s_t%s_%s.pdf",qcDir,sCond,tp,c(
			"ns_replicates","phi_replicates","replicates","filtering","logPhiDistribution","errorModel","errorProfile"
		))
		#discard those that don't exist or are unreadable
		reportFiles <- reportFiles[which(canRead(reportFiles))]
		#if more than one file is applicable, merge them
		if (length(reportFiles) > 0) {
			condensedFile <- sprintf("%s/%s_t%s_selectionQC.pdf",qcDir,sCond,tp)
			gsArgs <- c(
				"-dNOPAUSE",
				"-sDEVICE=pdfwrite",
				paste0("-sOUTPUTFILE=",condensedFile),
				"-dBATCH",
				reportFiles
			)
			retVal <- system2("gs",gsArgs)
			if (retVal == 0 && !args$retainSingles) {
				file.remove(reportFiles)
			}
		}
	}
}

