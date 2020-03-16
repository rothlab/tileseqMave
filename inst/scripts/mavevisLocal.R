#!/usr/bin/Rscript

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
p <- add_argument(p, "infile", help="input file. Must be CSV file in MaveDB format.")
p <- add_argument(p, "uniprot", help="Uniprot Accession for the underlying protein.")
p <- add_argument(p, "--pdb", help="PDB structures. Semicolon-separated list of pairings between PDB IDs and chain IDS.")
p <- add_argument(p, "--out", help="output pdf file. Defaults to the name of the input file with pdf extension.")
args <- parse_args(p)


#get input file argument
# infile <- getArg("infile",required=TRUE)
#and check that the file exists and can be read
if (!canRead(args$infile)) {
	stop("The given input file does not exist or cannot be read!")
}
if (!grepl("\\.csv$",args$infile)) {
	stop("The given input file is not a CSV file!")
}

#get the uniprot id argument and validate
# uniprot <- getArg("uniprot",required=TRUE)
uniprotRX <- "^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$"
if (!grepl(uniprotRX,args$uniprot)) {
	stop("The given Uniprot Accession is invalid!")
}
#TODO: Validate internet connection

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


# pdffile <- getArg("pdfOut",default=NULL)
pdffile <- if (is.na(args$out)) sub("\\.csv$",".pdf",args$infile) else args$out
# if (is.null(pdffile)) {
# 	pdffile <- sub("\\.csv",".pdf",args$infile)
# }


cat("Reading and parsing input data...")
indata <- read.csv(args$infile,comment.char="#")
if (!all(c("hgvs_pro","score") %in% colnames(indata))) {
	stop("Input file must be in MaveDB format!")
}
mutdata <- parseHGVS(indata$hgvs_pro,aacode=1)
mutdata[which(mutdata$variant=="*"),"type"] <- "nonsense"
data <- cbind(mutdata,indata[,-1])

cat("done!\n")

#derive WT sequence
cat("Deriving WT sequence...\n")
ancestrals <- with(data,tapply(ancestral,start,unique))
wt.aa <- ancestrals[as.character(1:max(data$start))]
wt.aa[[1]] <- "M"

if (any(is.na(wt.aa))) {
	warning("Unable to fully derive WT sequence! Defaulting to Uniprot sequence.")
	wt.aa <- yogitools::toChars(mavevis::getUniprotSeq(args$uniprot))
}

td <- new.trackdrawer(length(wt.aa),nox=TRUE)

if (!is.null(args$uniprot)) {
	cons <- calc.conservation(args$uniprot)
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
img.width <- length(wt.aa) * 0.13 + 4
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


cat("\nScript completed successfully!\n")



