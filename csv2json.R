#!/usr/bin/Rscript
options(stringsAsFactors=FALSE)

#load libraries
library(tileseqMave)
library(argparser)

#process command line arguments
p <- arg_parser(
	"Checks TileSeq parameter CSV file for validity and converts it to JSON format.",
	name="csv2json.R"
)
p <- add_argument(p, "infile", help="input file")
p <- add_argument(p, "--outfile", help="output file")
args <- parse_args(p)
outfile <- if (is.na(args$outfile)) NULL else args$outfile

csvParam2Json(args$infile,outfile)
