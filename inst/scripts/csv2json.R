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
#This is a command line wrapper for CSV2JSON
##############################################

options(stringsAsFactors=FALSE)

#load libraries
library(tileseqMave)
library(argparser)
library(yogilog)

#process command line arguments
p <- arg_parser(
	"Checks TileSeq parameter CSV file for validity and converts it to JSON format.",
	name="csv2json.R"
)
p <- add_argument(p, "infile", help="input file")
p <- add_argument(p, "--outfile", help="output file. Defaults to parameters.json in the same directory.")
p <- add_argument(p, "--logfile", help="log file. Defaults to csv2json.log in the same directory")
args <- parse_args(p)
outfile <- if (is.na(args$outfile)) sub("[^/]+$","parameters.json",args$infile) else args$outfile
logfile <- if (is.na(args$logfile)) sub("[^/]+$","csv2json.log",args$infile) else args$logfile

#set up logger and shunt it into the error handler
logger <- new.logger(logfile)
registerLogErrorHandler(logger)

#run the actual function
invisible(
	csvParam2Json(args$infile,outfile,logger)
)
