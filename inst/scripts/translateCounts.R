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
#This is a command line wrapper for translateHGVS
##############################################

options(stringsAsFactors=FALSE)

#load libraries
library(tileseqMave)
library(argparser)
library(yogilog)

#process command line arguments
p <- arg_parser(
	"Translates HGVS strings in the count files and adds a frequency column.",
	name="translateCounts.R"
)
p <- add_argument(p, "dataDir", help="data directory")
args <- parse_args(p)

#get data directory and make sure it ends with a "/"
dataDir <- args$infile
if (!grepl("/$",dataDir)) {
	dataDir <- paste0(dataDir,"/")
}

#set up logger and shunt it into the error handler
logfile <- paste0(dataDir,"translateCounts.log")
logger <- new.logger(logfile)
registerLogErrorHandler(logger)

logger$info("Reading parameters")
paramFile  <- paste0(dataDir,"parameters.json")
params <- parseParameters(paramFile)


logger$info("Locating count files")
countfiles <- list.files(paste0(dataDir,"counts"),pattern="counts_sample\\d+\\.csv",full.names=TRUE)

invisible(lapply(countfiles, function(filename) {

	logger$info("Translating",filename)

	#read data
	counts <- parseCountFile(filename)
	#setup HGVS builder
	builder <- new.hgvs.builder.p(aacode=3)
	#run translation
	hgvsp <- sapply(foo$HGVS,translateHGVS,params,builder)
	#add result columns
	counts$HGVS_pro <- hgvsp
	counts$frequency <- counts$count/attr(counts,"depth")
	counts <- counts[,c("HGVS","HGVS_pro","count","frequency")]

	#re-build file header
	header <- mapply(
		function(name,value) sprintf("#%s: %s",name,value), 
		name=c("Sample","Tile","Condition","Replicate","Timepoint","Read-depth"),
		value=sapply(
			c("sample","tile","condition","replicate","timepoint","depth"),
			function(a)attr(counts,a)
		)
	)

	logger$info("Saving to file",filename)
	#overwrite count file with results
	con <- file(filename,open="w")
	writeLines(header,con)
	write.csv(rbind(multi,single),con,row.names=FALSE,quote=FALSE)
	close(con)

}))

logger$info("Done")
