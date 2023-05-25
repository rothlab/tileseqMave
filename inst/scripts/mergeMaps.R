#!/usr/bin/env Rscript

# Copyright (C) 2021  Jochen Weile, Roth Lab
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

options(stringsAsFactors=FALSE)

library(argparser)
library(yogitools)
library(yogilog)

#process command line arguments
p <- arg_parser(
  "Merge VE maps",
  name="mergeMaps.R",hide.opts=TRUE
)
p <- add_argument(p, "--maps", help="Input maps in MaveDB format",nargs=Inf)
p <- add_argument(p, "--out", help="Output file",default="merged_map.csv")
p <- add_argument(p, "--weightedAverage", help="Use inverse-variance weighted average",flag=TRUE)
p <- add_argument(p, "--logfile", help="The desired log file location.",default="referenceSets.log")

args <- parse_args(p)

#set up logger and shunt it into the error handler
logger <- new.logger(args$logfile)
registerLogErrorHandler(logger)

#check that input files exist and are readable
readable <- canRead(args$maps)
if (!all(readable)) {
  stop("Unable to read file(s): ",paste(args$maps[!readable],collapse=", "))
}

logger$info("Reading input...")

#read, validate and index input files
maps <- lapply(args$maps,function(mfile) {
  m <- read.csv(mfile,comment.char="#")
  if (!all(c("hgvs_pro","score","se","df") %in% colnames(m))) {
    stop("File ",mfile," is not a valid MaveDB file!")
  }
  rownames(m) <- m
  return(m)
})

allVars <- Reduce(union,lapply(maps,rownames))

logger$info("Merging...")

out <- as.df(lapply(allVars, function(v) {
  scores <- na.omit(sapply(maps,`[`,v,"score"))
  ses <- na.omit(sapply(maps,`[`,v,"se"))
  dfs <- na.omit(sapply(maps,`[`,v,"df"))
  sds <- ses*sqrt(dfs)
  ws <- if (args$weightedAverage) {
    vs <- sds^2
    (1/vs)/sum(1/vs)
  } else {
    rep(1,length(scores))
  }
  if (length(scores) == 1) { #then there's only one table with that variant
    return(list(hgvs_pro=v,score=scores,se=ses,dfs=df)
  } else {
    out <- weightedAverage(scores,ses,dfs,ws)
    return(list(hgvs_pro=hgvs_pro,score=out[["mj"]],se=out[["sj"]]/sqrt(out[["dfj"]]),df=out[["dfj"]]))
  }
}))

write.csv(out,args$out,colnames=FALSE)

logger$info("Done!")