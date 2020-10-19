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
# This is a command line wrapper for scoring
#####################################################

options(
	stringsAsFactors=FALSE,
	ignore.interactive=TRUE
)
library(tileseqMave)
library(yogitools)
library(argparser)

#process command line arguments
p <- arg_parser(
	"Contamination diagnostic",
	name="contamination.R"
)
p <- add_argument(p, "dataDir", help="workspace data directory")
p <- add_argument(p, "--parameters", help="parameter file. Defaults to parameters.json in the data directory.")
args <- parse_args(p)

#ensure datadir ends in "/" and exists
dataDir <- args$dataDir
if (!grepl("/$",dataDir)) {
	dataDir <- paste0(dataDir,"/")
}
if (!dir.exists(dataDir)) {
	#logger cannot initialize without dataDirectory, so just a simple exception here.
	stop("Data folder does not exist!")
}
paramfile <- if (is.na(args$parameters)) paste0(dataDir,"parameters.json") else args$parameters


params <- parseParameters(paramfile)
#find counts folder
latest <- latestSubDir(parentDir=dataDir,pattern="_mut_call$|mut_count$")
mcDir <- latest[["dir"]]
qcDir <- paste0(latest[["label"]],latest[["timeStamp"]],"_QC")
dir.create(qcDir,showWarnings=FALSE)

marg <- read.csv(paste0(mcDir,"/marginalCounts.csv"))
marg$pos <- as.integer(gsub("\\D+","",marg$aaChange))
marg$tile <- params$pos2tile(marg$pos)

depths <- read.csv(paste0(mcDir,"/sampleDepths.csv"))
depths$margsum <- sapply(1:nrow(depths),function(i) with(depths[i,],{
	cn <- sprintf("%s.t1.rep%d.count",Condition,Replicate)
	rs <- which(marg$tile==Tile.ID)
	sum(marg[rs,cn])
}))

# depths <- depths[depths$Condition=="nonselect",]
pdf(paste0(qcDir,"/contamDiagn.pdf"),8.5,11)
op <- par(mfrow=c(3,2))
# tcID <- with(depths,paste0(Condition,Tile.ID,"t",Time.point))
plotcols <- rainbow(max(params$tiles[,1]))
invisible(lapply(params$conditions$names, function(cond) {
	subdepths <- depths[depths$Condition==cond,]
	with(subdepths,{
		plot(depth,margsum,
			xlab="Read depth",ylab=expression(Sigma~"Marginal Counts"),pch=20,
			main=cond
		)
	})
	tapply(1:nrow(subdepths),subdepths$Tile.ID,function(rs) {
		coords <- subdepths[rs,c("depth","margsum")]
		tile <- subdepths[rs[[1]],"Tile.ID"]
		lines(coords[,1],coords[,2],col=plotcols[tile])
		text(mean(coords[,1]),mean(coords[,2]),tile,col=plotcols[tile])
	})
}))
par(op)
invisible(dev.off())

cat("Done!\n")
