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

##############################################
#This is a command line tool for drawing a PRC curve against a reference set
##############################################


options(stringsAsFactors=FALSE)

library(argparser)
library(yogiroc)
library(yogilog)
library(yogitools)

#process command line arguments
p <- arg_parser(
  "Draw a PRC curve against a reference set",
  name="drawPRC.R"
)
p <- add_argument(p, "map", help="VE map file in MaveDB format")
p <- add_argument(p, "reference", help="reference variant set (from referenceSets.R")
p <- add_argument(p, "--predictors", help="comma-separated list of files with other predictors (in MaveDB format)")
p <- add_argument(p, "--predictorNames", help="comma-separated list of names for the above predictors")
p <- add_argument(p, "--predictorOrders", help="comma-separated list letters 'a' for ascending or 'd' for whether predictor scores are ascending or descending.")
p <- add_argument(p, "--outfile", help="The desired prefix for the output file name.")
p <- add_argument(p, "--posRanges", help="Positional ranges within the map to be plotted separately. E.g '1-24,25-66,67-")
p <- add_argument(p, "--logfile", help="The desired log file location.",default="prc.log")
args <- parse_args(p)
# args <- list(map="CHEK2_experimental.csv",reference="CHEK2_refVars_ClinvarPlusPlus.csv", outfile="b05/CHEK2_experimental_LLR_b05", bandwidth=0.5,kernel="epanechnikov",gauss=FALSE,spline=FALSE,posRange=NA,outlierSuppression=1,logfile="llr.log",iterations=10000)

#set up logger and shunt it into the error handler
logger <- new.logger(args$logfile)
registerLogErrorHandler(logger)

if (is.na(args$outfile)) {
  outprefix <- paste0(
    sub("\\.csv$","_",args$map),
    basename(sub("\\.csv$","_",args$reference)),
    "PRC"
  )
} else {
  outprefix <- args$outfile
}

#read reference set and check for proper format
refset <- read.csv(args$reference)
if (!all(c("hgvsp","referenceSet") %in% colnames(refset))) {
  stop("Reference file ",args$reference," is not a proper reference set file!")
}
refset$pos <- gsub("\\D+","",refset$hgvsp)|>as.integer()

#read map, check for proper format, index the table and add aa positions
map <- read.csv(args$map,comment.char="#")
if (!all(c("hgvs_pro","score") %in% colnames(map))) {
  stop("Map file ",args$map," is not following valid MaveDB format!")
}
map$pos <- gsub("\\D+","",map$hgvs_pro)|>as.integer()
rownames(map) <- map$hgvs_pro

maxPos <- max(map$pos,na.rm=TRUE)

if (!is.na(args$posRanges)) {
  posRanges <- strsplit(args$posRanges,",")[[1]] |> 
  lapply(function(xs)strsplit(xs,"-")[[1]]|>as.integer()) |>  
  lapply(function(xs) {
    if (length(xs) < 2) {
      xs <- c(xs,maxPos)
    } else if (is.na(xs[[1]])) {
      xs[[1]] <- 1
    }
    xs
  })
} else {
  posRanges <- list()
}

#also add global range to the end
posRanges <- c(posRanges,list(c(1,maxPos)))

#read predictor files
if (!is.na(args$predictors)) {
  predFiles <- strsplit(args$predictors,",")[[1]]
  #check if all the predictor files actually exist
  missing <- which(!sapply(predFiles,file.exists))
  if (length(missing) > 0) {
    stop("Predictor file(s) not found: ",paste(predFiles[missing],collapse=", "))
  }
  #read them, check them for proper format, index them and add aa positions
  predictors <- lapply(predFiles,function(pfile) {
    prd <- read.csv(pfile,comment.char="#")
    if (!all(c("hgvs_pro","score") %in% colnames(prd))) {
      stop("Predictor file ",pfile," is not following valid MaveDB format!")
    }
    prd$pos <- gsub("\\D+","",prd$hgvs_pro)|>as.integer()
    rownames(prd) <- prd$hgvs_pro
  })

  if (is.na(args$predictorNames)) {
    logger$warn("No predictor names provided. Using numbers instead...")
    prNames <- paste0("Predictor",seq_along(predictors))
  } else {
    prNames <- strsplit(args$predictorNames,",")[[1]]
    if (length(prNames) != length(predictors)) {
      stop("Predictor names list does not match number of predictors!")
    }
  }

  if (is.na(args$predictorOrders)) {
    logger$warn("No predictor orders provided, assuming ascending for each...")
    prOrders <- rep(TRUE,length(predictors))
  } else {
    prOrders <- strsplit(args$predictorOrders,",")[[1]]=="a"
    if (length(prOrders) != length(predictors)) {
      stop("Predictor orders list does not match number of predictors!")
    }
  }
} else {
  predictors <- list()
}


#iterate over ranges
pdf(paste0(outprefix,".pdf"),5,5)
lapply(posRanges, function(range) {

  rangeLabel <- sprintf("Position range %d-%d",range[[1]],range[[2]])
  refsubset <- refset[which(refset$pos >= range[[1]] & refset$pos <= range[[2]]),]
  scores <- map[refsubset$hgvsp,"score"]
  data <- data.frame(map=scores)
  dataOrder <- FALSE

  if (length(predictors) > 0) {
    prScores <- do.call(data.frame,lapply(predictors, function(prd) map[refsubset$hgvsp,"score"]))
    colnames(prScores) <- prNames
    data <- cbind(data,prScores)
    dataOrder <- c(dataOrder,prOrders)
  } 

  yr <- yogiroc::yr2(refsubset$referenceSet=="Positive", data, high=dataOrder)
  yogiroc::draw.prc.CI(yr,balanced=TRUE,main=rangeLabel)
  text(0,(ncol(data)+1)*7,sprintf("%d P/LP; %d B/LB",sum(refsubset$referenceSet=="Positive"),sum(refsubset$referenceSet=="Negative")),pos=4)
})
dev.off()|>invisible()

logger$info(paste0("Output file: ",outprefix,".pdf"))
logger$info("Done!")

