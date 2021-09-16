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

#' parse a count table
#' @param filename input file
#' @return data.frame with attributes containing header values
#' @export
parseCountFile <- function(filename) {
  op <- options(stringsAsFactors=FALSE)

  lines <- scan(filename,what="character",sep="\n",quiet=TRUE)
  #parse the header
  header <- do.call(c,lapply(
    strsplit(lines[grep("^#",lines)],":\\s*"),
    function(xs) setNames(xs[[2]],trimws(xs[[1]]))
  ))
  #Check that all required fields are present in the header
  requiredFields <- c(
    sample="#Sample",tile="#Tile",condition="#Condition",
    # replicate="#Replicate",timepoint="#Timepoint",depth="#Read-depth (pairs) after merging R1 and R2"
    replicate="#Replicate",timepoint="#Timepoint",
    raw.depth="#Raw read depth",
    depth="#Final read-depth",
    wtpairs="#Number of read pairs without mutations",
    mutpairs="#Total read pairs with mutations",
    unmapped="#Number of read pairs did not map to gene",
    mismapped="#Number of reads outside of the tile"
  )
  #otherwise throw error
  if (any(!(requiredFields %in% names(header)))) {
    missingFields <- requiredFields[which(!(requiredFields %in% names(header)))]
    stop(filename," is missing header field(s): ",paste(missingFields,collapse=", "))
  }
  #prepare metdata object
  metadata <- lapply(requiredFields,function(f) header[[f]])

  #if there are no mutations, return an empty table
  if (all(grepl("^#",lines))) {
    countTable <- data.frame(HGVS=character(),count=integer())
  } else {
    #parse main table
    countTable <- read.csv(textConnection(lines),comment.char="#")
  }
  #and attach metadata
  attributes(countTable) <- c(attributes(countTable),metadata)

  options(op)

  return(countTable)
}

#' parse coverage data file
#'
#' @param filename the input csv file
#' @param params the parameter sheet
#' @param tile the tile number associated with this file
#'
#' @return a vector detailing the number of rejected reads for each position
#' @export
#'
parseCoverageFile <- function(filename,params,tile) {
  
  coverageTable <- read.csv(filename)
  # coverageTable$m_either <- apply(coverageTable[,c("m_r1","m_r2")],1,max)
  coverageTable$rejected <- with(coverageTable,m_r1 + m_r2 - m_both - passed)
  
  tileStart <- params$tiles[params$tiles[,"Tile Number"]==tile,"Start NC in CDS"]
  tileEnd <- params$tiles[params$tiles[,"Tile Number"]==tile,"End NC in CDS"]
  posNames <- as.character(tileStart:tileEnd)

  rejections <- with(coverageTable,setNames(rejected,pos))
  rejections <- setNames(rejections[posNames],posNames)
  if (any(is.na(rejections))) {
    rejections[which(is.na(rejections))] <- 0
  }
  
  return(rejections)
}
