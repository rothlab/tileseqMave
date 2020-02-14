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
		replicate="#Replicate",timepoint="#Timepoint",depth="#Final read-depth"
	)
	#otherwise throw error
	if (any(!(requiredFields %in% names(header)))) {
		missingFields <- requiredFields[which(!(requiredFields %in% names(header)))]
		stop(filename," is missing header field(s): ",paste(missingFields,collapse=", "))
	}
	#prepare metdata object
	metadata <- lapply(requiredFields,function(f) header[[f]])

	#parse main table
	countTable <- read.csv(textConnection(lines),comment.char="#")
	#and attach metadata
	attributes(countTable) <- c(attributes(countTable),metadata)

	options(op)

	return(countTable)
}

#' Translate a single HGVS string from nucleotide to amino acid level
#' 
#' This function takes a codon-level HGVS string and translates it to amino acid level.
#' It does so in two different ways: A joint view, that describes the overall effect of all mutations
#' on the protein, as well as a codon-wise segmented way that describes changes on individual 
#' codons/amino-acids separately, so that they are more easily separable for Marginal Frequency analysis
#' 
#' @param hgvs a single cds-level HGVS string. May contain multiple mutations in-cis.
#' @param params a parameter object with CDS definitions
#' @param builder an optional protein HGVS builder object, if none is provided, a new one is created.
#' @param strictMode whether HGVS will be validated against the reference sequence before translation
#' @return a named vector containing the following strings: hgvsp: the translated HGVS string; 
#'   codonChanges: A simple string expression of the involved codon changes; codonHGVS: A codon-wise
#'   segmented HGVS string; aaChanges: A simple amino acid change string; and aaChangeHGVS: A codon-wise
#'   segmented HGVS string.
#' @export
translateHGVS <- function(hgvs, params, 
	builder=new.hgvs.builder.p(aacode=3),cbuilder=new.hgvs.builder.c(),
	strictMode=TRUE) {

	library(yogitools)
	library(hgvsParseR)

	if (!inherits(hgvs,"character") && length(hgvs) != 1) {
		stop("Input must be single HGVS string!")
	}
	if (!grepl("^c\\.",hgvs)) {
		stop("HGVS string must be at CDS level!")
	}

	#extract coding sequence and break into codons
	cdsSeq <- params$template$cdsSeq
	codonStarts <- seq(1,nchar(cdsSeq),3)
	# codonIndices <- sapply(codonStarts, function(i) c(i,i+1,i+2))
	codons <- sapply(codonStarts,function(i) substr(cdsSeq,i,i+2))

	#load DNA translation table
	data(trtable)

	#parse the DNA HGVS string
	breakdown <- parseHGVS(hgvs)
	#not every breakdown will have an end column, but we rely on it below
	#so we add a fake one if it's missing.
	if (!("end" %in% colnames(breakdown))) {
		breakdown$end <- NA
	}

	#check for unsupported mutation types
	unsupportedTypes <- c("duplication","inversion","conversion","amplification")
	if (any(breakdown$type %in% unsupportedTypes)) {
		offenders <- with(breakdown,unique(type[which(type %in% unsupportedTypes)]))
		stop("Unsupported variant type(s): ",paste(offenders, collapse=","), " in ", hgvs)
	}

	#########################
	# CHECK FOR FRAMESHIFTS #
	#########################

	#First, check for frameshifts. 
	fsCandidates <- which(breakdown$type %in% c("singledeletion","deletion","insertion","delins"))
	if (length(fsCandidates) > 0) {
		inFrame <- sapply(fsCandidates,function(i) {
			switch(
				breakdown[i,"type"],
				singledeletion = FALSE,
				deletion = {
					with(breakdown[i,],(end-start+1)%%3==0)
				},
				insertion = {
					nchar(breakdown[i,"variant"])%%3==0
				},
				delins = {
					with(breakdown[i,],(nchar(breakdown[i,"variant"])-(end-start+1))%%3==0)
				}
			)
		})
		frameshifts <- fsCandidates[!inFrame]
	} else {
		frameshifts <- integer(0)
	}

	#if any frameshifts were found, describe the first of them and we're done.
	if (length(frameshifts) > 0) {
		first <- frameshifts[[which.min(breakdown[frameshifts,"start"])]]
		fsStart <- breakdown[first,"start"]
		codonIdx <- (fsStart -1) %/%3 + 1
		aa <- trtable[[codons[[codonIdx]]]]
		fsout <- builder$frameshift(codonIdx,aa)
		fsSimple <- paste0(aa,codonIdx,"fs")
		fsCodon <- paste0(codons[[codonIdx]],codonIdx,"indel")
		return(c(hgvsp=fsout,codonChanges=fsCodon,codonHGVS=breakdown[first,"hgvs"],aaChanges=fsSimple,aaChangeHGVS=fsout))
	}

	#if there are no frameshifts, then we have a LOT more work to do...

	########################################
	# APPLY CHANGES TO EACH AFFECTED CODON #
	########################################

	#list affected codons of each mutation to identify potential overlap
	affectedCodons <- lapply(1:nrow(breakdown),function(i) {
		ncs <- if (!is.na(breakdown[i,"end"])) {   
			seq(breakdown[i,"start"],breakdown[i,"end"])
		} else {
			breakdown[i,"start"]
		}
		unique((ncs -1) %/%3 + 1)
	})

	#calculate cumulative AA changes by iterating over each codon affected at least once
	aaChanges <- as.df(lapply(Reduce(union,affectedCodons), function(codonIdx) {
		#extract relevant codon sequence
		codon <- codons[[codonIdx]]
		cStart <- (codonIdx-1)*3+1
		#pull up the relevant mutation entries that affect this codon
		relevantMuts <- which(sapply(affectedCodons,function(l) codonIdx %in% l))
		#iterate over the relevant entries and cumulatively apply them to the codon
		for (mi in relevantMuts) {
			#calculate start and end indices with respect to whole CDS and codon
			mStart <- breakdown[mi,"start"]
			mEnd <- breakdown[mi,"end"]
			mStartInner <- (mStart-1)%%3+1
			mEndInner <- (mEnd-1)%%3+1
			insSeq <- breakdown[mi,"variant"]
			switch(
				breakdown[mi,"type"],
				substitution={
					if (strictMode) {
						wtNc <- breakdown[mi,"ancestral"]
						if (substr(codon,mStartInner,mStartInner) != wtNc) {
							stop("Reference mismatch!")
						}
					}
					substr(codon,mStartInner,mStartInner) <- insSeq
				},
				#Reminder: deletion / insertion lengths can only be multiple of 3
				#since they would have been caught in the frameshift check otherwise.
				deletion={
					#does the deletion leave a codon prefix at the beginning?
					if (mStart > cStart) {
						#FIXME: What if another mutation changes the downstream sequence?
						preSeq <- substr(codon,1,mStartInner-1)
						sufLen <- 3-nchar(preSeq)
						sufSeq <- substr(cdsSeq,mEnd+1,mEnd+sufLen)
						codon <- paste0(preSeq,sufSeq)
					} else {
						#otherwise it either covers the whole codon, 
						#or it has already been fused with an upstream codon
						#so we must delete this codon
						codon <- ""
					}
				},
				insertion={
					#double-check that start and end are right next to each other
					if (mEnd != mStart+1) {
						stop("Invalid insertion point!")
					}
					preSeq <- substr(codon,1,mStartInner)
					sufSeq <- substr(codon,mEndInner,3) 
					#this is now actually more than just one codon!
					codon <- paste0(preSeq,insSeq,sufSeq)
				},
				delins={
					#check if we're starting mid-codon
					if (mStart > cStart) {
						preSeq <- substr(codon,1,mStartInner-1)
						sufLen <- 3-nchar(preSeq)
						sufSeq <- substr(insSeq,1,sufLen)
						codon <- paste0(preSeq,sufSeq)
					} else if (mEnd > cStart+2) {
						#if this is not the last codon in the replacement area, we replace everything
						offset <- cStart-mStart
						codon <- substr(insSeq,offset+1,offset+3)
					} else {
						#this is the last codon, which may have a suffix, as well as leftover insertion sequence
						sufSeq <- substring(codon,mEndInner+1,3)
						offset <- cStart-mStart
						codon <- paste0(substr(insSeq,offset+1,nchar(insSeq)),sufSeq)
					}
				},
				{#any other type
					stop("Unsupported HGVS type:",breakdown[mi,"type"])
				}
			)
		}
		wtcodon <- codons[[codonIdx]]
		wtaa <- trtable[[wtcodon]]
		mutaa <- if (codon == "") {
			"-"
		} else if (nchar(codon) > 3) {
			paste(sapply(seq(1,nchar(codon),3),function(i) trtable[[substr(codon,i,i+2)]]),collapse="")
		} else {
			trtable[[codon]]
		}
		list(pos=codonIdx,wtaa=wtaa,mutaa=mutaa,wtcodon=wtcodon,mutcodon=codon)
	}))
	
	#sort aa changes by position
	aaChanges <- aaChanges[order(aaChanges$pos),]

	#######################################################
	# Build codon-centric change strings for later output #
	#######################################################
	#Helper function to build codon-specific HGVS
	codonWiseHGVS <- function(wt,pos,mut) {
		cstart <- pos*3-2
		diffs <- sapply(1:3,function(i)substr(wt,i,i)!=substr(mut,i,i))
		ndiff <- sum(diffs)
		if (mut == "") {
			return(cbuilder$deletion(cstart,cstart+2))
		} else if (nchar(mut) == 3 && ndiff == 1) { #SNV
			offset <- which(diffs)
			wtbase <- substr(wt,offset,offset)
			mutbase <- substr(mut,offset,offset)
			snvpos <- cstart+offset-1
			return(cbuilder$substitution(snvpos,wtbase,mutbase))
		} else { 
			return(cbuilder$delins(cstart,cstart+nchar(mut)-1,mut))
		}
	}
	#Helper function to single-AA specific HGVS
	aaWiseHGVS <- function(wt,pos,mut) {
		if(mut=="-") {
			return(builder$deletion(pos,wt,pos,wt))
		} else if (nchar(mut) > 1) {
			return(builder$delins(pos,wt,pos,wt,toChars(mut)))
		} else if (wt==mut) {
			return(builder$synonymous(pos,wt))
		} else {
			return(builder$substitution(pos,wt,mut))
		}
	}
	#simple codon change description string
	codonChangeStr <- paste(with(aaChanges,paste0(wtcodon,pos,mutcodon)),collapse="|")
	#HGVS-compliant codon change description string
	codonChangeHGVS <- if (nrow(aaChanges) > 1) {
		do.call(cbuilder$cis,with(aaChanges,mapply(codonWiseHGVS,wtcodon,pos,mutcodon,SIMPLIFY=FALSE)))
	} else {
		with(aaChanges,codonWiseHGVS(wtcodon,pos,mutcodon))
	}
	#simple AA change description string
	aaChangeStr <- paste(with(aaChanges,paste0(wtaa,pos,mutaa)),collapse="|")
	#HGVS-compliant AA change description string
	aaChangeHGVS <- if (nrow(aaChanges) > 1) {
		do.call(builder$cis,with(aaChanges,mapply(aaWiseHGVS,wtaa,pos,mutaa,SIMPLIFY=FALSE)))
	} else {
		with(aaChanges, aaWiseHGVS(wtaa,pos,mutaa))
	}

	##############################
	# BUILD OVERALL PROTEIN HGVS #
	##############################

	#Check for runs in the amino acid changes and group them accordingly
	if (nrow(aaChanges) > 1) {
		#chainHeads are mutations which are not immediately to the right of another one
		chainHeads <- which(c(TRUE,sapply(2:nrow(aaChanges),function(i) aaChanges$pos[[i]]-aaChanges$pos[[i-1]] > 1)))
		#find the head for each chain member
		assignedHead <- sapply(1:nrow(aaChanges),function(i) chainHeads[[max(which(chainHeads <= i))]])
		#break down the table according to chain assignment
		aaChangeGroups <- tapply(1:nrow(aaChanges),assignedHead,function(is) aaChanges[is,])
	} else {
		#if there is only one mutation, it's all one group :P
		aaChangeGroups <- list(aaChanges)
	}

	#iterate over groups and build HGVS strings
	aaHGVS <- lapply(aaChangeGroups, function(acs) {
		if (nrow(acs) > 1) {
			#this is a chain
			insSeq <- gsub("-","",paste0(acs$mutaa,collapse=""))
			fromPos <- acs$pos[[1]]
			fromAA <- acs$wtaa[[1]]
			toPos <- acs$pos[[nrow(acs)]]
			toAA <- acs$wtaa[[nrow(acs)]]
			#if it's all "-" it's a deletion
			if (insSeq == "") {
				builder$deletion(fromPos,fromAA,toPos,toAA)
			} else {
				builder$delins(fromPos,fromAA,toPos,toAA,toChars(insSeq))
			}
		} else {
			#singletons can be synonymous, deletions, delins or substitutions (which include nonsense)
			with(acs,{
				if (mutaa == wtaa) {
					builder$synonymous(pos,wtaa)
				} else if (mutaa == "-") {
					builder$deletion(pos,wtaa,pos,wtaa)
				} else if (nchar(mutaa) > 1) {
					builder$delins(pos,wtaa,pos,wtaa,toChars(mutaa))
				} else {
					builder$substitution(pos,wtaa,mutaa)
				}
			})
			
		}
	})

	#collapse in-cis multi-mutants
	finalHGVS <- if (length(aaHGVS) > 1) {
		do.call(builder$cis,aaHGVS)
	} else {
		aaHGVS[[1]]
	}

	#and finally, return the result
	return(c(
		hgvsp=finalHGVS,
		codonChanges=codonChangeStr,codonHGVS=codonChangeHGVS,
		aaChanges=aaChangeStr,aaChangeHGVS=aaChangeHGVS
	))

}

