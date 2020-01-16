
#' Translate a single HGVS string from nucleotide to amino acid level
#' 
#' @param hgvs a single cds-level HGVS string. May contain multiple mutations in-cis.
#' @param params a parameter object with CDS definitions
#' @param builder an optional protein HGVS builder object, if none is provided, a new one is created.
#' @param strictMode whether HGVS will be validated against the reference sequence before translation
#' @return translated HGVS string
#' @export
translateHGVS <- function(hgvs, params, 
	builder=new.hgvs.builder.p(aacode=3),
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
	if (any(breakdown$type %in% c("duplication","inversion","conversion","amplification"))) {
		stop("Unsupported variant type:",hgvs)
	}

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
					with(breakdown[i,],(end-start+1+nchar(breakdown[i,"variant"]))%%3==0)
				}
			)
		})
		frameshifts <- fsCandidates[!inFrame]
	} else {
		frameshifts <- integer(0)
	}

	#if any frameshifts were found, describe the first of them and we're done.
	if (length(frameshifts) > 0) {
		fsStart <- min(breakdown[frameshifts,"start"])
		codonIdx <- (fsStart -1) %/%3 + 1
		aa <- trtable[[codons[[codonIdx]]]]
		return(builder$frameshift(codonIdx,aa))
	}

	#if there are no frameshifts, then we have a LOT more work to do...

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
		wtaa <- trtable[[codons[[codonIdx]]]]
		mutaa <- if (codon == "") {
			"-"
		} else if (nchar(codon) > 3) {
			paste(sapply(seq(1,nchar(codon),3),function(i) trtable[[substr(codon,i,i+2)]]),collapse="")
		} else {
			trtable[[codon]]
		}
		list(pos=codonIdx,wtaa=wtaa,mutaa=mutaa)
	}))

	#Check for runs in the amino acid changes and group them accordingly
	aaChanges <- aaChanges[order(aaChanges$pos),]
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
	return(finalHGVS)

}

