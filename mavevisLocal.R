#!/usr/bin/Rscript
options(stringsAsFactors=FALSE)

library(mavevis)
library(hgvsParseR)
library(yogitools)

#get input file argument
infile <- getArg("infile",required=TRUE)
#and check that the file exists and can be read
if (!canRead(infile)) {
	stop("The given input file does not exist or cannot be read!")
}
if (!grepl("\\.csv$",infile)) {
	stop("The given input file is not a CSV file!")
}

#get the uniprot id argument and validate
uniprot <- getArg("uniprot",required=TRUE)
uniprotRX <- "^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$"
if (!grepl(uniprotRX,uniprot)) {
	stop("The given Uniprot Accession is invalid!")
}
#TODO: Validate internet connection

#get the PDB argument, parse and validate it.
pdbArg <- getArg("pdb",default=NULL)
if (!is.null(pdbArg)) {
	pdb <- strsplit(pdbArg,";")[[1]]
	pdbIds <- sapply(strsplit(pdb,"#"),`[[`,1)
	pdbChains <- sapply(strsplit(pdb,"#"),`[[`,2)

	if (!all(grepl("^[0-9][A-Za-z0-9]{3}$",pdbIds))) {
		stop("One or more of the given PDB IDs is invalid!")
	}
	if (!all(grepl("^[A-Z]{1}$",pdbChains))) {
		stop("One or more of the PDB chain identifiers is invalid!")
	}
}


pdffile <- getArg("pdfOut",default=NULL)
if (is.null(pdffile)) {
	pdffile <- sub("\\.csv",".pdf",infile)
}


cat("Reading and parsing input data...")
indata <- read.csv(infile)
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
	wt.aa <- yogitools::toChars(mavevis::getUniprotSeq(uniprot))
}

if (!is.null(pdbArg)) {
	td <- new.trackdrawer(length(wt.aa),nox=TRUE)

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
} else {
	td <- NULL
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



