
options(stringsAsFactors=FALSE)
library(yogilog)
library(mavevis)
library(yogitools)
library(tileseqMave)
library(argparser)

#process command line arguments
p <- arg_parser(
	"Processes list of raw datasets using the legacy pipeline",
	name="processDatasets.R"
)
p <- add_argument(p, "datasets", help="input csv file containing list of datasets to process")
p <- add_argument(p, "--pseudoObservations", 
	help="number of pseudo-observations for Baldi&Long regularization", 
	default=2
)
p <- add_argument(p, "--conservative", help="toggle conservative mode",flag=TRUE)
args <- parse_args(p)

# inputFile <- getArg("datasets",default="workspace/input/datasets.csv")
inputFile <- args$datasets
if (!canRead(inputFile)) {
	stop("Unable to read input file",inputFile,"!")
}
datasets <- read.csv(inputFile)

# pseudoObservations <- as.integer(getArg("pseudoObservations",default=2))
pseudoObservations <- as.integer(args$pseudoObservations)
if (is.na(pseudoObservations)) {
	stop("pseudoObservations must be integer number!")
}
# conservativeMode <- as.logical(getArg("conservativeMode",default=TRUE))
conservativeMode <- args$conservative
if (is.na(conservativeMode)) {
	conservativeMode <- FALSE
}


flip <- function(xs) sapply(xs,function(x) if (is.na(x)) NA else if (x > 1) 1/x else x)

drawGenopheno <- function(outdir,uniprot) {
	scores <- read.csv(paste0(outdir,"mavedb_scores_perAA.csv"))
	flippedScores <- flip(scores$score)
	# wtseq <- getUniprotSeq(uniprot)
	mutations <- parseHGVS(scores$hgvs_pro,aacode=1)
	mutations$type[which(mutations$variant=="*")] <- "nonsense"
	
	ancestrals <- with(mutations,tapply(ancestral,start,unique))
	wt.aa <- ancestrals[as.character(1:max(mutations$start))]
	wt.aa[[1]] <- "M"
	# wt.aa <- toChars(wtseq)
	img.width <- length(wt.aa) * 0.06 + 2.5
	pdf(paste0(outdir,"genophenogram.pdf"),width=img.width,height=2.5)
	genophenogram(
		wt.aa,
		mutations$start,
		mutations$variant,
		scores$score,
		1,0,
		error=scores$se,
		grayBack=TRUE,
		img.width=img.width
	)
	invisible(dev.off())

	pdf(paste0(outdir,"genophenogram_flipped.pdf"),width=img.width,height=2.5)
	genophenogram(
		wt.aa,
		mutations$start,
		mutations$variant,
		flippedScores,
		1,0,
		error=scores$se,
		grayBack=TRUE,
		img.width=img.width
	)
	invisible(dev.off())
}

for (i in 1:nrow(datasets)) {
	with(datasets[i,],{
		#make sure the output directory exists
		dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
		#reset the log file and create logger
		logfile <- paste0(outdir,"legacyTileSeq.log")
		if (file.exists(logfile)) {
			file.remove(logfile)
		}
		logger <- new.logger(logfile)
		logger$info("\n\nProcessing",countfile)
		logger$info("Output directory:",outdir,"\n\n")
		#run analysis
		analyzeLegacyTileseqCounts(countfile,regionfile,outdir,
			inverseAssay=(uniprot=="Q96IV0"),
			logger=logger,
			pseudoObservations=pseudoObservations,
			conservativeMode=conservativeMode
		)
		#draw genophenogram
		drawGenopheno(outdir,uniprot)
	})
}

