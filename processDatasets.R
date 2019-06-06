
options(stringsAsFactors=FALSE)
library(yogilog)
library(mavevis)
library(yogitools)
library(tileseqMave)

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

inputFile <- getArg("datasets",default="workspace/input/datasets.csv")
datasets <- read.csv(inputFile)

pseudoObservations <- as.integer(getArg("pseudoObservations",default=2))
conservativeMode <- as.logical(getArg("conservativeMode",default=TRUE))

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

