
options(stringsAsFactors=FALSE)
library(yogilog)
library(mavevis)
library(yogitools)

drawGenopheno <- function(outdir,uniprot) {
	scores <- read.csv(paste0(outdir,"mavedb_scores_perAA.csv"))
	wtseq <- getUniprotSeq(uniprot)
	mutations <- parseHGVS(scores$hgvs_pro,aacode=1)
	mutations$type[which(mutations$variant=="*")] <- "nonsense"

	wt.aa <- toChars(wtseq)
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
}

inputFile <- getArg("datasets",default="workspace/input/datasets_noCBS.csv")
datasets <- read.csv(inputFile)

for (i in 1:nrow(datasets)) {
	with(datasets[i,],{
		analyzeLegacyTileseqCounts(countfile,regionfile,outdir,
			logger=new.logger(paste0(outdir,"legacyTileSeq.log"))
		)
		drawGenopheno(outdir,uniprot)
	})
}


