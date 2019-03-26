
options(stringsAsFactors=FALSE)
library(yogilog)
library(mavevis)

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

datasets <- as.data.frame(rbind(
	c(
		countfile = "workspace/input/HMGCR/rawData_S.txt",
		regionfile = "workspace/input/HMGCR/regions.txt",
		outdir = "workspace/output/HMGCR/",
		uniprot = "P04035"
	),
	c(
		countfile = "workspace/input/MTHFR/rawData_A222V_Q20_fol12.5.txt",
		regionfile = "workspace/input/MTHFR/regions.txt",
		outdir = "workspace/output/MTHFR/A222V/fol12/",
		uniprot = "P42898"
	),
	c(
		countfile = "workspace/input/MTHFR/rawData_A222V_Q20_fol25.txt",
		regionfile = "workspace/input/MTHFR/regions.txt",
		outdir = "workspace/output/MTHFR/A222V/fol25/",
		uniprot = "P42898"
	),
	c(
		countfile = "workspace/input/MTHFR/rawData_A222V_Q20_fol100.txt",
		regionfile = "workspace/input/MTHFR/regions.txt",
		outdir = "workspace/output/MTHFR/A222V/fol100/",
		uniprot = "P42898"
	),
	c(
		countfile = "workspace/input/MTHFR/rawData_A222V_fol200.txt",
		regionfile = "workspace/input/MTHFR/regions.txt",
		outdir = "workspace/output/MTHFR/A222V/fol200/",
		uniprot = "P42898"
	),
	c(
		countfile = "workspace/input/MTHFR/rawData_WT_Q20_fol12.5.txt",
		regionfile = "workspace/input/MTHFR/regions.txt",
		outdir = "workspace/output/MTHFR/WT/fol12/",
		uniprot = "P42898"
	),
	c(
		countfile = "workspace/input/MTHFR/rawData_WT_Q20_fol25.txt",
		regionfile = "workspace/input/MTHFR/regions.txt",
		outdir = "workspace/output/MTHFR/WT/fol25/",
		uniprot = "P42898"
	),
	c(
		countfile = "workspace/input/MTHFR/rawData_WT_Q20_fol100.txt",
		regionfile = "workspace/input/MTHFR/regions.txt",
		outdir = "workspace/output/MTHFR/WT/fol100/",
		uniprot = "P42898"
	),
	c(
		countfile = "workspace/input/MTHFR/rawData_WT_fol200.txt",
		regionfile = "workspace/input/MTHFR/regions.txt",
		outdir = "workspace/output/MTHFR/WT/fol200/",
		uniprot = "P42898"
	)
))

for (i in 1:nrow(datasets)) {
	with(datasets[i,],{
		analyzeLegacyTileseqCounts(countfile,regionfile,outdir,
			logger=new.logger(paste0(outdir,"legacyTileSeq.log"))
		)
		drawGenopheno(outdir,uniprot)
	})
}


