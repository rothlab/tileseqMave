

differenceMap <- function(condMapFile,refMapFile,outFile,pdfFile) {
	
	op <- options(stringsAsFactors=FALSE)
	
	library(mavevis)
	library(hgvsParseR)
	library(yogitools)
	library(hash)

	cat("Reading and parsing input data...")
	condMap <- read.csv(condMapFile)
	refMap <- read.csv(refMapFile)
	
	condIdx <- hash(condMap$hgvs_pro,1:nrow(condMap))
	refIdx <- hash(refMap$hgvs_pro,1:nrow(refMap))

	deltaMap <- to.df(do.call(rbind,lapply(union(keys(condIdx),keys(refIdx)),function(mut) {
		if (!has.key(mut,condIdx) || !has.key(mut,refIdx)) {
			return(NULL)
		}
		ic <- condIdx[[mut]]
		ir <- refIdx[[mut]]
		condScore <- condMap[ic,"score"]
		refScore <- refMap[ir,"score"]
		condSD <- condMap[ic,"sd"]
		refSD <- refMap[ir,"sd"]
		condDF <- with(condMap[ic,],((1/se)*sd)^2)
		refDF <- with(refMap[ir,],((1/se)*sd)^2)
		deltaSD <- sqrt(condSD^2 + refSD^2)
		list(
			hgvs_pro=condMap[ic,"hgvs_pro"],
			condScore=condScore,
			condSD=condSD,
			refScore=refScore,
			refSD=refSD,
			deltaScore=condScore-refScore,
			deltaSD=deltaSD,
			deltaSE=deltaSD/sqrt(condDF+refDF)
		)
	})))

	write.csv(deltaMap,outFile,row.names=FALSE)

	deltaMap <- cbind(parseHGVS(deltaMap$hgvs_pro,aacode=1),deltaMap[,-1])

	ancestrals <- with(deltaMap,tapply(ancestral,start,unique))
	wt.aa <- ancestrals[as.character(1:max(deltaMap$start))]
	wt.aa[[1]] <- "M"


	img.width <- length(wt.aa) * 0.13 + 4
	img.height <- 4.5 

	pdf(pdfFile,width=img.width,height=img.height)
	genophenogram(
		wt.aa,
		deltaMap$start,
		deltaMap$variant,
		deltaMap$deltaScore,
		0,-0.5,
		error=deltaMap$deltaSE,
		grayBack=TRUE,
		img.width=img.width
	)
	invisible(dev.off())
	
	options(op)
	
}


#refMapFile <- "/home/jweile/projects/tileseqMave/workspace/output/CBS/merged/lowB6/CBS_low_imputation_refined_mavedb.csv"
#condMapFile <- "/home/jweile/projects/tileseqMave/workspace/output/CBS/merged/highB6/CBS_high_imputation_refined_mavedb.csv"
#outFile <- "/home/jweile/projects/tileseqMave/workspace/output/CBS/merged/highB6/CBS_delta.csv"
#pdfFile <- "/home/jweile/projects/tileseqMave/workspace/output/CBS/merged/highB6/CBS_delta.pdf"
#differenceMap(condMapFile,refMapFile,outFile,pdfFile)#
