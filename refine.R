#!/usr/bin/Rscript

#################################################################################################
# refine.R performs refinement based on the output of Joe's imputation webtool
#
# Usage: Rscript refine.R in=<path-to-input-file> out=<path-to-output-file> pdf=<path-to-pdf-output>
# If no output parameters are provided, they will be generated automatically in the same 
# folder as the input file.
#
# DEPENDENCIES: yogitools and hgvsParseR
# To install these dependencies use devtools:
# install.packages(devtools)
# devtools::install_github("jweile/yogitools")
# devtools::install_github("jweile/hgvsParseR")
# devtools::install_github("jweile/mavevis")
# The MaveVis dependency requires clustalo, dssp and freesasa to be installed as well!
##################################################################################################

options(stringsAsFactors=FALSE)
library(yogitools)
library(hgvsParseR)
library(mavevis)

#define helper function for joining datasets weighted by stdev
join.datapoints <- function(ms,sds) {
	ws <- (1/sds)/sum(1/sds)
	mj <- sum(ws*ms)
	vj <- sum(ws*(sds^2+ms^2)) -mj^2
	c(mj=mj,sj=sqrt(vj))
}
#define helper function to calculate RMSD
calcRMSD <- function(a,b) {
	sqrt(mean((a-b)^2,na.rm=TRUE))
}

#Read command line arguments
infile <- getArg("in",required=TRUE)
if (!file.exists(infile)) {
	stop(infile," does not exist!")
}
outfile <- getArg("out",default=sub("\\.csv$","_refined_mavedb.csv",infile))
pdffile <- getArg("pdf",default=sub("\\.csv$","_refined_mavedb.pdf",infile))

#Read input file
cat("Reading input file...")
indata <- read.csv(infile)
cat("done\n")

#Check that all required columns are present
required <- c("aa_ref","aa_pos","aa_alt","fitness","fitness_sd","fitness_imputed")
if (!all(required %in% colnames(indata))) {
	missing <- setdiff(required,colnames(indata))
	stop("Input file ",infile," is missing the following column(s): ",paste(missing,collapse=", "))
}

#Construct the HGVS strings
cat("Constructing HGVS descriptors...")
hgvsp <- new.hgvs.builder.p(aacode=3)
hgvss <- sapply(1:nrow(indata),function(i) with(indata[i,],{
	if (aa_alt=="_") {
		hgvsp$substitution(aa_pos,aa_ref,"*")
	} else if (aa_alt==aa_ref) {
		hgvsp$synonymous(aa_pos,aa_ref)
	} else {
		hgvsp$substitution(aa_pos,aa_ref,aa_alt)
	}
}))
cat("done\n")

cat("Stretching predicted values to [0:1] interval...")
#find the numerical interval that the imputation filled
imputedRange <- range(indata$fitness_imputed,na.rm=TRUE)
#stretch the predictions to the full [0;1] range
indata$stretched <- (indata$fitness_imputed - imputedRange[[1]])/(imputedRange[[2]]-imputedRange[[1]])
cat("done\n")

cat("Calculating prediction RMSD...")
#determine the RMSD of the prediction
rmsd <- with(indata,calcRMSD(stretched,fitness))
cat("done\n")

cat("Running refinement...")
#and refine based on that RMSD
refined <- as.df(lapply(1:nrow(indata),function(i) with(indata[i,],{
	if (is.na(stretched)) {
		c(mj=fitness,sj=fitness_sd)
	} else if (is.na(fitness)) {
		c(mj=stretched,sj=rmsd)
	} else {
		join.datapoints(ms=c(stretched,fitness),sds=c(rmsd,fitness_sd))
	}
})))
cat("done\n")


cat("Writing output table...")
#build output table
numrep <- sapply(indata$num_replicates,function(x) if(is.na(x)) 1 else x)
outdata <- cbind(hgvss,refined,refined$sj/sqrt(numrep))
colnames(outdata) <- c("hgvs_pro","score","sd","se")
write.csv(outdata,outfile,row.names=FALSE)
cat("done\n")

cat("Drawing genophenogram...")
#derive WT sequence
ancestrals <- with(indata,tapply(aa_ref,aa_pos,unique))
wt.aa <- ancestrals[as.character(1:max(indata$aa_pos))]
wt.aa[[1]] <- "M"

#build genophenogram
img.width <- length(wt.aa) * 0.06 + 2.5
pdf(pdffile,width=img.width,height=2.5)
genophenogram(
	wt.aa,
	indata$aa_pos,
	sapply(indata$aa_alt,function(x) if (x == "_") "*" else x),
	outdata$score,
	1,0,
	error=outdata$se,
	grayBack=TRUE,
	img.width=img.width
)
invisible(dev.off())
cat("done\n")


cat("\nScript completed successfully!\n")
