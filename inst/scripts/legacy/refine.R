#!/usr/bin/Rscript

#################################################################################################
# refine.R performs refinement based on the output of Joe's imputation webtool
#
# Usage: Rscript refine.R in=<path-to-input-file> [org=<path-to-detailed-scores>] out=<path-to-output-file> pdf=<path-to-pdf-output>
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
library(hash)

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
orgfile <- getArg("org",default=NULL)
if (!is.null(orgfile)) {
	if (!file.exists(orgfile)) {
		stop(orgfile," does not exist!")
	}
}
outfile <- getArg("out",default=sub("\\.csv$","_refined_mavedb.csv",infile))
pdffile <- getArg("pdf",default=sub("\\.csv$","_refined_mavedb.pdf",infile))

#whether or not to use RMSD as the confidence gauge for imputed values
useRMSD <- as.logical(getArg("useRMSD",default=TRUE))

#cutoff for how many experimental values must exist in a column to include imputation output
columnCutoff <- as.integer(getArg("columnCutoff",default=3))

#Read input file
cat("Reading input file...")
indata <- read.csv(infile)
#Check that all required columns are present
required <- c(
	"aa_ref","aa_pos","aa_alt","num_replicates",
	"fitness","fitness_sd","fitness_imputed","fitness_imputed_se"
)
if (!all(required %in% colnames(indata))) {
	missing <- setdiff(required,colnames(indata))
	stop("Input file ",infile," is missing the following column(s): ",paste(missing,collapse=", "))
}

if (!is.null(orgfile)) {
	org <- read.csv(orgfile,comment.char="#")
	#Check that all required columns are present
	required <- c(
		"hgvs_pro","hgvs_nt","score","se"
	)
	if (!all(required %in% colnames(org))) {
		missing <- setdiff(required,colnames(org))
		stop("Input file ",orgfile," is missing the following column(s): ",paste(missing,collapse=", "))
	}
	hgvsLookup <- with(org,hash(hgvs_pro,strsplit(hgvs_nt," ")))
}


cat("done\n")
# with(indata,topoScatter(
# 	fitness_sd/sqrt(2),fitness_imputed_se,
# 	xlab="experiment se",ylab="imputation se",
# 	resolution=80
# ))
# abline(0,1,col="gray")

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

if (exists("hgvsLookup")) {
	hgvsc <- lapply(hgvss, function(h) {
		if (has.key(h,hgvsLookup)) {
			hgvsLookup[[h]]
		} else {
			"c.?"
		}
	})
}

#Determine which columns exceed the cutoff for imputation to be trusted
valsPerCol <- with(indata,tapply(fitness,aa_pos,function(xs)sum(!is.na(xs))))
trustedColumns <- as.integer(names(which(valsPerCol > columnCutoff)))
untrustedColumns <- as.integer(names(which(valsPerCol <= columnCutoff)))
if (length(untrustedColumns) > 0) {
	cat("Exluding imputation for columns",paste(untrustedColumns,collapse=", "),"not meeting cutoff.")
}



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
	impse <- if (is.na(fitness_imputed_se)) rmsd else fitness_imputed_se
	#if no imputation exists, or the column isn't trusted, use the experimental value (or NA)
	if (is.na(stretched) | !(aa_pos %in% trustedColumns)) {
		c(mj=fitness,sj=fitness_sd)
	#otherwise, if no experimental value exists, use imputation
	} else if (is.na(fitness)) {
		if (useRMSD) {
			c(mj=stretched,sj=rmsd)
		} else {
			c(mj=stretched,sj=impse)
		}
	#or, by default, form the confidence-weighted average between both
	} else {
		if (useRMSD) {
			out <- join.datapoints(ms=c(stretched,fitness),sds=c(rmsd,fitness_sd/sqrt(num_replicates)))
			c(mj=out[["mj"]],sj=out[["sj"]]*sqrt(num_replicates))
		} else {
			imp.se.adj <- impse*sqrt(10)/sqrt(num_replicates)
			exp.se <- fitness_sd/sqrt(num_replicates)
			out <- join.datapoints(ms=c(stretched,fitness),sds=c(imp.se.adj,exp.se))
			c(mj=out[["mj"]],sj=out[["sj"]]*sqrt(num_replicates))
		}
	}
})))
cat("done\n")

#build output table
numrep <- sapply(indata$num_replicates,function(x) if(is.na(x)) 1 else x)
outdata <- cbind(hgvss,refined,refined$sj/sqrt(numrep))
colnames(outdata) <- c("hgvs_pro","score","sd","se")

if (exists("hgvsc")) {
	outdata2 <- do.call(rbind,lapply(1:nrow(outdata), function(i) {
		cbind(hgvs_nt=hgvsc[[i]],outdata[i,])
	}))
}

cat("Writing output table...")
if (exists("outdata2")) {
	write.csv(outdata2,outfile,row.names=FALSE)
} else {
	write.csv(outdata,outfile,row.names=FALSE)
}
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
