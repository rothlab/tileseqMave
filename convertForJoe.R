#!/usr/bin/Rscript

#################################################################################################
# convertForJoe.R converts a MaveDB file into Joe's imputation input format.
#
# Usage: Rscript convertForJoe.R in=<path-to-input-file> out=<path-to-output-file> df=<num>
# If no output parameter is provided, one will be generated automatically in the same 
# folder as the input file.
# The "df" parameters represents the degrees of freedom to uses. If not provided, it defaults to 2.
#
# DEPENDENCIES: yogitools and hgvsParseR
# To install these dependencies use devtools:
# install.packages(devtools)
# devtools::install_github("jweile/yogitools")
# devtools::install_github("jweile/hgvsParseR")
##################################################################################################

options(stringsAsFactors=FALSE)
library(yogitools)
library(hgvsParseR)

#Read command line arguments
infile <- getArg("in",required=TRUE)
if (!file.exists(infile)) {
	stop(infile," does not exist!")
}
outfile <- getArg("out",default=sub("\\.csv$","_joe.txt",infile))
df <- as.integer(getArg("df",default=2))
if (is.na(df)) {
	stop("df must be an integer!")
}

#Read input file
cat("Reading input file...")
indata <- read.csv(infile)
cat("done\n")

#Check that all required columns are present
required <- c("hgvs_pro","score","sd","se")
if (!all(required %in% colnames(indata))) {
	missing <- setdiff(required,colnames(indata))
	stop("Input file ",infile," is missing the following column(s): ",paste(missing,collapse=", "))
}

#parse HGVS variant descriptor strings
cat("Parsing HGVS variant descriptors...")
vardata <- parseHGVS(indata$hgvs_pro,aacode=1)
vardata[which(vardata$variant=="*"),"type"] <- "nonsense"
cat("done\n")

#Check to make sure the HGVS strings were all valid
if (any(vardata$type == "invalid")) {
	#If not, complain!
	warning(sprintf(
		"There were %d invalid HGVS strings!",
		sum(vardata$type == "invalid")
	))
}

#Build output table
cat("Building output table...")
outdata <- data.frame(
	aa_pos=vardata$start,
	aa_ref=vardata$ancestral,
	aa_alt=apply(vardata, 1, function(row) {
		switch(row[["type"]],
			synonymous = row[["ancestral"]],
			nonsense = "_",
			substitution = row[["variant"]],
			default = stop("Unrecognized variant type!")
		)
	}),
	#we intentionally use the inverse of standard error as the quality score here instead of
	# quality_score=indata$phiPrime,
	quality_score=1/indata$se,
	num_replicates=df,
	fitness_input=indata$score,
	fitness_input_sd=indata$sd
)
cat("done\n")

#write results to file
cat("Writing results to file...")
write.table(outdata,outfile,sep="\t",row.names=FALSE,quote=FALSE)
cat("done\n")

cat("\nScript completed successfully!\n")
