#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)

library(maveLLR)
library(argparser)
library(yogiroc)
library(yogilog)

#process command line arguments
p <- arg_parser(
  "Automatically assemble reference variant sets via Clinvar and Gnomad",
  name="calcLLR.R"
)
p <- add_argument(p, "map", help="VE map file in MaveDB format")
p <- add_argument(p, "reference", help="reference variant set (from referenceSets.R")
p <- add_argument(p, "--outfile", help="The desired prefix for the output file name.")
p <- add_argument(p, "--gauss",flag=TRUE, help="Use simple gaussians to determine llr instead of kernel density estimation")
p <- add_argument(p, "--spline",flag=TRUE, help="Use apply spline monotonization to gauss LLR function. (Only used if --gauss is used)")
p <- add_argument(p, "--bandwidth", help=paste("Kernel bandwidth. Only used when --gauss isn't used.",
                        "Valid options: Any number, or the following auto-selection methods:",
                        "'nrd0' (Silverman's rule of thumb), 'nrd' (Scott's rule of thumb), 'bcv' (biased cross validation),",
                        "'SJ' (Sheather and Jones), 'ucv' (unbiased cross validation)"),default="ucv")#default as string to avoid type enforcement
p <- add_argument(p, "--kernel", help=paste("Kernel type. Only used when --gauss isn't used.",
                        "Valid options: 'gaussian', 'epanechnikov', 'rectangular', 'triangular', 'biweight', 'cosine'",
                        "'optcosine','tricube','triweight','laplace','gamma','gamma_biased'"),default="epanechnikov")
p <- add_argument(p, "--posRange", help="Positional range within the map to be tested. Must be two integer numbers separated by dash, e.g. '1-189'")
p <- add_argument(p, "--logfile", help="The desired log file location.",default="llr.log")
args <- parse_args(p)

#set up logger and shunt it into the error handler
logger <- new.logger(args$logfile)
registerLogErrorHandler(logger)

#set outfile prefix
if (is.na(args$outfile)) {
  outprefix <- sub("\\.\\w+$","_llr",basename(args$reference))
} else {
  outprefix <- args$outfile
}

#handle bandwidth datatype conversion if necessary
if (!is.na(as.numeric(args$bandwidth))) {
  args$bandwidth <- as.numeric(args$bandwidth)
}

#handle positional range 
if (!is.na(args$posRange)) {
  posSplit <- strsplit(args$posRange,"-")[[1]]
  if (length(posSplit) == 2) {
    posSplit <- as.integer(posSplit)
  } else {
    stop("--posRange must be two integers separated by a dash")
  }
} else {
  posSplit <- NA
}

#import the map
map <- read.csv(args$map,comment.char="#")
rownames(map) <- map$hgvs_pro

#and filter for positional range if applicable
if (!any(is.na(posSplit))) {
  mpos <- as.integer(gsub("\\D+","",map$hgvs_pro))
  map <- map[which(mpos > posSplit[[1]] & mpos < posSplit[[2]]),]
}

#and import the reference set
refSets <- read.csv(args$reference)

#extract positive and negative reference sets
posRef <- unique(refSets[refSets$referenceSet=="Positive","hgvsp"])
negRef <- unique(refSets[refSets$referenceSet=="Negative","hgvsp"])
#and their scores
posScores <- na.omit(setNames(map[posRef,"score"],posRef))
negScores <- na.omit(setNames(map[negRef,"score"],negRef))

#build the LLR function
if (args$gauss) {
  llrObj <- buildLLR.gauss(posScores,negScores,args$spline)
} else {
  llrObj <- buildLLR.kernel(posScores,negScores,bw=args$bandwidth,kernel=args$kernel)
}

#draw the LLR function
pdf(paste0(outprefix,".pdf"),5,5)
drawDensityLLR(map$score,llrObj$llr,llrObj$posDens,llrObj$negDens,posScores,negScores)
dev.off()

#calculate all LLRs and their confidence intervals
mapLLR <- llrObj$llr(map$score)
left <- with(map,llrObj$llr(qnorm(.025,map$score,map$se)))
right <- with(map,llrObj$llr(qnorm(.975,map$score,map$se)))
out <- data.frame(map,llr=mapLLR,llrCI=sprintf("[%.03f;%.03f]",left,right))

#write result to file
outfile <- paste0(outprefix,".csv")
write.csv(out,outfile,row.names=FALSE)

#generate a PRC curve
yrobj <- yr2(
  truth=c(rep(TRUE,length(posScores)),rep(FALSE,length(negScores))),
  scores=data.frame(map=c(posScores,negScores)),high=FALSE
)
pdf(paste0(outprefix,"_prc.pdf"),7,7)
draw.prc.CI(yrobj,balanced=TRUE,names="map")
dev.off()



