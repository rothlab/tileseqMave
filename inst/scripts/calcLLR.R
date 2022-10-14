#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)

library(maveLLR)
library(argparser)
library(yogiroc)
library(yogilog)
library(pbmcapply)

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
p <- add_argument(p, "--outlierSuppression", help="Strength of outlier suppression to use. Between 0 and 1",default=0.0001)
p <- add_argument(p, "--logfile", help="The desired log file location.",default="llr.log")
p <- add_argument(p, "--iterations", help="The number of bootstrap iterations for determining LLR confidence interval",default=10000)
p <- add_argument(p, "--cores", help="The number of processes to run in parallel for multi-core processing. Warning: This also multiplies the amount of RAM used!",default=4)
p <- add_argument(p, "--printTransitions", help="Print evidence code transitions.",flag=TRUE)
args <- parse_args(p)
# args <- list(map="CHEK2_experimental.csv",reference="CHEK2_refVars_ClinvarPlusPlus.csv", outfile="b05/CHEK2_experimental_LLR_b05", bandwidth=0.5,kernel="epanechnikov",gauss=FALSE,spline=FALSE,posRange=NA,outlierSuppression=1,logfile="llr.log",iterations=10000)

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

logger$info("Reading input")
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

if (length(posScores) < 2 || length(negScores) < 2) {
  stop("Insufficient data!\nOnly ",length(posScores),
    " positive and ",length(negScores),
    " negative reference variant(s) observed in map!"
  )
}

logger$info("Calculating densties")
#build the LLR function
if (args$gauss) {
  llrObj <- buildLLR.gauss(posScores,negScores,args$spline)
} else {
  llrObj <- buildLLR.kernel(posScores,negScores,bw=args$bandwidth,kernel=args$kernel,outlierSuppression=args$outlierSuppression)
}

logger$info("Drawing LLR plot")
#draw the LLR function
pdf(paste0(outprefix,".pdf"),5,5)
drawDensityLLR(map$score,llrObj$llr,llrObj$posDens,llrObj$negDens,posScores,negScores)
invisible(dev.off())

#calculate all LLRs and their confidence intervals
logger$info("Bootstrapping confidence intervals at ",args$iterations," iterations")
mapLLR <- llrObj$llr(map$score) 
bootstrapScores <- Map(rnorm, args$iterations, map$score, map$se)
scoreLLR <- do.call("cbind", pbmclapply(bootstrapScores, llrObj$llr, mc.cores=args$cores))
left <- apply(scoreLLR, MARGIN=2, FUN=quantile, probs=0.025)
right <- apply(scoreLLR, MARGIN=2, FUN=quantile, probs=0.975)
out <- data.frame(map,llr=mapLLR,llrCI=sprintf("[%.03f;%.03f]",left,right))

#write result to file
outfile <- paste0(outprefix,".csv")
logger$info("Writing output to file ",outfile)
write.csv(out,outfile,row.names=FALSE)

#generate a PRC curve
logger$info("Drawing PRC curve")
yrobj <- yr2(
  truth=c(rep(TRUE,length(posScores)),rep(FALSE,length(negScores))),
  scores=data.frame(map=c(posScores,negScores)),high=FALSE
)
pdf(paste0(outprefix,"_prc.pdf"),7,7)

draw.prc.CI(yrobj,balanced=TRUE)
invisible(dev.off())

#print category transitions
if (args$printTransitions) {

  logger$info("Category transitions: ")
  optiLLR <- function(prior=0.1,posterior=0.9) {
    log10(posterior*(1-prior)/(prior*(1-posterior)))*4/3
  }
  llrThresholds <- function(LLRpvst=optiLLR(0.1),X=2) {
    c(
      patho.vstrong=LLRpvst,
      patho.strong=LLRpvst/X,
      patho.moderate=LLRpvst/(X^2),
      patho.support=LLRpvst/(X^3),
      benign.support=-LLRpvst/(X^3),
      benign.strong=-LLRpvst/(X^1)
    )
  }
  #TODO: Export this from maveLLR instead of re-implementing abve
  llrTs <- llrThresholds()
  llrTs <- c(llrTs[1:4],none=0,llrTs[5:6])

  xs <- seq(min(map$score,na.rm=TRUE),max(map$score,na.rm=TRUE),length.out=1000)
  ys <- llrObj$llr(xs)
  layers <- sapply(ys,function(y) {
    if (y > 0) {
      head(names(which(llrTs < y)),1)
    } else if (y < 0) {
      tail(names(which(llrTs > y)),1)
    } else {
      "none"
    }
  })
  for (i in 2:length(layers)) {
    if (layers[[i]] != layers[[i-1]]) {
      cat(sprintf("%.03f : %s -> %s\n",xs[[i]],layers[[i-1]],layers[[i]]))
    }
  }

}


logger$info("Done!")
