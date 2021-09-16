#!/usr/bin/env Rscript

options(
  stringsAsFactors=FALSE,
  ignore.interactive=TRUE
)

library(yogitools)
library(tileseqMave)
library(argparser)
library(hgvsParseR)
library(pbmcapply)

#process command line arguments
p <- arg_parser(
  "Suggests synonymous / nonsense modes based on correlation with PROVEAN.",
  name="proveanNormalize.R"
)
p <- add_argument(p, "scorefile", help="tileseqMave complete score file")
p <- add_argument(p, "proveanfile", help="provean data file")
p <- add_argument(p, "--parameters", help="parameter file. Defaults to parameters.json")
p <- add_argument(p, "--output", help="output file")
p <- add_argument(p, "--selection", help="selection name",default="select")
p <- add_argument(p, "--precision", help="numerical precision in log(phi) space",default=0.02)
p <- add_argument(p, "--cores", default=6L, help="number of CPU cores to use in parallel for multi-threading")
args <- parse_args(p)

# args <- list(
#   scorefile="CHEK2v2_2020-05-11-17-36-33_scores/select_t1_complete.csv",
#   proveanfile="O96017_provean_output.tsv",
#   parameters="parameters.json",
#   output="provean",
#   selection="select",
#   precision=0.03,
#   cores=6L
# )

#parse parameter file
paramfile <- if (is.na(args$parameters)) paste0(dataDir,"parameters.json") else args$parameters
params <- parseParameters(paramfile)

#setup output name base
outbase <- if (is.na(args$output)) paste0(dataDir,"provean") else args$output

#read provean file
provean <- read.delim(args$proveanfile)

#build hgvs strings for provean variants
hgvsp <- new.hgvs.builder.p(aacode=3)
provean$hgvs <- sapply(1:nrow(provean), function(i) with(provean[i,],{
  if (RESIDUE_REF==RESIDUE_ALT) {
    hgvsp$synonymous(POSITION,RESIDUE_REF)
  } else {
    hgvsp$substitution(POSITION,RESIDUE_REF,RESIDUE_ALT)
  }
}))
#and index by hgvs
rownames(provean) <- provean$hgvs


#read complete score file
allscores <- read.csv(args$scorefile,comment.char="#")

#split aa change information for easy access
aac <- allscores$aaChange
allscores$fromAA <- substr(aac,1,1)
allscores$toAA <- substr(aac,nchar(aac),nchar(aac))
allscores$pos <- as.integer(substr(aac,2,nchar(aac)-1))

#find default modes for each region
defaultModes <- do.call(rbind,lapply(1:nrow(params$regions), function(regi) {

  #otherwise find all high-quality variants in this region
  rStart <- params$regions[regi,"Start AA"]
  rEnd <- params$regions[regi,"End AA"]
  scoresFiltered <- if (!all(is.na(allscores$logPhi.sd))) {
    with(allscores,allscores[which(is.na(filter) & pos >= rStart & pos <= rEnd & logPhi.sd < params$scoring$sdThreshold),])
  } else {
    with(allscores,allscores[which(is.na(filter)) & pos >= rStart & pos <= rEnd,])
  }

  #calculate medians and stdev for syn and non
  synMed <- with(scoresFiltered,median(
    logPhi[which(type == "synonymous")]
  ,na.rm=TRUE))
  synSD <- with(scoresFiltered,sd(
    logPhi[which(type == "synonymous")]
  ,na.rm=TRUE))
  nonMed <- with(scoresFiltered,median(
    logPhi[which(type == "nonsense")]
  ,na.rm=TRUE))
  nonSD <- with(scoresFiltered,sd(
    logPhi[which(type == "nonsense")]
  ,na.rm=TRUE))

  #apply any potential individual overrides #check if overrides are provided
  ovr <- getNormOverrides(params,args$selection,"1",regi)
  if (!is.na(ovr[["syn"]])) {
    synMed <- ovr[["syn"]]
  }
  if (!is.na(ovr[["non"]])) {
    nonMed <- ovr[["non"]]
  }

  return(c(synMed=synMed,synSD=synSD,nonMed=nonMed,nonSD=nonSD))

}))

#define the search space based on quantiles around the syn/non distributinos
ranges <- c(
  setNames(lapply(1:nrow(defaultModes),function(i) with(as.data.frame(defaultModes)[i,],{
    # seq(synMed-synSD/2,synMed+synSD/2,args$precision)
    seq(qnorm(0.20,synMed,synSD),qnorm(0.80,synMed,synSD),args$precision)
  })),paste0("syn",1:nrow(defaultModes))),
  setNames(lapply(1:nrow(defaultModes),function(i) with(as.data.frame(defaultModes)[i,],{
    # seq(nonMed-nonSD/2,nonMed+nonSD/2,args$precision)
    seq(qnorm(0.20,nonMed,nonSD),qnorm(0.80,nonMed,nonSD),args$precision)
  })),paste0("non",1:nrow(defaultModes)))
)
#add the medians themselves to the ranges
ranges <- mapply(function(a,b) sort(c(a,b)) ,a=ranges,b=c(defaultModes[,"synMed"],defaultModes[,"nonMed"]))


# searchSpace <- do.call(expand.grid,ranges)
searchSpace <- do.call(data.frame,lapply(ranges, sample, 10000,replace=TRUE))


#annotate each variant with its appropriate region
allscores$region <- sapply(allscores$pos, function(pos) {
  with(params$regions, `Region Number`[`Start AA` <= pos & `End AA` >= pos])
})
#filter out poorly measured variants
filterIdx <- with(allscores,which(
  is.na(filter) & 
  logPhi.sd < params$scoring$sdThreshold & 
  type=="missense"
))
filteredScores <- allscores[filterIdx,]

#helper function for weighed averaging across datapoints
join.datapoints <- function(ms,sds,dfs) {
  #weights
  ws <- (1/sds)/sum(1/sds)
  #weighted mean
  mj <- sum(ws*ms)
  #weighted joint variance
  vj <- sum(ws*(sds^2+ms^2)) -mj^2
  #return output
  return(c(mj=mj,sj=sqrt(vj),dfj=sum(dfs)))
}

##############
#Check default settings as baseline
##############

#Calculate scores for default settings
nregions <- nrow(params$regions)
jointScores <- do.call(rbind,tapply(1:nrow(filteredScores),filteredScores$hgvsp,function(js) {
  c(join.datapoints(
    filteredScores[js,"score"], filteredScores[js,"score.sd"],
    rep(params$numReplicates[[args$selection]],length(js))
  ),region=unique(filteredScores[js,"region"]))
}))

#extract corresponding provean scores
proveanScores <- provean[rownames(jointScores),"SCORE"]
#calculate Correlation
default.cor <- cor(na.omit(cbind(jointScores[,"mj"],proveanScores)))[1,2]


cat(sprintf("PCC with default modes: %.03f",default.cor),"\n")
#draw scatterplot
pdf(paste0(outbase,"_default_cor.pdf"),4,4)
plot(na.omit(cbind(jointScores[,"mj"],proveanScores)),
  xlab="score",ylab="PROVEAN",pch=".",
  main=sprintf("PCC = %.03f",default.cor)
)
invisible(dev.off())


#iterate over search space samples
pcc <- do.call(c,pbmclapply(1:nrow(searchSpace), function(i) {

  #calculate new scores for the current set of modes
  synModes <- unlist(searchSpace[i,])[filteredScores$region]
  nonModes <- unlist(searchSpace[i,])[filteredScores$region+nregions]
  newscores <- (filteredScores$logPhi - nonModes) / (synModes - nonModes)
  newsd <- filteredScores$logPhi.sd / (synModes - nonModes)

  #collapse by AA outcome
  jointScores <- do.call(rbind,tapply(1:nrow(filteredScores),filteredScores$hgvsp,function(js) {
    join.datapoints(
      newscores[js], newsd[js],
      rep(params$numReplicates[[args$selection]],length(js))
    )
  }))
  #lookup corresponding provean scores
  cor(na.omit(cbind(jointScores[,"mj"],proveanScores)))[1,2]

},mc.cores=args$cores))

#write table of all results to file
result <- cbind(searchSpace,pcc=pcc)
write.csv(as.data.frame(result),paste0(outbase,"_PCCs.csv"),row.names=FALSE)

#draw histogram of PCC distribution
pdf(paste0(outbase,"_pccDistribution.pdf"),4,4)
hist(pcc,breaks=100,col="gray",border=NA)
invisible(dev.off())

#' Delegation function to draw score distribution plots with a given filter setting
#' 
#' @param aaScores the score table
#' @param seCutoff the stderr cutoff to apply
#' @return NULL
drawDistributions <- function(aaScores,defaults,best,seCutoff=Inf,reg=NA) {
  #extract filtered scores
  synScores <- with(aaScores,logPhi[grepl("=$",hgvsp) & se < seCutoff ])
  stopScores <- with(aaScores,logPhi[grepl("Ter$",hgvsp) & se < seCutoff])
  misScores <- with(aaScores,logPhi[!grepl("Ter$|=$",hgvsp) & se < seCutoff])
  allScores <- c(synScores,stopScores,misScores)

  #calculate plot ranges to nearest integers
  left <- floor(quantile(allScores,0.01,na.rm=TRUE))
  right <- ceiling(quantile(allScores,0.99,na.rm=TRUE))
  farleft <- floor(min(c(0,allScores),na.rm=TRUE))
  farright <- ceiling(max(c(1,allScores),na.rm=TRUE))
  #set histogram bins based on range (with extra space on the right)
  breaks <- seq(farleft,farright+0.1,0.1)
  #fill bins
  synHist <- hist(synScores,breaks=breaks,plot=FALSE)
  stopHist <- hist(stopScores,breaks=breaks,plot=FALSE)
  misHist <- hist(misScores,breaks=breaks,plot=FALSE)

  #draw top plot for syn/stop
  op <- par(mar=c(2,4,2,1)+.1)
  xs <- barplot(
    rbind(synHist$density,stopHist$density),
    beside=TRUE,col=c("darkolivegreen3","firebrick3"),
    border=NA,ylab="density",space=c(0,0),
    main=if (is.infinite(seCutoff)) {
      paste("Region",reg,"; Unfiltered")
    } else {
      bquote("Region"~.(reg)~";"~sigma < .(seCutoff))
    }
  )
  grid(NA,NULL)
  #add x-axis
  pips <- left:right
  pipx <- colMeans(xs[,which(breaks %in% pips)])
  axis(1,at=pipx,labels=pips)
  #draw legend
  legend("left",
    c("nonsense","synonymous","missense"),
    fill=c("firebrick3","darkolivegreen3","gray"),
    border=NA,bty="n"
  )
  #draw bottom plot (for missense)
  par(mar=c(1,4,0,1)+.1)
  xs <- barplot(
    -misHist$density,
    border=NA,ylab="density",space=0,
  )
  grid(NA,NULL)
  #establish user coordinate function on barplot
  x0 <- xs[which(breaks==0),1]
  x1 <- xs[which(breaks==1),1]
  uCoord <- function(x)  x0 + x*(x1-x0)
  #draw medians
  abline(v=uCoord(defaults[["syn"]]),col="darkolivegreen4",lwd=2,lty="dotted")
  abline(v=uCoord(defaults[["non"]]),col="firebrick3",lwd=2,lty="dotted")
  text(
    uCoord(defaults[["syn"]]),
    -max(misHist$density)/3,
    sprintf("default = %.03f",defaults[["syn"]]),
    col="darkolivegreen4"
  )
  text(
    uCoord(defaults[["non"]]),
    -2*max(misHist$density)/3,
    sprintf("default = %.03f",defaults[["non"]]),
    col="firebrick3"
  )

  abline(v=uCoord(best[["syn"]]),col="darkolivegreen4",lwd=4,lty="dotted")
  abline(v=uCoord(best[["non"]]),col="firebrick3",lwd=4,lty="dotted")
  text(
    uCoord(best[["syn"]]),
    -max(misHist$density)/5,
    sprintf("best = %.03f",best[["syn"]]),
    col="darkolivegreen4"
  )
  text(
    uCoord(best[["non"]]),
    -4*max(misHist$density)/5,
    sprintf("best = %.03f",best[["non"]]),
    col="firebrick3"
  )
  par(op)

}


#Less severely filtered set for visualizations of score distributions
filteredScores2 <- allscores[is.na(allscores$filter),]
#calculate logphi values per AA change
jointLogPhi <- as.data.frame(do.call(rbind,tapply(1:nrow(filteredScores2),filteredScores2$hgvsp,function(js) {
  c(join.datapoints(
    filteredScores2[js,"logPhi"], filteredScores2[js,"logPhi.sd"],
    rep(params$numReplicates[[args$selection]],length(js))
  ),region=unique(filteredScores2[js,"region"]))
})))
colnames(jointLogPhi) <- c("logPhi","se","df","region")
jointLogPhi$hgvsp <- rownames(jointLogPhi)

# maxi <- which.max(pcc)


#plot the best mode combination over the distributinos
pdf(paste0(outbase,"_marginalPCCs.pdf"),4,10)
layout(cbind(1:nregions))
for (regi in 1:nregions) {
  plot(NA,type="n",xlim=c(-1,.3),ylim=c(0,.3),
    xlab=expression("mode parameters"~(log(phi))),
    ylab="marginal mean PCC",
    main=paste("Region",regi)
  )
  syn1means <- tapply(result[,"pcc"],result[,paste0("syn",regi)], mean)
  syn1sd <- tapply(result[,"pcc"],result[,paste0("syn",regi)], sd)
  lines(as.numeric(names(syn1means)),syn1means,col="darkolivegreen3")
  errorBars(as.numeric(names(syn1means)),syn1means,syn1sd,col="darkolivegreen3")

  non1means <- tapply(result[,"pcc"],result[,paste0("non",regi)], mean)
  non1sd <- tapply(result[,"pcc"],result[,paste0("non",regi)], sd)
  lines(as.numeric(names(non1means)),non1means,col="firebrick3")
  errorBars(as.numeric(names(non1means)),non1means,non1sd,col="firebrick3")
}
invisible(dev.off())

#extract the parameters with the best marginals
bestModes <- do.call(rbind,lapply(1:nregions, function(regi) {
  syn1means <- tapply(result[,"pcc"],result[,paste0("syn",regi)], mean)
  non1means <- tapply(result[,"pcc"],result[,paste0("non",regi)], mean)
  c(syn=as.numeric(names(which.max(syn1means))), non=as.numeric(names(which.max(non1means))))
}))


pdf(paste0(outbase,"_logPhiDistr_bestModes.pdf"),7,14)
layout(cbind(1:(2*nregions)))
for (regi in 1:nregions) {
  drawDistributions(
    jointLogPhi[jointLogPhi$region==regi,],
    defaults=c(syn=as.vector(defaultModes[regi,"synMed"]),non=as.vector(defaultModes[regi,"nonMed"])),
    # best=c(syn=searchSpace[maxi,paste0("syn",regi)],non=searchSpace[maxi,paste0("non",regi)]),
    best=bestModes[regi,],
    reg=regi
  )
}
invisible(dev.off())



#calculate new scores for the best set of modes
synModes <- bestModes[,"syn"][filteredScores$region]
nonModes <- bestModes[,"non"][filteredScores$region]
newscores <- (filteredScores$logPhi - nonModes) / (synModes - nonModes)
newsd <- filteredScores$logPhi.sd / (synModes - nonModes)

#collapse by AA outcome
jointScores <- do.call(rbind,tapply(1:nrow(filteredScores),filteredScores$hgvsp,function(js) {
  join.datapoints(
    newscores[js], newsd[js],
    rep(params$numReplicates[[args$selection]],length(js))
  )
}))
#lookup corresponding provean scores
best.cor <- cor(na.omit(cbind(jointScores[,"mj"],proveanScores)))[1,2]


cat(sprintf("PCC with optimal modes: %.03f",best.cor),"\n")
#draw scatterplot
pdf(paste0(outbase,"_best_cor.pdf"),4,4)
plot(na.omit(cbind(jointScores[,"mj"],proveanScores)),
  xlab="score",ylab="PROVEAN",pch=".",
  main=sprintf("PCC = %.03f",best.cor)
)
invisible(dev.off())

colnames(jointScores) <- c("score","se","df")
bestScores <- data.frame(hgvs_pro=rownames(jointScores),jointScores)
write.csv(bestScores,paste0(outbase,"_optimized_map.csv"),row.names=FALSE)

cat("Done!\n")
