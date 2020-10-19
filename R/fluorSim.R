options(stringsAsFactors=FALSE)
library(sn)

target <- rbind(
  mean=log10(c(WT=123867.23,del1=81996.05,del2=51562.45,pool=114515.45
  )),
  mode=c(WT=5.12,del1=4.96,del2=4.72,pool=5.16
  ),
  sd=c(
    WT=54544.7/(log(10)*123867.23),
    del1=36065.51/(log(10)*81996.05),
    del2=26505.39/(log(10)*51562.45),
    pool=62350.01/(log(10)*114515.45)
  )
)

params <- rbind(
  xi=c(WT=5.31, del1=5.12, del2=4.92, pool=5.3),
  omega=c(WT=0.28,del1=0.28,del2=0.32,pool=0.35),
  alpha=c(WT=-2,del1=-2,del2=-2,pool=-2)
)
# 
# plotcolors=c(WT="black",del1="firebrick3",del2="firebrick2",pool="blue")
# stats <- sapply(colnames(params),function(currTarget){
#   param <- c(xi=5.3,omega=0.3,alpha=-2)
#   samples <- sn::rsn(n=100000,dp=params[,currTarget])
#   m <- mean(samples)
#   s <- sd(samples)
#   histo <- hist(samples,breaks=1000,plot=FALSE)
#   mode <- histo$mids[[which.max(histo$counts)]]
#   curve(
#     dsn(x,dp=params[,currTarget]),
#     from=2,to=6,
#     add=currTarget!="WT",
#     col=plotcolors[[currTarget]],
#     xlab=expression(log[10](fluorescence)),
#     ylab="density"
#   )
#   c(mean=m,sd=s,mode=mode)
# })
# legend("left",c("WT","del1","del2","pool"),lty=1,col=plotcolors)
# target[c(1,3,2),]
# stats

#helper function to floor numbers at zero and turn them to integers
posInt <- function(xs) sapply(xs, function(x) if (x < 0) 0 else round(x))

#helper function to plot bin log-phi trajectories
plotTrails <- function(rows) {
  plotColors=c(WT="black",del1="purple",del2="firebrick2")
  plot(NA,type="n",
       xlim=c(1,length(bin.names)),ylim=range(logPhis),
       axes=FALSE,ylab=expression(log[10](phi)),xlab="bins",
       main=sprintf("%.01e reads",total.reads)
  )
  axis(2)
  axis(1,at=1:length(bin.names),bin.names)
  for (row in rows) {
    lines(1:length(bin.names),logPhis[row,],col=plotColors[[true.class[[row]]]])
  }
}

bin.options <- list(
  c(0.1,0.1,0.6,0.1,0.1),
  rep(1/6,6),
  rep(1/5,5),
  rep(1/4,4),
  rep(1/3,3)
)
read.options <- c(1e6)
all.options <- expand.grid(bin.options,read.options)

for (it in 1:nrow(all.options)) {
  
  bin.widths <- all.options[it,1][[1]]
  total.reads <- all.options[it,2]

  #simulated true fitness levels
  prot.length <- 860
  N <- prot.length*3*32
  true.class <- sample(c("WT","del1","del2"),N,replace=TRUE,prob=c(0.3,0.1,0.6))
  nonselect.logfreq <- rnorm(N,mean=-3.7,sd=0.4)
  nonselect.marCount <- posInt((10^nonselect.logfreq)*total.reads)
  
  all.fluor.points <- do.call(c,lapply(1:N, function(i) {
    xi <- params["xi",true.class[[i]]]
    omega <- params["omega",true.class[[i]]]
    alpha <- params["alpha",true.class[[i]]]
    sn::rsn(n=nonselect.marCount[[i]],xi=xi,omega=omega,alpha=alpha)
  }))
  
  # hist(all.fluor.points,breaks=500,
  #   border=NA,col="gray",freq=FALSE,
  #   xlab=expression(log[10](fluorescence)),
  #   main="Fluorescence of pool"
  # )
  
  bin.borders <- c(0,cumsum(bin.widths))
  bin.names <- sapply(1:length(bin.widths),function(i) sprintf("[%.01f;%.01f]",bin.borders[[i]],bin.borders[[i+1]]))
  bin.qs <- quantile(all.fluor.points,bin.borders)
  b <- length(bin.names)
  
  probs <- yogitools::as.df(lapply(1:N, function(i) {
    sapply(1:b,function(bin) {
      left <- sn::psn(x=bin.qs[[bin]],dp=params[,true.class[[i]]])
      right <- sn::psn(x=bin.qs[[bin+1]],dp=params[,true.class[[i]]])
      prob <- right-left
      prob
    })
  }))
  colnames(probs) <- bin.names
  
  readCounts <- yogitools::as.df(lapply(1:N, function(i) {
    n <- nonselect.marCount[[i]]*b/(b+1)
    lambda <- n * unlist(probs[i,])
    rpois(ncol(probs)+1,c(lambda,n/(b+1)))
  }))
  colnames(readCounts) <- c(bin.names,"all")
  
  logPhis <- t(apply(readCounts, 1, function(row) {
    row <- row+0.1
    log10(row/row[["all"]])[-length(row)]
  }))
  
  linfits <- t(apply(logPhis,1,function(row) {
    x <- 1:b
    coefficients(lm(row ~ x))
  }))
  
  pdf(sprintf("%.00ereads_%s.pdf",total.reads,paste(round(bin.widths*100),collapse="-")),10,5)
  op <- par(mfrow=c(1,2))
  plotColors=c(WT="black",del1="purple",del2="firebrick2")
  plot(linfits,pch=".",
       col=plotColors[true.class],
       xlab="intercept",ylab="slope"
  )
  plotTrails(1:100)
  par(op)
  dev.off()
  
}
