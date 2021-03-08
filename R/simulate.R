# Copyright (C) 2018  Jochen Weile, Roth Lab
#
# This file is part of tileseqMave.
#
# tileseqMave is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# tileseqMave is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with tileseqMave.  If not, see <https://www.gnu.org/licenses/>.


#' Simulate Pool
#'
#' This simulates a pool of variant-bearing clones. The output is a list containing
#' two elements: (i) a table of all underlying codon changes and their true marginal
#' frequencies, and (ii) a function that samples (multi-variant-bearing) clone populations
#' on demand
#'
#' @param params a parameter sheet to follow (in terms of CDS sequence etc.)
#' @param poolCV the coefficient of variation for the distribution of marginal frequencies
#' @param poolLambda the mean of the poisson distibution used to sample the number of 
#'   codon-changes per clone
#'
#' @return
#' @export
#'
#' @examples
simulatePool <- function(params,poolCV=0.2,poolLambda=1) {
  library("hash")
  #generate all possible codons
  allCodons <- apply(do.call(expand.grid,
    replicate(3,yogitools::toChars("ACGT"),simplify=FALSE)
  ),1,paste,collapse="")
  #import translation table
  data(trtable)
  #generate all possible codon changes
  allCodonChanges <- do.call(rbind,lapply(2:params$template$proteinLength,function(pos) {
    from <- params$template$cdsCodons[[pos]]
    to <- setdiff(allCodons,from)
    data.frame(id=paste0(from,pos,to),from=from,pos=pos,to=to,
      fromAA=trtable[[from]],toAA=sapply(to,function(toc)trtable[[toc]]),
      tile=params$pos2tile(pos),region=params$pos2reg(pos)
    )
  }))
  #sample the underlying frequency of each variant
  freqs <- rnorm(nrow(allCodonChanges),1,poolCV)
  allCodonChanges$freq <- freqs-min(freqs)
  
  #sampling function
  sampleVariants <- function(n,region) {
    nmut <- rpois(n,poolLambda)
    with(allCodonChanges[allCodonChanges$region==region,],{
      prob <- freq/sum(freq)
      table(sapply(nmut,function(m) {
        if (m == 0) "WT" else paste(sample(id,m,prob=prob),collapse="|")
      }))
      # Thought a hash-based approach would be faster, but... no
      # counts <- hash::hash()
      # for (m in nmut) {
      #   mutid <- if (m == 0) "WT" else paste(sample(id,m,prob=prob),collapse="|")
      #   if (hash::has.key(mutid,counts)){
      #     counts[[mutid]] <- counts[[mutid]]+1
      #   } else {
      #     counts[[mutid]] <- 1
      #   }
      # }
      # hash::values(counts)
    })
  }
  
  list(allCodonChanges=allCodonChanges,sampleVariants=sampleVariants)
}


#' Simulate Clone fitness
#' 
#' This function generates assay fitness values for clones and their variants.
#' It outputs a list with two elements: (i) a table listing the fitness effects of
#' All possible nucleotide changes, and (ii) a function that calculates the true fitness
#' of any given clone given its variant combination. Currently, this does not simulate
#' genetic interactions, but simply uses multiplicativity of fitness. It also currently 
#' does not simulate hypercomplementation effects.
#'
#' @param params The parameter sheet
#' @param pool the output of simulatePool()
#' @param lofFreq the frequency of loss-of-function effects to simulate
#' @param hypoFreq the frequency of hypomorphic effects to simulate
#'
#' @return
#' @export
#'
#' @examples
simulateFitness <- function(params,pool,lofFreq=0.3,hypoFreq=0.1,hypoSD=0.1) {
 
  allAAChanges <- expand.grid(pos=2:params$template$proteinLength,toAA=sort(unique(pool$allCodonChanges$toAA)))
  allAAChanges$fromAA <- sapply(allAAChanges$pos, function(p) substr(params$template$proteinSeq,p,p))
  allAAChanges$id <- with(allAAChanges,paste0(pos,toAA))
  allAAChanges <- allAAChanges[,c("id","fromAA","pos","toAA")]
  rownames(allAAChanges) <- allAAChanges$id
  
  allAAChanges$trueEffect <- sapply(1:nrow(allAAChanges), function(i) with(allAAChanges[i,],{
    if (fromAA==toAA) 1 
    else if (toAA=="*") 0
    else {
      rand <- runif(1,0,1)
      if (rand < lofFreq) 0
      else if (rand < lofFreq+hypoFreq) rnorm(1,0.5,hypoSD)
      else 1
    }
  }))
  
  calcTrueFitness <- function(codonChanges) {
    #import translation table
    data(trtable)
    sapply(strsplit(codonChanges,"\\|"), function(ccs) {
      if (length(ccs)==1 && ccs=="WT") return (1)
      # from <- substr(ccs,1,3)
      pos <- as.integer(substr(ccs,4,nchar(ccs)-3))
      to <- substr(ccs,nchar(ccs)-2,nchar(ccs))
      # fromAA <- sapply(from,function(fromc)trtable[[fromc]])
      toAA <- sapply(to,function(toc)trtable[[toc]])
      ids <- paste0(pos,toAA)
      fits <- allAAChanges[ids,"trueEffect"]
      prod(fits)
    })
  }
  
  list(allAAChanges=allAAChanges,calcTrueFitness=calcTrueFitness)
  
}


#' Simulate selection assay
#' 
#' This simulates the effects of a growth-based fitness assay. It samples a pool of cells
#' of the desired size and applies fitness effects with a pre-determined amount of noise. 
#' It then outputs the clone frequencies in the pre- and post-selection pool
#'
#' @param params the parameter sheet
#' @param pool the output of simulatePool()
#' @param trueFit the output of simulateFitness()
#' @param n The size of the pre-selection pool (number of cells)
#' @param sd the standard deviation of the noise affecting the growth assay.
#' @param t the selection time in generations
#'
#' @return a list of vectors containing the clone frequencies in the pre- and post-selection pool
#' @export
#'
#' @examples
simulateAssay <- function(params,plasmidPools,trueFit,n,sd=0.01,t=3) {
  
  #sample pre-selection clones
  logInfo("Sampling pre-selection pool")
  # presel <- do.call(c,lapply(params$regions$`Region Number`, function(reg) {
  #   pool$sampleVariants(n,reg)
  # }))
  presel <- do.call(c,lapply(plasmidPools, function(ppool) {
    setNames(rpois(length(ppool),n*ppool/sum(ppool)),names(ppool))
  }))
  if (any(presel == 0)) {
    presel <- presel[presel > 0]
  }
  
  logInfo("Applying selection")
  #create random noise from normal distribution
  noise <- rnorm(length(presel),0,sd)
  #calculate true fitness values
  f <- trueFit$calcTrueFitness(names(presel))
  #calculate post-selection clones given fitness and noise
  postsel <- presel * t * 2^(f+noise/sqrt(presel))
  
  #return the pre- and post-selection clones
  list(presel=presel,postsel=postsel)
  
}

#' simulate TileSeq
#'
#' @param assayReps a list of clone pool replicates. These can be pre- or post-selection
#'   or even all WT
#' @param params the parameter sheet
#' @param depth the sequencing depth to simulate
#' @param varcallErr the variant calling error to simulate
#'
#' @return
#' @export
#'
#' @examples
simulateTileSeq <- function(assayReps,params,depth=1e5,varcallErr=3e-6) {
  
  #create list of possible SNVs per each tile
  snvList <- reachableChanges(params)
  snvList$tile <- params$pos2tile(snvList$pos)
  snvsPerTile <- tapply(1:nrow(snvList),snvList$tile,function(rows) with(snvList[rows,],{
    do.call(c,mapply(
      paste0,
      from=wtcodon,
      pos=pos,
      to=strsplit(mutcodons,"\\|")
    ))
  }))
  
  #number of nucleotides in each tile
  tileNClen <- apply(params$tiles,1,function(ti)ti[[5]]-ti[[4]]+1)

  #iterate over replicate/condition combos
  allReads <- do.call(c,lapply(assayReps,function(ar) {
    #the poisson lambda at the given depth
    l <- depth*ar/sum(ar)
    #separated mutations in each molecule
    muts <- strsplit(names(ar),"\\|")
    #iterate over tiles
    readsPerTile <- lapply(1:nrow(params$tiles), function(ti) {
      #get the actual tile ID
      tileID <- params$tiles[ti,"Tile Number"]
      #determine which mutations are visible in the tile
      tilemuts <- sapply(muts,function(ms) {
        #wt remains wt
        if (length(ms)==1 && ms == "WT") {
          return("WT")
        } else {
          #extract those muations that are visible in the current tile
          pos <- as.integer(substr(ms,4,nchar(ms)-3))
          tileMuts <- ms[params$pos2tile(pos)==tileID]
          #if none are visible, it will look like wt
          if (length(tileMuts)==0) {
            return("WT")
          } else {
            return(paste(tileMuts,collapse="|"))
          }
        }
      })
      #sample reads from poisson distribution over available molecules
      reads <- rpois(length(ar),l)
      
      #simulate errors
      nSeqErr <- rpois(1,tileNClen[[ti]]*varcallErr*depth)
      #sample random errors
      seqErrs <- sample(snvsPerTile[[ti]],nSeqErr,replace=TRUE)
      #pick some reads that will "host" the errors
      hostReads <- sample(tilemuts,nSeqErr,replace=TRUE,prob=reads/sum(reads))
      hostWithErr <- mapply(
        function(hr,er){
          if (hr=="WT") er else paste(hr,er,sep="|")
        },
        hostReads,seqErrs
      )
      
      #remove the host reads and replace with error-laden versions
      #also combine reads that look identical due to limited tile range
      apparentReads <- tapply(
        c(reads,rep(-1,nSeqErr),rep(1,nSeqErr)),
        c(tilemuts,hostReads,hostWithErr),
        sum
      )
      
      return(apparentReads[apparentReads>0])
      
    })
    setNames(readsPerTile,paste0("tile",params$tiles[,"Tile Number"]))
  }))
  
  return(allReads)
}

#' Export count data to file
#'
#' @param reads the read information
#' @param name the sample name (matching the sample sheet)
#' @param outdir the output directory
#' @param params the parameter sheet
#' @param hgvsb a (coding sequence) HGVS builder object
#'
#' @return
#' @export
#'
#' @examples
exportSeqSample <- function(reads, name, outdir, params, hgvsb=new.hgvs.builder.c()) {
  
  wtCount <- reads[["WT"]]
  
  aacs <- strsplit(names(reads),"\\|")
  hgvscs <- sapply(aacs,function(cisaacs) {
    
    hgvss <- sapply(cisaacs,function(aac) {
      
      if (aac=="WT") {
        return("WT")
      }
      n <- nchar(aac)
      aapos <- as.integer(substr(aac,4,n-3))
      ncPos <- (aapos-1)*3 + 1:3
      ncAnc <- substr(rep(aac,3),1:3,1:3)
      ncVar <- substr(rep(aac,3),n-3+1:3,n-3+1:3)
      chIdx <- which(sapply(1:3,
        function(i) substr(aac,i,i) != substr(aac,n-3+i,n-3+i)
      ))
      if (length(chIdx)==1) {
        #then it's an SNV
        hgvsb$substitution(ncPos[[chIdx]],ncAnc[[chIdx]],ncVar[[chIdx]])
      } else {
        #it's an MNV
        hgvsb$delins(
          ncPos[min(chIdx)],
          ncPos[max(chIdx)],
          paste(ncVar[min(chIdx):max(chIdx)],collapse="")
        )
      }
    })
    if (length(hgvss)==1) {
      return(hgvss)
    } else {
      return(do.call(hgvsb$cis,as.list(hgvss)))
    }
  })
  
  outData <- data.frame(HGVS=hgvscs,count=reads)
  #remove WT, as this information goes in the header instead
  outData <- outData[-which(hgvscs=="WT"),]
  
  sInfo <- params$samples[params$samples[,1]==name,]
  tile <- sInfo$`Tile ID`
  tileInfo <- params$tiles[params$tiles[,1]==tile,]
  header <- paste0(
    "#Sample:",name,"\n",
    "#Tile:",tile,"\n",
    "#Tile Starts:",tileInfo[["Start AA"]],"\n",
    "#Tile Ends:",tileInfo[["End AA"]],"\n",
    "#Condition:",sInfo$Condition,"\n",
    "#Replicate:",sInfo$Replicate,"\n",
    "#Timepoint:",sInfo$`Time point`,"\n",
    "#Final read-depth:",sum(reads),"\n"
  )
  filename <- paste0(outdir,"counts_sample_",name,".csv")
  con <- file(filename, open="w")
  cat(header,file=con)
  write.csv(outData,con,row.names=FALSE)
  close(con)
  
}

#' Export simulated experiment
#'
#' @param reads list of sample read counts with names corresponding to samples
#' @param params the parameter object
#' @param outdir the target output folder to export to
#'
#' @return nothing
#' @export
#'
#' @examples
exportExperiment <- function(reads,params,outdir) {
  sampleParams <- yogitools::extract.groups(names(reads),"^rep(\\d+)\\.(.+)\\.tile(\\d+)$")
  sampleParams <- data.frame(
    rep=as.integer(sampleParams[,1]),
    cond=sampleParams[,2],
    tile=as.integer(sampleParams[,3])
  )
  sampleNames <- lapply(1:nrow(sampleParams),function(i) with(sampleParams[i,],{
    with(params$samples,{
      `Sample ID`[which(Condition==cond & Replicate==rep & `Tile ID`==tile)]
    })
  }))
  hgvsb <- new.hgvs.builder.c()
  
  mapply(function(snames,sreads){
    lapply(snames,function(sname) exportSeqSample(sreads, sname, outdir, params, hgvsb))
  },sampleNames,reads)
  
  return(invisible(NULL))
}


simulateExperiment <- function(workdir=tempdir(),paramFile,
                               poolCV=0.2,poolLambda=1,
                               lofFreq=0.3,hypoFreq=0.1,nPlasmids=1e5,
                               nClones=1e5,assaySD=0.01,assayTime=3,
                               seqDepth=5e6,varcallErr=3e-6) {
  
  params <- parseParameters(paramFile)
  
  if (!grepl("/$",workdir)) {
    workdir <- paste0(workdir,"/")
  }
  
  timestamp <- format(Sys.time(),"%Y-%m-%d-%H-%M-%S")
  
  logInfo("Generating variant pool")
  pool <- simulatePool(params,poolCV,poolLambda)
  
  logInfo("Generating mock fitness landscape")
  trueFit <- simulateFitness(params,pool,lofFreq,hypoFreq)
  
  
  plasmidPools <- lapply(params$regions$`Region Number`, function(ri) {
    pool$sampleVariants(n = nPlasmids, region = ri)
  })
  
  assayReps <- do.call(c,lapply(getSelects(params),function(sCond) {
    
    quads <- findQuads(params,sCond,"1")
    repNames <- paste0("rep",1:ncol(quads$repMatrix))
    
    logInfo("Simulating assay for",sCond)
    condReps <- do.call(c,setNames(lapply(repNames, function(repi) {
      setNames(
        simulateAssay(params,plasmidPools,trueFit,nClones,sd=assaySD,t=assayTime),
        c(quads$condQuad[["nonselect"]],quads$condQuad[["select"]])
      )
    }),repNames))
    wtReps <- setNames(
      replicate(length(repNames),c(WT=nClones),simplify=FALSE),
      paste0(repNames,".",quads$condQuad[["nonWT"]])
    )
    
    return(c(condReps,wtReps))
    
  }))
  
  logInfo("Simulating TileSeq")
  sampleReads <- simulateTileSeq(assayReps,params,depth=seqDepth,varcallErr=varcallErr) 
  
  logInfo("Exporting simulated data")
  outdir <- paste0(workdir,"simul_",timestamp,"_mut_count/")
  dir.create(outdir,recursive = TRUE, showWarnings = FALSE)
  exportExperiment(sampleReads,params,outdir)
  
  buildJointTable(dataDir = workdir, inDir = outdir, paramFile=paramFile)
  
  # libraryQC(dataDir=workdir,inDir=outdir, paramFile=paramFile)
  
  calcEnrichment(workdir, outdir, paramFile = paramFile)
  
  # selectionQC(workdir, outdir, paramFile = paramFile)
  
  margfile <- paste0(outdir,"marginalCounts.csv")
  # enrfile <- paste0(workdir,"simul_",timestamp,"_scores/",sCond,"_t1_enrichment.csv")
  margs <- read.csv(margfile, comment.char="#")

  rownames(margs) <- margs$codonChange
  # rownames(pool$allCodonChanges) <- pool$allCodonChanges$id
  margComp <- cbind(
    true=pool$allCodonChanges$freq/sum(pool$allCodonChanges$freq), 
    observed=margs[pool$allCodonChanges$id,"nonselect.t1.rep1.frequency"]
  )
  plot(margComp)
  
  return(list(countDir=outdir,pool=pool,trueFitness=trueFit))
  
}


