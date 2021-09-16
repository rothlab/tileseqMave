#!/usr/bin/env Rscript

options(
  stringsAsFactors=FALSE,
  ignore.interactive=TRUE
)

library(TmCalculator)
library(tileseqMave)
library(argparser)
library(yogitools)
library(pbmcapply)
library(hash)

data(trtable)

#process command line arguments
p <- arg_parser(
  "Popcode oligo designer",
  name="popcodeOligos.R"
)
p <- add_argument(p, "fasta", help="input FASTA file. Must contain prefix, ORF, and suffix entries")
p <- add_argument(p, "--length", help="Oligo length to aim for",default=33L)
p <- add_argument(p, "--wiggle", help="Maximum wiggle-room for length",default=5L)
p <- add_argument(p, "--out", help="output filename base.",default="popcodeOligos")
p <- add_argument(p, "--method", help="melting temperature calculation method. NN, GC or Wallace",default="NN")
args <- parse_args(p)

stopifnot(!is.null(args$fasta), !is.na(args$fasta))
oligo.length <- args$length
stopifnot(oligo.length > 10)
wiggle <- args$wiggle
stopifnot(wiggle >= 0, wiggle < oligo.length/2)
calcTm <- switch(match.arg(args$method,c("NN","GC","Wallace")), NN=Tm_NN, GC=Tm_GC, Wallace=Tm_Wallace)
outfile <- args$out
stopifnot(!is.null(outfile), !is.na(outfile))

#helper function to read FASTA files
readMultiFASTA <- function(infile) {
  lines <- readLines(infile)
  headerLines <- grep("^>",lines)
  ids <- sub("^>","",lines[headerLines])
  bodyStarts <- headerLines+1
  bodyEnds <- c(headerLines[-1]-1,length(lines))
  setNames(lapply(1:length(ids),function(i) paste(lines[bodyStarts[[i]]:bodyEnds[[i]]],collapse="")),ids)
}

#helper function to write FASTA files
writeMultiFasta <- function(oligos,filename) {
  lines <- do.call(c,lapply(1:length(oligos),function(i) {
    c(
      paste0(">",names(oligos)[[i]]),
      oligos[[i]]
    )
  }))
  writeLines(lines,filename)
}

seqs <- readMultiFASTA(args$fasta)
seqs$construct <- do.call(paste0,seqs[c("prefix","ORF","suffix")])

cat("Exploring possible oligos...\n")
codon.starts <- nchar(seqs$prefix) + seq(1+3,nchar(seqs$ORF),3)
n.choices <- length(codon.starts) * (2*wiggle+1)^2
oligo.choices <- as.df(do.call(c,pbmclapply(codon.starts, function(codon.start) {
  do.call(c,lapply(-wiggle:wiggle, function(left.offset) {
    lapply(-wiggle:wiggle, function(right.offset) {
      start <- codon.start + 1 - floor(oligo.length/2) + left.offset
      end <- codon.start + 1 + floor(oligo.length/2) + right.offset
      sequence <- substr(seqs$construct,start,end)
      list(
        codon.start.in.construct=codon.start,
        codon.start.in.oligo=codon.start-start+1,
        oligo.start=start,
        oligo.end=end,
        sequence=sequence,
        tm.tot=calcTm(sequence),
        tm.left=calcTm(substr(sequence,1,codon.start-start)),
        tm.right=calcTm(substr(sequence,codon.start-start+4,nchar(sequence)))
      )
    })
  }))
},mc.cores=8)))

cat("Optimizing melting temperatures...\n")
median.tms <- apply(oligo.choices[,c("tm.left","tm.right")],2, median)
deviation <- rowSums(abs(oligo.choices[,c("tm.left","tm.right")]-t(replicate(n.choices,median.tms))))
min.idx <- tapply(deviation,oligo.choices$codon.start.in.construct, which.min)
best.oligos <- oligo.choices[(1:length(codon.starts)-1) * (2*wiggle+1)^2 + min.idx,]



cat("Plotting offset distribution...\n")

png(paste0(outfile,".png"),width=600,height=300)
op <- par(mfrow=c(1,2),cex=.9)
  
offset.counts <- table(best.oligos$codon.start.in.oligo-1-floor((oligo.length-3)/2))

barplot(
  offset.counts,
  xlab="Offset",
  ylab="Frequency",
  col="gold2",
  border="gold4"
)
grid(nx=NA,ny=NULL)

hist(
  best.oligos$tm.tot,
  breaks=20,
  col="steelblue3",
  border="steelblue4",
  xlab="Melting temperature (C)",
  main=""
)
grid(nx=NA,ny=NULL)
tmm <- median(best.oligos$tm.tot)
abline(v=tmm,lty="dashed",col="darkgray")
text(tmm, (par("usr")[4]-par("usr")[3])/2, paste(format(tmm,digits=4),"C"))

par(op)
invisible(dev.off())


cat("inserting NNKs...\n")
oligos <- do.call(c, lapply(1:nrow(best.oligos), function(i) {
  with(best.oligos[i,], {
    codon <- substr(sequence, codon.start.in.oligo, codon.start.in.oligo+2)
    # aa <- as.character(translate(DNAString(codon)))
    aa <- trtable[[codon]]
    outseq <- sequence
    substr(outseq, codon.start.in.oligo, codon.start.in.oligo+2) <- "NNK"
    outseqs <- outseq
    names(outseqs) <- paste(aa,i,"X",sep="")
    # cat(i,"\t",substr(sequence,1,1),"\t",substr(sequence,nchar(sequence),nchar(sequence)),"\n")
    outseqs
  })
}))

outTable <- data.frame(id=names(oligos),sequence=oligos,
  tm.tot=best.oligos$tm.tot,tm.left=best.oligos$tm.left,tm.right=best.oligos$tm.right
)

cat("Writing output...\n\n")

#write spreadsheet table
write.csv(outTable,paste0(outfile,".csv"),row.names=FALSE)
#write fasta file
writeMultiFasta(oligos,paste0(outfile,".fasta"))

cat("Done!\n")
