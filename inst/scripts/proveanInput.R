options(stringsAsFactors=FALSE)

library(yogitools)
library(tileseqMave)
library(argparser)

#process command line arguments
p <- arg_parser(
	"Generates an input file for the PROVEAN webtool based on the ORF sequence.",
	name="mavevisLocal.R"
)
p <- add_argument(p, "paramfile", help="parameter file")
p <- add_argument(p, "--output", help="output file")
args <- parse_args(p)

#read parameter sheet
params <- parseParameters(args$paramfile)

#extract wt sequence from parameter file
wt.aa <- toChars(params$template$proteinSeq)
#make list of all possible AAs
aas <- toChars("AVLIMFYWRHKDESTNQGCP")
#create a table of all position / AA combinations
featable <- expand.grid(pos=1:length(wt.aa),mut.aa=aas)
#use that table to create the provean input
provean.in <- data.frame(
	protein=params$template$uniprot,
	pos=featable$pos,
	wt=wt.aa[featable$pos],mut=featable$mut.aa
)
#write results to file
outfile <- if (is.na(args$output)) paste0(params$template$uniprot,"_provean_input.txt") else args$output
write.table(provean.in,outfile,
	sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE
)

#Done
cat("Results written to ",outfile,"\n")
