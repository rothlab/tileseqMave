#!/usr/bin/env Rscript

#parse CLI arguments
library(argparser)
p <- arg_parser(
	"Create links to Rscripts in target directory",
	name="linkBinaries.R"
)
p <- add_argument(p, "targetDir", help="target directory", default="~/.local/bin/")
args <- parse_args(p)

#list of scripts to link
scripts <- c(
	"csv2json.R","joinCounts.R","runLibraryQC.R",
	"calcEnrichment.R","runSelectionQC.R","scaleScores.R","mavevisLocal.R",
	"condenseQC.R","popcodeOligos.R","colorizeStructure.R","differenceMap.R",
	"proveanInput.R"
)

#find scripts folder in the local library installation
scriptsFolder <- system.file("scripts/",
	package = "tileseqMave",
	mustWork = TRUE
)

#function to create symlinks
linkScript <- function(scriptName,targetDir="~/.local/bin/") {
	scriptPath <- paste0(scriptsFolder,scriptName)
	if (!file.exists(scriptPath)) {
		stop(scriptName, "not found!")
	}
	file.symlink(from=scriptPath,to=paste0(targetDir,scriptName))
}

#create target dir if it doesn't exist
if (!dir.exists(args$targetDir)) {
	dir.create(args$targetDir,recursive=TRUE)
}

#run the function
lapply(scripts,linkScript,targetDir=args$targetDir)
