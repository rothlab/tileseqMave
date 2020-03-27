#!/usr/bin/env Rscript

#target bin directory
targetDir <- "~/.local/bin/"

#list of scripts to link
scripts <- c(
	"csv2json.R","joinCounts.R","runLibraryQC.R",
	"runScoring.R","runSelectionQC.R","mavevisLocal.R"
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

#run the function
lapply(scripts,linkScript,targetDir=targetDir)
