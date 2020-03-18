# tileseqMave

Analysis functions for TileSEQ.

## Requirements and Installation

tileseqMave requires R &ge; 3.4.4 and has a number of R package dependencies, which will be installed automatically. It also uses the output of [tilseq_mutcount](https://github.com/RyogaLi/tilseq_mutcount), which should be installed separately.

To install tileseqMave, use R devtools:

```R
install.packages("devtools")
devtools::install_github("jweile/tileseqMave")
```

The tileseqMave R package comes with a set of wrapper scripts, which can be symlinked to your preferred `bin/` directory, such that they become available in your `$PATH`. Here's a quick set of R commands to automatically create these symlinks:

```R
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

```

## Overview of the pipeline

The TileSEQ pipeline breaks down into multiple modules:

1. csv2json: This component reads and validates the parameter sheet which contains all relevant options and parameters of the experiment to be analyzed.
2. fastq2count: This component processes the FASTQ files from a sequencing run and produces variant count files. It is implemented in the [tilseq_mutcount](https://github.com/RyogaLi/tilseq_mutcount) python module.
3. joinCounts: Combines the individual count files from each sequencing sample into a single table, tabulating the read counts and relative frequencies for each unique variant in each condition and replicate. It also calculates the protein-level consequences of each variant and produces marginal counts and frequencies.
4. libraryQC: Produces a number of diagnostic graphs for quality control purposes at the variant library level.
5. scoring: Calculates enrichment ratios and fitness scores for each variant in each selection condition. Also performs error regularization and filtering.
6. selectionQC: Produces additional diagnostic graphs at the selection assay level.

## The parameter sheet

tileseqMave uses a central input document called the parameter sheet, which attempts to capture all relevant aspects of the TileSEQ experiment from which the sequencing data was generated. The parameter sheet is provided as a comma-separated-values (CSV) file, following a strict template. The easiest way to create such a file is to use a spreadsheet editor, such as Google Sheets, LibreOffice or Excel. An example spreadsheet can be found [here](https://docs.google.com/spreadsheets/d/1tIblmIFgOApPNzWN2KUwj8BKzBiJ1pOL7R4AOUGrqvE/edit?usp=sharing).

The parameter sheet has the following sections:

1. **Project name**: This can be any name you like
2. **Template construct**: Here we define the sequencing template, its sequence and what it represents.
    * Gene name: The official HGNC gene name
    * Sequence: The full nucleotide sequence of the template including the coding sequence (CDS) and its flanking priming sequences.
    * CDS start: The position in the above sequence at which the coding sequence (CDS) begins. This must be an ATG codon! (Numbering starts at 1, not at 0; don't @ me, nerds. :P )
    * CDS end: The equivalent sequence position at which the CDS ends. This must be in-frame with the start position, i.e. end-start+1 must be divisible by 3.
    * Uniprot Accession: The UniprotKB accession of the protein encoded by the template gene.
3. **Assay**: A summary of the underlying selection assay
    * Assay Type: This can be any free-text label you like. Pick a simple name, like "Y2H", "Yeast complementation" or "LDL uptake via FACS"
    * Selection: This field indicates whether the assay performs a positive or negative selection. So only the values "Positive" and "Negative" are allowed. Positive selection indicates that the assay causes damaging variants to be depleted in the pool, where as negative selection indicates that the assay causes neutral variants to be depleted. (Most assays will be positive).
4. **Conditions and replicates**: This is a custom table that lists the different experimental conditions and their intended number of replicates and time points.
    * List of conditions: In this table row, provide a list of condition identifiers. These should be short and must not contain any special characters.
    * Number of replicates: For each of the defined conditions in the previous row, provide the number of technical replicates here.
    * Number of time points: For each of the defined conditions in the first row, provide the number of timepoints. (At the time of writing, this feature has not been implemented yet).
5. **Mutagenesis regions**: This is a fixed table defining the regions which underwent separate PopCode mutagenesis. It has the following columns:
    * Region number: A simple numbering of the regions to serve as identifiers thereof.
    * Start AA: The first amino acid position that counts as within the region.
    * End AA: The last amino acid position that counts as within the region.
6. **Sequencing tiles**: This also a fixed table, defining the eponymous TileSeq sequencing tiles in the experiment. It has the following columns:
    * Tile number: A simple numbering of the tiles to serve as identifiers thereof.
    * Start AA: The first amino acid position that counts as within the tile.
    * End AA: The last amino acid position that counts as within the tile.
7. **Condition definitions**: This section is used to define the meanings of the different conditions as to how they relate to each other. Currently, two kinds of relationships are supported: `is_selection_for`, and `is_wt_control_for`. For example, if we have defined three conditions in section 4; sel, non and ctrl, then sel could be the selection condition for non and ctrl could be the wt control for either of them. These relationships are defined using a fixed table with three columns:
    * Condition 1: The ID of the first condition in the relationship. Must have been previously declared in the list of conditions.
    * Relationship: The name of the relationship, either `is_selection_for` or `is_wt_control_for`
    * Condition 2: The ID of the section condition in the relationship. Must have been previously declared in the list of conditions.
8. **Time point definitions**: This fixed table is used to define time points in the experiment. At least one time point must be defined. The table has the following columns:
    * Time point name: The name of the timepoint. This should be a short identifier without spaces or special characters
    * Time: The numerical part of time point definition.
    * Unit: The time unit, such as "s" for seconds, "m" for minutes, "h" for hours, "d" for days, etc.
9. **Sequencing samples**: This final table contains the information of what the sequencing samples represent. It has the following columns:
    * Sample ID: This is the ID of the sample used in the name of the corresponding pair of FASTQ files. This can be a number or an alphanumerical label. No special characters (except "-" signs) are allowed.
    * Tile ID: The sequencing tile to which this sample belongs. This is a cross-reference to the "Tile Number" in the "Sequencing tiles" table and must have a matching entry there.
    * Condition: The condition to which this sample belongs. This is a cross-reference to the list of condition names and must have a matching entry there.
    * Time point: The time point to which this sample belongs. This is a cross-reference to the Time point definitions and must have a matching entry there.
    * Replicate: The replicate to which this sample belongs. This must be an integer number and must be within the range of replicates defined for the appropriate condition.

### Converting the parameter sheet
To convert the parameter sheet to JSON format you can use the csv2json.R script


