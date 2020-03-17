# tileseqMave

Analysis functions for TileSEQ.

## Overview

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

1. Project name: This can be any name you like
2. Template construct: Here we define the sequencing template, its sequence and what it represents.
    * Gene name: The official HGNC gene name
    * Sequence: The full nucleotide sequence of the template including the coding sequence (CDS) and its flanking priming sequences.
    * CDS start: The position in the above sequence at which the coding sequence (CDS) begins. This must be an ATG codon! (Numbering starts at 1, not at 0; don't @ me, nerds. :P )
    * CDS end: The equivalent sequence position at which the CDS ends. This must be in-frame with the start position, i.e. end-start+1 must be divisible by 3.
    * Uniprot Accession: The UniprotKB accession of the protein encoded by the template gene.
3. Assay: A summary of the underlying selection assay
    * Assay Type: This can be any free-text label you like. Pick a simple name, like "Y2H", "Yeast complementation" or "LDL uptake via FACS"
    * Selection: This field indicates whether the assay performs a positive or negative selection. So only the values "Positive" and "Negative" are allowed. Positive selection indicates that the assay causes damaging variants to be depleted in the pool, where as negative selection indicates that the assay causes neutral variants to be depleted. (Most assays will be positive).
4. Conditions and replicates: This is a custom table that lists the different experimental conditions and their intended number of replicates and timepoints.



