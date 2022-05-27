#!/usr/bin/env Rscript

# Copyright (C) 2021  Jochen Weile, Roth Lab
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

options(stringsAsFactors=FALSE)

library(yogilog)
library(httr)
library(RJSONIO)
library(yogitools)
library(argparser)

#process command line arguments
p <- arg_parser(
  "Automatically assemble reference variant sets via Clinvar and Gnomad",
  name="referenceSets.R"
)
p <- add_argument(p, "geneName", help="HGNC gene name (used to query clinvar and infer ensemblID)")
p <- add_argument(p, "--outfile", help="The desired output file location. Defaults to <geneName>_refVars.csv")
p <- add_argument(p, "--ensemblID", help="equivalent ensembl ID used to query gnomad.")
p <- add_argument(p, "--trait", help="trait name regular expression to filter clinvar.")
p <- add_argument(p, "--starsMin", help="minimum quality stars to filter clinvar.",default=2L)
p <- add_argument(p, "--homMin", help="minimum homozygous count to filter gnomad.",default=1L)
p <- add_argument(p, "--mafMin", help="minimum allele frequency to filter gnomad.",default=1e-3)
p <- add_argument(p, "--excludeLBLP", flag=TRUE, help="Excludes likely benign and likely pathogenic cases")
p <- add_argument(p, "--overrideCache", flag=TRUE, help="Override previously cached results.")
p <- add_argument(p, "--maxVars", help="Maximum number of clinvar variants to query",default=2000L)
p <- add_argument(p, "--logfile", help="The desired log file location.",default="referenceSets.log")

args <- parse_args(p)

if (is.na(args$outfile)) {
  outfile <- paste0(args$geneName,"_refVars.csv")
} else {
  outfile <- args$outfile
}

#set up logger and shunt it into the error handler
logger <- new.logger(args$logfile)
registerLogErrorHandler(logger)

#' Find cache file location by name
#' 
#' Finds the location for a cache file. The file does not necessary need to exist yet,
#' as this function is meant to be used determine to a location for both storage and retrieval.
#' 
#' Depending on the execution context, the storage location may differ. The cache location can 
#' be controlled with the environment variable \code{$MAVE_CACHE}.  If the variable is not set, 
#' a directory ".mavecache/" will be created in the user's home directory to be used as the 
#' storage location.
#' 
#' @param name the name of the file
#' @return the full path to the file
#' @export
#' @examples
#' file <- getCacheFile("P12345_alignment.fasta")
#' 
getCacheFile <- function(name) {
  cache.loc <- Sys.getenv("MAVE_CACHE",unset=NA)
  if (is.na(cache.loc)) {
    cache.loc <- paste0(Sys.getenv("HOME"),"/.mavecache/")
  }
  if (!file.exists(cache.loc)) {
    dir.create(cache.loc,showWarnings=FALSE,recursive=TRUE)
  }
  paste0(cache.loc,name)
}

#' Decode HTML strings
#' 
#' Decode HTML strings by converting special character escape sequences (e.g. \code{&gt;})
#' to their corresponding UTF-8 characters.
#' 
#' @param str input string
#' @return the decoded string
#' 
htmlDecode <- function(str) {
  str <- gsub("&gt;",">",str)
  #implement more as needed
  return(str)
}

subsets <- function(N,m) {
  if (m >= N) {
    return(list(1:N))
  }
  froms <- seq(1,N,m)
  mapply(function(from,to){from:to},froms,c(froms[-1],N))
}

fetchEnsemblID <- function(gene) {
  
  apiBase <- "https://rest.ensembl.org/xrefs/symbol/homo_sapiens/"

  library(httr)
  library(RJSONIO)
  library(yogitools)

  set_config(config(ssl_verifypeer = 0L))

  url <- paste0(apiBase,gene)
  htr <- GET(url,query=list(
    `content-type`="application/json"
  ))
  if (http_status(htr)$category == "Success") {
    #parse returned content as JSON
    returnData <- fromJSON(content(htr,as="text",encoding="UTF-8"))
    values <- do.call(c,lapply(returnData,function(entry) {
      if (entry[["type"]]=="gene" && grepl("ENSG",entry[["id"]])) {
        entry[["id"]]
      } else NULL
    }))
    return(values)

  } else {
    stop("Ensembl server message: ",http_status(htr)$message)
  }

}

#' Fetch missense variants in ClinVar for given gene
#' 
#' Makes a query to the NCBI webservice to fetch the ClinVar entries for the given gene
#' and filters them down to missense variants
#' 
#' @param gene gene name as in ClinVar (e.g. 'CALM1')
#' @param stagger logical determining whether HTTP requests should be staggered by
#'   a third of a second to avoid rejection by the server. Defaults to TRUE.
#' @param overrideCache logical determining whether local cache should be overridden
#' @param logger a yogilogger object to write log messages to. Defaults to NULL
#' @return a \code{data.frame} with the following columns:
#' \itemize{
#'   \item \code{hgvsc} The HGVS variant descriptor at the coding sequence (DNA) level
#'   \item \code{hgvsp} The HGVS variant descriptor at the amino acid sequence (Protein) level
#'   \item \code{clinsig} The ClinVar clinical significance string (e.g "Likely pathogenic")
#' }
#' @export
#' 
fetchClinvar <- function(gene,stagger=TRUE,overrideCache=FALSE,logger=NULL,maxVars=2000) {

  op <- options(stringsAsFactors=FALSE); on.exit(options(op))

  if (!is.null(logger)) {
    stopifnot(inherits(logger,"yogilogger"))
    library("yogilog")
  }

  cacheFile <- getCacheFile(paste0("clinvar_",gene,".csv"))

  if (!file.exists(cacheFile) || overrideCache) {

    library(httr)
    library(RJSONIO)
    library(yogitools)

    set_config(config(ssl_verifypeer = 0L))

    #This is a two-step process: First we have to query Clinvar for a list of matching
    # DB entries. Then we can make a query for the details of the matched entries.

    searchBase <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    summaryBase <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"

    clinvarIds <- character()

    if (stagger) {
      #stick to slightly less than 3 queries per second
      Sys.sleep(0.34)
    }

    if (!is.null(logger)) {
      logger$info("Querying Clinvar for ",gene)
    }
    #make HTTP GET request to find matching entries
    htr <- GET(searchBase,query=list(
      db="clinvar",
      term=paste0(gene,"[gene] AND single_gene[prop] AND ( clinsig_benign[prop] OR clinsig_pathogenic[prop] )"),
      retmax=maxVars,
      retmode="json"
    ))
    if (http_status(htr)$category == "Success") {
      #parse returned content as JSON
      returnData <- fromJSON(content(htr,as="text",encoding="UTF-8"))

      if (is.null(returnData) || length(returnData) < 1) {
        stop("Server returned no data.")
      }

      if (as.numeric(returnData$esearchresult$count) > maxVars) {
        warning("More than ",maxVars," results! Excess skipped.")
      }

      #extract the list of matches
      clinvarIds <- returnData$esearchresult$idlist

    } else {
      stop("Clinvar server message: ",http_status(htr)$message)
    }

    if (length(clinvarIds) < 1) {
      stop("No results!")
    }

    #Subdivide ids into smaller chunks if needed
    cidSets <- lapply(subsets(length(clinvarIds),200),function(js) clinvarIds[js])

    if (!is.null(logger)) {
      logger$info("Fetching detailed data from Clinvar")
    }

    results <- do.call(rbind,lapply(cidSets, function(clinvarIds) {

      if (stagger) {
        #stick to slightly less than 3 queries per second
        Sys.sleep(0.34)
      }

      #make HTTP get request for details of the matches
      cat(".")
      htr <- GET(summaryBase,query=list(
        db="clinvar",
        id=paste(clinvarIds,collapse=","),
        retmode="json"
      ))
      if (http_status(htr)$category == "Success") {
        #parse JSON
        returnData <- fromJSON(content(htr,as="text",encoding="UTF-8"))

        if (is.null(returnData) || length(returnData) < 1) {
          stop("Server returned no data.")
        }

        #the first result is simply a repetition of the entries, so we discard it.
        #the rest are the actual details on the entries.
        results <- as.df(lapply(returnData$result[-1], function(vset) {
          #extract the variant description
          # varStr <- htmlDecode(vset$variation_set[[1]]$cdna_change)
          varStr <- htmlDecode(vset$variation_set[[1]]$variation_name)
          trait <- htmlDecode(vset$trait_set[[1]]$trait_name)
          #and the clinical significance statement
          clinsig <- vset$clinical_significance[["description"]]
          quality <- vset$clinical_significance[["review_status"]]
          return(list(var=varStr,clinsig=clinsig,trait=trait,quality=quality))
        }))

        return(results)

      } else {
        stop("server message: ",http_status(htr)$message)
      }

    }))
    cat("\n")

    #filter the results down to only missense variants:
    #no must have protein-level HGVS string indicating AA change
    missense <- results[grepl("\\(p\\.\\w{3}\\d+\\w{3}\\)",results$var),]
    #must not be nonsense
    missense <- missense[!grepl("\\(p\\.\\w{3}\\d+Ter\\)",missense$var),]
    #must not be deletion
    missense <- missense[!grepl("\\(p\\.\\w{3}\\d+del\\)",missense$var),]
    #must not be frameshift
    missense <- missense[!grepl("\\(p\\.\\w{3}\\d+(\\w{3})?fs\\)",missense$var),]

    if (nrow(missense) > 0) {
      #extract the HGVS descriptors at the coding and protein levels.
      missense$hgvsc <- extract.groups(missense$var,"(c\\.\\d+[ACGT]>[ACGT]|c\\.\\d+_\\d+delins[ACGT]+)")[,1]
      missense$hgvsp <- extract.groups(missense$var,"(p\\.\\w{3}\\d+\\w{3})")[,1]
      #remove duplicates
      missense <- missense[!duplicated(missense$hgvsc),]
      output <- missense[,c("hgvsc","hgvsp","clinsig","trait","quality")]

    } else {
      output <- as.data.frame(matrix(nrow=0,ncol=5,dimnames=list(NULL,c(
        "hgvsc","hgvsp","clinsig","trait","quality"
      ))))
    }

    if (!is.null(logger)) {
      logger$info("Caching Clinvar data for ",gene)
    }
    write.table(output,cacheFile,sep=",")

  } else {
    if (!is.null(logger)) {
      logger$info("Retrieving cached data for ",gene)
    }
    output <- read.csv(cacheFile)
  }

  return(output)
}

clinvarStars <- function(qualities) {
  sapply(qualities, function(quality) {
    if (grepl("no assertion",quality)) return(0)
    if (grepl("criteria provided",quality) && grepl("conflicting|single submitter",quality)) return(1)
    if (grepl("multiple submitters",quality)) return(2)
    if (grepl("expert panel",quality)) return(3)
    if (grepl("guideline",quality)) return(4)
    else return(NA)
  })
}



#' Fetch missense variants in GnomAD for given gene
#' 
#' Makes a query to the ExAC webservice to fetch the GnomAD entries for the given gene
#' and filters them down to missense variants.
#' 
#' @param ensemblID The Ensembl gene identifier (e.g. \code{ENSG00000130164})
#' @param overrideCache logical determining whether local cache should be overridden
#' @param logger a yogilogger object to write log messages to. Defaults to NULL
#' @return a \code{data.frame} with the following columns:
#' \itemize{
#'   \item \code{hgvsc} The HGVS variant descriptor at the coding sequence (DNA) level
#'   \item \code{hgvsp} The HGVS variant descriptor at the amino acid sequence (Protein) level
#'   \item \code{maf} The minor allele frequency.
#' }
#' @export
#' 
fetchExac <- function(ensemblID,overrideCache=FALSE,logger=NULL) {

  op <- options(stringsAsFactors=FALSE); on.exit(options(op))

  if (!is.null(logger)) {
    stopifnot(inherits(logger,"yogilogger"))
    library("yogilog")
  }

  cacheFile <- getCacheFile(paste0("gnomad_",ensemblID,".csv"))

  if (!file.exists(cacheFile) || overrideCache) {

    library(httr)
    library(RJSONIO)
    library(yogitools)

    exacURL <- "http://exac.hms.harvard.edu/rest/gene/variants_in_gene/"

    if (!is.null(logger)) {
      logger$info("Querying GnomAD for ",ensemblID)
    }
    #make HTTP GET request
    htr <- GET(paste0(exacURL,ensemblID))
    if (http_status(htr)$category == "Success") {
      #parse JSON
      returnData <- fromJSON(content(htr,as="text",encoding="UTF-8"))

      #filter down to only canonical missense variants
      isMissense <- sapply(returnData,`[[`,"category")=="missense_variant"
      isCanonical <- sapply(returnData,`[[`,"CANONICAL")=="YES"
      missense.idx <- which(isCanonical & isMissense)

      missense.gnomad <- as.df(lapply(returnData[missense.idx],function(entry) {
        #extract hgvs and allele frequency data
        with(entry,list(
          hgvsc = HGVSc,
          hgvsp = HGVSp,
          maf = if (is.null(allele_freq)) NA else allele_freq,
          hom = if (is.null(hom_count)) 0 else hom_count
        ))
      }))

      #the annotation of missense is in gnomad is not correct, so we have to filter again!
      missense.gnomad <- missense.gnomad[grepl("^p\\.\\w{3}\\d+\\w{3}$",missense.gnomad$hgvsp),]
      #eliminate potential duplicates
      missense.gnomad <- missense.gnomad[!duplicated(missense.gnomad$hgvsc),]
    } else {
      stop("server message: ",http_status(htr)$message)
    }

    #Write results to cache
    if (!is.null(logger)) {
      logger$info("Caching GnomAD data for ",ensemblID)
    }
    write.table(missense.gnomad,cacheFile,sep=",")


  } else {
    if (!is.null(logger)) {
      logger$info("Retrieving cached data for ",ensemblID)
    }
    #eliminate potential duplicates
    missense.gnomad <- missense.gnomad[!duplicated(missense.gnomad$hgvsc),]
    missense.gnomad <- read.csv(cacheFile)
  }

  return(missense.gnomad)

}

#' Fetch missense variants in GnomAD for given gene
#' 
#' Makes a query to the ExAC webservice to fetch the GnomAD entries for the given gene
#' and filters them down to missense variants.
#' 
#' @param ensemblID The Ensembl gene identifier (e.g. \code{ENSG00000130164})
#' @param overrideCache logical determining whether local cache should be overridden
#' @param logger a yogilogger object to write log messages to. Defaults to NULL
#' @return a \code{data.frame} with the following columns:
#' \itemize{
#'   \item \code{hgvsc} The HGVS variant descriptor at the coding sequence (DNA) level
#'   \item \code{hgvsp} The HGVS variant descriptor at the amino acid sequence (Protein) level
#'   \item \code{maf} The minor allele frequency.
#' }
#' @export
#' 
fetchGnomad <- function(ensemblID,overrideCache=FALSE,logger=NULL) {

  op <- options(stringsAsFactors=FALSE); on.exit(options(op))

  if (!is.null(logger)) {
    stopifnot(inherits(logger,"yogilogger"))
    library("yogilog")
  }

  cacheFile <- getCacheFile(paste0("gnomad_",ensemblID,".csv"))

  if (!file.exists(cacheFile) || overrideCache) {

    library(httr)
    library(RJSONIO)
    library(yogitools)

    apiURL <- "https://gnomad.broadinstitute.org/api"

    if (!is.null(logger)) {
      logger$info("Querying GnomAD for ",ensemblID)
    }
    #build query (GraphQL)
    query <- sprintf("
{
  gene(gene_id: \"%s\", reference_genome: GRCh37) {
    variants(dataset: gnomad_r2_1_controls) {
      hgvsc
      hgvsp
      exome {
        homozygote_count
        af
      }
      genome {
        homozygote_count
        af
      }
      consequence
      consequence_in_canonical_transcript
    }
  }
}
",ensemblID)
    #make HTTP GET request
    htr <- POST(apiURL,encode="json",body=list(
      query=query,
      variables=list()
    ))
    if (http_status(htr)$category == "Success") {
      #parse JSON
      returnData <- fromJSON(content(htr,as="text",encoding="UTF-8"))

      #check for errors
      if (is.null(returnData$errors)) {
        vardata <- returnData[["data"]][["gene"]][["variants"]]

        #filter down to only canonical missense variants
        isMissense <- sapply(vardata,`[[`,"consequence")=="missense_variant"
        # isCanonical <- sapply(vardata,`[[`,"consequence_in_canonical_transcript")=="YES"
        missense.idx <- which(isMissense)

        missense.gnomad <- as.df(lapply(vardata[missense.idx],function(entry) {
          #extract hgvs and allele frequency data
          with(entry,list(
            hgvsc = hgvsc,
            hgvsp = if (!is.null(hgvsp)) hgvsp else "",
            maf = if (!is.null(exome)) {
              exome[["af"]]
            } else if (!is.null(genome)) {
              genome[["af"]]
            } else NA,
            hom = if (!is.null(exome)) {
              as.integer(exome[["homozygote_count"]])
            } else if (!is.null(genome)) {
              as.integer(genome[["homozygote_count"]])
            } else NA
          ))
        }))

        #the annotation of missense is in gnomad is not correct, so we have to filter again!
        missense.gnomad <- missense.gnomad[grepl("^p\\.\\w{3}\\d+\\w{3}$",missense.gnomad$hgvsp),]
        #eliminate potential duplicates
        missense.gnomad <- missense.gnomad[!duplicated(missense.gnomad$hgvsc),]
      } else {
        stop(returnData$errors[[1]])
      }
    } else {
      stop("server message: ",http_status(htr)$message)
    }

    #Write results to cache
    if (!is.null(logger)) {
      logger$info("Caching GnomAD data for ",ensemblID)
    }
    write.table(missense.gnomad,cacheFile,sep=",")


  } else {
    if (!is.null(logger)) {
      logger$info("Retrieving cached data for ",ensemblID)
    }
    missense.gnomad <- read.csv(cacheFile)
  }

  return(missense.gnomad)

}



#' Build reference set from Clinvar and GnomAD
#'  
#' @param geneName The plain gene name which will be used to query Clinvar
#' @param ensemblID The Ensembl gene identifier (e.g. \code{ENSG00000130164})
#' @param trait a trait name to filter by (or NA by default for no filter)
#' @param homMin minimum homozygous count in Gnomad to consider benign
#' @param mafMin minimum allele frequency in Gnomad to consider benign
#' @param starsMin minimum number of (quality) stars required in Clinvar
#' @param includeLikely whether to include likely benign and likely pathogenic variants
#' @param overrideCache logical determining whether local cache should be overridden
#' @param logger a yogilogger object to write log messages to. Defaults to NULL
#' @return a \code{data.frame} with the following columns:
#' \itemize{
#'   \item \code{hgvsc} The HGVS variant descriptor at the coding sequence (DNA) level
#'   \item \code{hgvsp} The HGVS variant descriptor at the amino acid sequence (Protein) level
#'   \item \code{maf} The minor allele frequency.
#'   \item \code{hom} The number of homozygous counts
#'   \item \code{referenceSet} Whether the variant is the positive or negative reference set
#' }
#' @export
#' 
buildReferenceSet <- function(geneName,ensemblID,
  trait=NA,homMin=1,mafMin=1e-3,starsMin=2,includeLikely=TRUE,
  stagger=TRUE,overrideCache=FALSE,logger=NULL,maxVars=2000) {

  op <- options(stringsAsFactors=FALSE); on.exit(options(op))

  clinvar <- fetchClinvar(geneName,stagger,overrideCache,logger,maxVars)
  gnomad <- fetchGnomad(ensemblID,overrideCache,logger)
  rownames(gnomad) <- gnomad$hgvsc

  #filter by clinvar trait
  if (!is.na(trait)) {
    clinvar <- clinvar[grep(trait,clinvar$trait,ignore.case=TRUE),]
  } else {
    #check if there's more than one trait
    traitTable <- table(clinvar$trait)
    if (length(traitTable > 1) && !is.null(logger)) {
      traitSummary <- paste(sprintf("%dx %s",traitTable,names(traitTable)),collapse=", ")
      logger$warn("More than one trait present in Clinvar results:",traitSummary)
    }
  }

  #filter by clinvar stars
  clinvar <- clinvar[clinvarStars(clinvar$quality) >= starsMin,]

  #split into positive and negative
  clinvarPatho <- clinvar[grep("pathogenic",clinvar$clinsig,ignore.case=includeLikely),1:2]
  clinvarBenign <- clinvar[grep("benign",clinvar$clinsig,ignore.case=includeLikely),1:2]

  #integrate MAF information where available
  clinvarPatho$maf=gnomad[clinvarPatho$hgvsc,"maf"]
  clinvarPatho$hom=gnomad[clinvarPatho$hgvsc,"hom"]
  clinvarBenign$maf=gnomad[clinvarBenign$hgvsc,"maf"]
  clinvarBenign$hom=gnomad[clinvarBenign$hgvsc,"hom"]

  #filter gnomad down based on given criteria
  if (homMin > 0) {
    gnomadBenign <- gnomad[with(gnomad,hom >= homMin | maf >= mafMin),]
  } else {
    gnomadBenign <- gnomad[with(gnomad,maf >= mafMin),]
  }
  #and exclude AA changes any that are explicitly pathogenic in clinvar
  gnomadBenign <- gnomadBenign[!(gnomadBenign$hgvsp %in% clinvarPatho$hgvsp),]
  #also, we don't need to include any cases from gnomad that are already benign in clinvar (i.e. redundant)
  gnomadBenign <- gnomadBenign[!(gnomadBenign$hgvsc %in% clinvarBenign$hgvsc),]

  #compile final output
  referenceSets <- data.frame(
    rbind(clinvarPatho,clinvarBenign,gnomadBenign),
    referenceSet=c(
      rep("Positive",nrow(clinvarPatho)),
      rep("Negative",nrow(clinvarBenign)+nrow(gnomadBenign))
    ),
    source=c(
      rep("Clinvar",nrow(clinvarPatho)+nrow(clinvarBenign)),
      rep("GnomAD",nrow(gnomadBenign))
    )
  )
  #reject any un-indexable variants
  referenceSets <- referenceSets[!is.na(referenceSets$hgvsc),]
  rownames(referenceSets) <- referenceSets$hgvsc

  return(referenceSets)
}


#if no ensembl id is provided by user, try to look it up via the gene name
if (is.na(args$ensemblID)) {
  logger$info("Fetching ensembl ID...")
  values <- fetchEnsemblID(args$geneName)
  if (length(values) ==0) {
    stop("Unable to find ensembl ID for gene name. Please provide one.")
  } else if (length(values) > 1) {
    stop("Gene name matches multiple ensembl IDs, please pick one one:",paste(values,", "))
  } else {
    args$ensemblID <- values[[1]]
  }
}

#run the actual process
referenceSets <- with(args,
  buildReferenceSet(geneName=geneName,ensemblID=ensemblID,
    trait=trait,homMin=homMin,mafMin=mafMin,starsMin=starsMin,includeLikely=!excludeLBLP,
    stagger=TRUE,overrideCache=overrideCache,logger=logger,maxVars=maxVars)
)

#write output to file
logger$info("Writing results to file ",outfile)
write.csv(referenceSets,outfile,row.names=FALSE)

logger$info("Success!")