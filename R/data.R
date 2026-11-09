#' Codon translation table
#'
#' A lookup table for translating DNA codons into amino acids, following the
#' standard genetic code. It is stored as a \code{hash} object (see the
#' \pkg{hash} package), loaded with \code{data(trtable)}, and queried by codon,
#' for example \code{trtable[["ATG"]]}.
#'
#' @format A \code{hash} object with 64 entries. The keys are the 64 DNA codons
#'   (for example "ATG"), and the values are the corresponding one-letter amino
#'   acid codes, with "*" for stop codons.
#'
#' @keywords datasets
"trtable"