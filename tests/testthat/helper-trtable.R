# Shared helpers for the trtable tests.

getInternalTrtable <- function() {
  get("trtable", envir = asNamespace("tileseqMave"), inherits = FALSE)
}

# hash object -> named character vector, sorted by codon
trtableAsVector <- function(h) {
  v <- vapply(hash::keys(h), function(k) h[[k]], character(1))
  v[order(names(v))]
}

testParamFile <- function() {
  system.file("testdata", "parameters.json", package = "tileseqMave")
}

# Leaked copies in the global environment would hide a missing internal table.
skipIfTrtableInGlobalEnv <- function() {
  testthat::skip_if(
    exists("trtable", envir = globalenv(), inherits = FALSE),
    "trtable already exists in the global environment"
  )
}

# TRUE if the function's code contains a data(trtable) call (including nested functions)
callsDataTrtable <- function(fnName) {
  f <- get(fnName, envir = asNamespace("tileseqMave"))
  any(grepl("data(trtable)", deparse(body(f)), fixed = TRUE))
}