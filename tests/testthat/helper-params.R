# Parse a parameter sheet CSV into a params object via the real parser.

loadFixtureParams <- function(csvPath) {
  jsonFile <- tempfile(fileext = ".json")
  csvParam2Json(csvPath, jsonFile)
  list(params = parseParameters(jsonFile), jsonFile = jsonFile)
}

mixedFixturePath <- function() {
  test_path("fixtures", "paramtest_nonselect_timepoint.csv")
}

uniformFixturePath <- function() {
  system.file("testdata/paramtest.csv", package = "tileseqMave", mustWork = TRUE)
}

# The enrichment files calcEnrichment should write: one per select
# condition and time point that the condition itself declares.
expectedEnrichmentFiles <- function(params) {
  unlist(lapply(getSelects(params), function(sCond) {
    sprintf("%s_t%s_enrichment.csv", sCond, getTimepointsFor(sCond, params))
  }))
}