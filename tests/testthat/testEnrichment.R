library(tileseqMave)

runEnrichment <- function(csvPath) {
  fx <- loadFixtureParams(csvPath)
  countDir <- tempfile("counts_")
  outDir <- tempfile("scores_")
  workDir <- tempfile("work_")
  dir.create(workDir)
  simulateCountFiles(fx$params, countDir)
  calcEnrichment(
    workDir, inDir = countDir, outDir = outDir, paramFile = fx$jsonFile,
    mc.cores = 1, nbs = 0, bcOverride = TRUE
  )
  list(params = fx$params, outDir = outDir)
}

test_that("mean.sd.count rejects a time point the select condition does not declare", {
  params <- loadFixtureParams(mixedFixturePath())$params
  # the time point check runs before the count tables are touched
  expect_error(
    tileseqMave:::mean.sd.count("select", NULL, NULL, "1", params),
    "Unrecognized time point 1 for selection condition select"
  )
})

test_that("calcEnrichment scores only time points declared for the select condition", {
  res <- runEnrichment(mixedFixturePath())
  expect_setequal(
    list.files(res$outDir, pattern = "_enrichment\\.csv$"),
    expectedEnrichmentFiles(res$params)
  )
  # explicit: the global time point the select does not declare is absent
  expect_false(file.exists(file.path(res$outDir, "select_t1_enrichment.csv")))
})

test_that("calcEnrichment output is unchanged in structure on a uniform time point sheet", {
  res <- runEnrichment(uniformFixturePath())
  expect_setequal(
    list.files(res$outDir, pattern = "_enrichment\\.csv$"),
    expectedEnrichmentFiles(res$params)
  )
})