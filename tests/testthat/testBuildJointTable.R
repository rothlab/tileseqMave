test_that("buildJointTable fails with a clear message when no count files exist", {
  d <- tempfile("tsm"); dir.create(d)
  inDir <- file.path(d, "x_mut_count"); dir.create(inDir)
  paramFile <- system.file("testdata", "parameters.json", package = "tileseqMave")
  expect_error(
    buildJointTable(d, inDir = inDir, outDir = file.path(d, "out"),
                    paramFile = paramFile),
    "No count files found"
  )
  unlink(d, recursive = TRUE)
})