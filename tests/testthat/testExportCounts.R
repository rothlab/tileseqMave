test_that("toMavedbCountTable converts allCounts and marginalCounts tables", {
  allCounts <- data.frame(
    hgvsc = "c.1A>T",
    hgvsp = "p.Met1Leu",
    codonChanges = "ATG1TTG",
    codonHGVS = "c.[1A>T]",
    aaChanges = "M1L",
    aaChangeHGVS = "p.[Met1Leu]",
    sample.count = 5,
    sample.frequency = 0.5,
    sample.effectiveDepth = 10,
    stringsAsFactors = FALSE
  )

  allCountsMaveDB <- toMavedbCountTable(allCounts)

  expect_named(allCountsMaveDB, c(
    "hgvs_nt", "hgvs_pro", "sample.count", "sample.frequency",
    "sample.effectiveDepth"
  ))
  expect_equal(allCountsMaveDB[1, ], data.frame(
    hgvs_nt = "c.1A>T",
    hgvs_pro = "p.Met1Leu",
    sample.count = 5,
    sample.frequency = 0.5,
    sample.effectiveDepth = 10,
    stringsAsFactors = FALSE
  ))

  marginalCounts <- data.frame(
    hgvsc = c("c.1A>T", "c.4G>C"),
    hgvsp = c("p.Met1Leu", "p.Val2Leu"),
    codonChange = c("ATG1TTG", "GTG2CTG"),
    aaChange = c("M1L", "V2L"),
    sample.count = c(5, 7),
    sample.frequency = c(0.5, 0.7),
    sample.effectiveDepth = c(10, 10),
    stringsAsFactors = FALSE
  )

  marginalMaveDB <- toMavedbCountTable(marginalCounts)

  expect_named(marginalMaveDB, c(
    "hgvs_nt", "hgvs_pro", "sample.count", "sample.frequency",
    "sample.effectiveDepth"
  ))
  expect_equal(nrow(marginalMaveDB), 2)
  expect_equal(marginalMaveDB$sample.count, c(5, 7))
})

test_that("exportCounts writes both tables and rejects an empty conversion", {
  workspace <- tempfile("export-counts-")
  input <- file.path(workspace, "input")
  output <- file.path(workspace, "output")
  dir.create(input, recursive = TRUE)

  allCounts <- data.frame(
    hgvsc = "c.1A>T", hgvsp = "p.Met1Leu", codonChanges = "ATG1TTG",
    codonHGVS = "c.[1A>T]", aaChanges = "M1L", aaChangeHGVS = "p.[Met1Leu]",
    sample.count = 5, stringsAsFactors = FALSE
  )
  marginalCounts <- data.frame(
    hgvsc = "c.1A>T", hgvsp = "p.Met1Leu", codonChange = "ATG1TTG",
    aaChange = "M1L", sample.count = 5, sample.frequency = 0.5,
    sample.effectiveDepth = 10, stringsAsFactors = FALSE
  )
  write.csv(allCounts, file.path(input, "allCounts.csv"), row.names = FALSE)
  write.csv(marginalCounts, file.path(input, "marginalCounts.csv"), row.names = FALSE)

  script <- system.file("scripts", "exportCounts.R", package = "tileseqMave")
  parameters <- system.file("testdata", "parameters.json", package = "tileseqMave")
  rscript <- file.path(R.home("bin"), "Rscript")
  args <- c(
    script, "--workspace", workspace, "--input", input,
    "--output", output, "--parameters", parameters
  )
  result <- system2(rscript, shQuote(args), stdout = TRUE, stderr = TRUE)
  status <- attr(result, "status")
  if (is.null(status)) {
    status <- 0L
  }

  expect_equal(status, 0L)
  expect_true(any(grepl("Script completed successfully", result, fixed = TRUE)))
  expect_equal(
    read.csv(file.path(output, "allCounts_mavedb.csv")),
    toMavedbCountTable(allCounts)
  )
  marginal <- read.csv(file.path(output, "marginalCounts_mavedb.csv"))
  expect_equal(marginal, toMavedbCountTable(marginalCounts))
  expect_gt(ncol(marginal), 0)

  invalidInput <- file.path(workspace, "invalid")
  invalidOutput <- file.path(workspace, "invalid-output")
  dir.create(invalidInput)
  write.csv(
    data.frame(codonChanges = "ATG1TTG", aaChanges = "M1L"),
    file.path(invalidInput, "invalidCounts.csv"),
    row.names = FALSE
  )
  args <- c(
    script, "--workspace", workspace, "--input", invalidInput,
    "--output", invalidOutput, "--parameters", parameters
  )
  result <- suppressWarnings(
    system2(rscript, shQuote(args), stdout = TRUE, stderr = TRUE)
  )
  status <- attr(result, "status")

  expect_gt(status, 0L)
  expect_false(file.exists(file.path(invalidOutput, "invalid_mavedb.csv")))
  expect_false(any(grepl("Script completed successfully", result, fixed = TRUE)))
})
