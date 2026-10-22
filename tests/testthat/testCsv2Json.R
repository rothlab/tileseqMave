library(tileseqMave)
context("parameter parsing")

test_that("CSV2JSON", {

  # infile <- "test/test_that/paramtest.csv" # nolint: commented_code_linter.
  infile <- system.file(
    "testdata/paramtest.csv",
    package = "tileseqMave",
    mustWork = TRUE
  )
  outfile <- tempfile(fileext = ".json")
  csvParam2Json(infile, outfile)

  parameters <- parseParameters(outfile)
  expect_type(parameters, "list")

})

test_that(
  "allows nonselect t0 with multiple select timepoints",
  {

    infile <- test_path(
      "fixtures",
      "paramtest_nonselect_timepoint.csv"
    )

    outfile <- tempfile(fileext = ".json")

    expect_no_error(
      csvParam2Json(infile, outfile)
    )

    expect_true(file.exists(outfile))

    parameters <- parseParameters(outfile)
    expect_type(parameters, "list")
  }
)

test_that(
  "missing tile is still rejected for condition-specific timepoints",
  {

    infile <- test_path(
      "fixtures",
      "paramtest_missing_tile.csv"
    )

    outfile <- tempfile(fileext = ".json")

    expect_error(
      csvParam2Json(infile, outfile),
      "Missing samples for tiles"
    )
  }
)