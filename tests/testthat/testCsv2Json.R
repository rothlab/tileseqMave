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

test_that("getTimepointsFor returns per-condition time points", {
  params <- list(samples = data.frame(
    Condition    = c("select", "select", "nonselect", "wtSel"),
    `Time point` = c("1", "2", "0", "0"),
    check.names  = FALSE
  ))
  expect_equal(getTimepointsFor("select", params), c("1", "2"))
  expect_equal(getTimepointsFor("nonselect", params), "0")
})

test_that("getTimepointsFor works when all conditions share one time point", {
  params <- list(samples = data.frame(
    Condition    = c("select", "nonselect"),
    `Time point` = c("1", "1"),
    check.names  = FALSE
  ))
  expect_equal(getTimepointsFor("select", params), "1")
  expect_equal(getTimepointsFor("nonselect", params), "1")
})

test_that("getTimepointsFor on an undeclared condition is NULL on mixed sheets", {
  # documents the behavior the calcEnrichment length(tps) == 0 guard relies on
  params <- list(samples = data.frame(
    Condition    = c("select", "select", "nonselect"),
    `Time point` = c("1", "2", "0"),
    check.names  = FALSE
  ))
  expect_null(getTimepointsFor("nope", params))
})

test_that("condition-specific time points are preserved when parsing a mixed sheet", {
  params <- loadFixtureParams(mixedFixturePath())$params

  expect_equal(getSelects(params), "select")
  expect_equal(getTimepointsFor("select", params), c("2", "3", "4"))
  expect_equal(getTimepointsFor("wtSelect", params), c("2", "3", "4"))
  expect_equal(getTimepointsFor("nonselect", params), "1")
  expect_equal(getTimepointsFor("wtNonselect", params), "1")
  # the precondition behind the calcEnrichment crash
  expect_true("1" %in% params$timepoints$`Time point name`)
  expect_false("1" %in% getTimepointsFor("select", params))
})