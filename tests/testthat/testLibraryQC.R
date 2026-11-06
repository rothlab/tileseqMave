test_that("libraryQC no longer calls data(trtable) and has the internal table available", {
  expect_false(callsDataTrtable("libraryQC"))
  expect_true(exists("trtable", envir = asNamespace("tileseqMave"), inherits = FALSE))
})