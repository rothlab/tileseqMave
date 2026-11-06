test_that("simulatePool finds trtable and does not leak it", {
  skipIfTrtableInGlobalEnv()
  expect_false(callsDataTrtable("simulatePool"))

  params <- parseParameters(testParamFile())
  pool <- simulatePool(params)
  expect_false(exists("trtable", envir = globalenv(), inherits = FALSE))

  cc <- pool$allCodonChanges
  expect_true(all(c("from", "to", "fromAA", "toAA") %in% colnames(cc)))
  expect_false(anyNA(cc$toAA))
  expect_true("*" %in% cc$toAA)

  trtable <- getInternalTrtable()
  expect_equal(unname(cc$toAA),
               unname(vapply(as.character(cc$to), function(cd) trtable[[cd]], character(1))))
})

# First attempt: simulateFitness has not been run on the test parameter file before.
test_that("simulateFitness / calcTrueFitness can look up amino acids", {
  skipIfTrtableInGlobalEnv()
  expect_false(callsDataTrtable("simulateFitness"))

  set.seed(1)
  params <- parseParameters(testParamFile())
  pool <- simulatePool(params)
  fit <- simulateFitness(params, pool)

  ids <- as.character(pool$allCodonChanges$id[1:3])
  res <- fit$calcTrueFitness(c("WT", ids))
  expect_length(res, 4)
  expect_equal(res[[1]], 1)
  expect_false(anyNA(res))
  expect_true(all(res >= 0 & res <= 1))
  expect_false(exists("trtable", envir = globalenv(), inherits = FALSE))
})