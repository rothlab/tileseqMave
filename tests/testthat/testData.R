test_that("the internal translation table is the standard genetic code", {
  trtable <- getInternalTrtable()
  expect_length(hash::keys(trtable), 64)
  expect_equal(trtable[["ATG"]], "M")
  expect_equal(trtable[["TGG"]], "W")
  for (stopCodon in c("TAA", "TAG", "TGA")) {
    expect_equal(trtable[[stopCodon]], "*")
  }
  validAAs <- c(strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]], "*")
  expect_true(all(trtableAsVector(trtable) %in% validAAs))
})

test_that("the internal copy is stored in the namespace itself", {
  expect_true(exists("trtable", envir = asNamespace("tileseqMave"), inherits = FALSE))
})

test_that("the internal copy matches the dataset shipped in data/", {
  shipped <- tryCatch(getExportedValue("tileseqMave", "trtable"), error = function(e) NULL)
  skip_if(is.null(shipped), "lazy-loaded dataset not available in this session")
  expect_equal(trtableAsVector(getInternalTrtable()), trtableAsVector(shipped))
})