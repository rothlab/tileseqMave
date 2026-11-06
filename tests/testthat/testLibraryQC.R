test_that("libraryQC no longer calls data(trtable) and has the internal table available", {
  expect_false(callsDataTrtable("libraryQC"))
  expect_true(exists("trtable", envir = asNamespace("tileseqMave"), inherits = FALSE))
})

test_that("drawNucleotideBias does not depend on its tile argument", {
  # minimal parameter sheet: only the tile numbers are read
  params <- list(tiles = cbind(`Tile Number` = c(1, 2)))
  # tile 1 has data (three SNVs, one MNV); tile 2 has none, so both the plotting
  # branch and the "no data" branch of the loop run
  sm <- data.frame(
    from  = c("AAA", "AAA", "CCC", "GGG"),
    pos   = c(1L, 1L, 2L, 3L),
    to    = c("AAG", "AAC", "CCA", "GTT"),
    freq  = c(0.01, 0.02, 0.03, 0.04),
    depth = 1000,
    tile  = 1,
    stringsAsFactors = FALSE
  )

  base <- tempfile("nucbias")
  makeInfo <- function(name) {
    dir <- file.path(base, name)
    dir.create(dir, recursive = TRUE)
    list(nsCond = "ns", tp = "1", outDir = paste0(dir, "/"), pdftag = "test")
  }
  pdfName <- "ns_t1_nucleotide_bias.pdf"

  # 1. new-style call, without the argument
  drawNucleotideBias(makeInfo("a"), params, sm)
  # 2. legacy call with an arbitrary value
  drawNucleotideBias(makeInfo("b"), params, sm, tile = 999)
  # 3. an argument that fails if the function ever evaluates it
  drawNucleotideBias(makeInfo("c"), params, sm,
                     tile = stop("tile must never be evaluated"))

  files <- file.path(base, c("a", "b", "c"), pdfName)
  expect_true(all(file.exists(files)))
  expect_true(all(file.size(files) > 0))
  headers <- vapply(files, function(f) rawToChar(readBin(f, "raw", 5)), character(1))
  expect_true(all(headers == "%PDF-"))
  # same input gives the same drawing; the date strings in the PDF have constant length
  expect_equal(file.size(files[[2]]), file.size(files[[1]]))
  expect_equal(file.size(files[[3]]), file.size(files[[1]]))

  graphics.off()
  unlink(base, recursive = TRUE)
})