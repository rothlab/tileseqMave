library(tileseqMave)
context("parameter parsing")

test_that("CSV2JSON", {

	# infile <- "test/test_that/paramtest.csv"
	infile <- system.file("testdata/paramtest.csv",
		package = "tileseqMave",
		mustWork = TRUE
	)
	outfile <- tempfile(fileext=".json")
	csvParam2Json(infile,outfile)

	parameters <- parseParameters(outfile)	

	expect_type(parameters,"list")

})

