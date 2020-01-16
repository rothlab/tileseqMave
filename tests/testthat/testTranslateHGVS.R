library(tileseqMave)
library(hgvsParseR)

context("HGVS translation")

test_that("translation works", {

	# paramfile <- "inst/testdata/parameters.json"
	paramfile <- system.file("testdata/parameters.json",
		package = "tileseqMave",
		mustWork = TRUE
	)
	
	parameters <- parseParameters(paramfile)

	builder <- new.hgvs.builder.p(aacode=3)

	expect_error(
		translateHGVS("c.1T>G",parameters,builder),
		"Reference mismatch!"
	)

	expect_equal(
		translateHGVS("c.1A>G",parameters,builder),
		"p.Met1Val"
	)

	expect_equal(
		translateHGVS("c.[1A>T;3G>T]",parameters,builder),
		"p.Met1Phe"
	)

	expect_equal(
		translateHGVS("c.[1A>T;10G>T]",parameters,builder),
		"p.[Met1Leu;Glu4Ter]"
	)

	expect_equal(
		translateHGVS("c.1_3delinsTTT",parameters,builder),
		"p.Met1Phe"
	)

	expect_equal(
		translateHGVS("c.2_3del",parameters,builder),
		"p.Met1fs"
	)

	expect_equal(
		translateHGVS("c.4_6del",parameters,builder),
		"p.Pro2del"
	)

	expect_equal(
		translateHGVS("c.4_9del",parameters,builder),
		"p.Pro2_Ser3del"
	)

	expect_equal(
		translateHGVS("c.4_9delinsTTTTTT",parameters,builder),
		"p.Pro2_Ser3delinsPhePhe"
	)

	expect_equal(
		translateHGVS("c.4_9delinsTTTTTTAAA",parameters,builder),
		"p.Pro2_Ser3delinsPhePheLys"
	)

	expect_equal(
		translateHGVS("c.[5_10delinsTTTTTT;11A>G]",parameters,builder),
		"p.Pro2_Glu4delinsLeuPheTrp"
	)

})

