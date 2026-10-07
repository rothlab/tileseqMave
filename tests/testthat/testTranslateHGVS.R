library(tileseqMave)
library(hgvsParseR)

context("HGVS translation")

test_that("translation works", {

	options(stringsAsFactors=FALSE)

	# paramfile <- "inst/testdata/parameters.json"
	paramfile <- system.file("testdata/parameters.json",
		package = "tileseqMave",
		mustWork = TRUE
	)
	
	parameters <- parseParameters(paramfile)
	cdsSeq <- parameters$template$cdsSeq

	builder <- new.hgvs.builder.p(aacode=3)

	expect_error(
		translateHGVS("c.1T>G",cdsSeq,builder),
		"Reference mismatch!"
	)

	expect_equal(
		translateHGVS("c.1A>G",cdsSeq,builder)[["hgvsp"]],
		"p.Met1Val"
	)

	expect_equal(
		translateHGVS("c.[1A>T;3G>T]",cdsSeq,builder)[["hgvsp"]],
		"p.Met1Phe"
	)

	expect_equal(
		translateHGVS("c.[3G>T;4C>T]",cdsSeq,builder)[["hgvsp"]],
		"p.Met1_Pro2delinsIleSer"
	)

	expect_equal(
		translateHGVS("c.[1A>T;10G>T]",cdsSeq,builder)[["hgvsp"]],
		"p.[Met1Leu;Glu4Ter]"
	)

	expect_equal(
		translateHGVS("c.1_3delinsTTT",cdsSeq,builder)[["hgvsp"]],
		"p.Met1Phe"
	)

	expect_equal(
		translateHGVS("c.2_3del",cdsSeq,builder)[["hgvsp"]],
		"p.Met1fs"
	)

	expect_equal(
		translateHGVS("c.4_6del",cdsSeq,builder)[["hgvsp"]],
		"p.Pro2del"
	)

	expect_equal(
		translateHGVS("c.4_9del",cdsSeq,builder)[["hgvsp"]],
		"p.Pro2_Ser3del"
	)

	expect_equal(
		translateHGVS("c.4_9delinsTTTTTT",cdsSeq,builder)[["hgvsp"]],
		"p.Pro2_Ser3delinsPhePhe"
	)

	expect_equal(
		translateHGVS("c.4_9delinsTTTTTTAAA",cdsSeq,builder)[["hgvsp"]],
		"p.Pro2_Ser3delinsPhePheLys"
	)

	expect_equal(
		translateHGVS("c.[5_10delinsTTTTTT;11A>G]",cdsSeq,builder)[["hgvsp"]],
		"p.Pro2_Glu4delinsLeuPheTrp"
	)

	expect_equal(
		translateHGVS("c.5_8delinsA",cdsSeq,builder)[["hgvsp"]],
		"p.Pro2_Ser3delinsHis"
	)
	
	expect_equal(
	  translateHGVS("c.0_1insA",cdsSeq,builder)[["codonChanges"]],
	  "silent"
	)
	
	expect_equal(
	  translateHGVS("c.[0_1insA;3_5delinsAAA]",cdsSeq,builder)[["hgvsp"]],
	  "p.Met1_Pro2delinsIleAsn"
	)
	
	

})

