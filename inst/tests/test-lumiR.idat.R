context("lumiR.idat suite")

test_that("iDAT files only", {
	library(lumi)
	options(width=80)
	
	path <- system.file("extdata", package="lumidat")
	files <- c("5356583020_A_Grn.idat", "5356583020_B_Grn.idat")
	manifestfile <- system.file("extdata", "HumanHT-12_V3_0_R1_99999999.txt", package="lumidat")
	x <- lumiR.idat(files, path, manifestfile=manifestfile, probeID="ProbeID")
	expect_that(x, is_a("LumiBatch"))
	expect_equivalent(dim(exprs(x)), c(90,2))
	
	res <- capture.output(x)[-c(7,12:19)]
	expected.res <- c("Summary of data information:", "\t Data File Information:", 
	"\t\tIllumina Inc. GenomeStudio version 1.8.0", "\t\tNormalization = none", 
	"\t\tArray Content = HumanHT-12_V3_0_R1_99999999.bgx.xml", "\t\tError Model = none", 
	"\t\tLocal Settings = en-US", "\t\t", "", "Major Operation History:", 
	"2         convertNuID = FALSE, lib.mapping = NULL, dec = \".\", parseColumnName = parseColumnName, ", 
	"3               checkDupId = checkDupId, QC = QC, columnNameGrepPattern = columnNameGrepPattern, ", 
	"4                        inputAnnotation = FALSE, annotationColumn = \"SYMBOL\", verbose = verbose)", 
	"5                            lumiQ(x.lumi = x.lumi, detectionTh = detectionTh, verbose = verbose)", 
	"  lumiVersion", "1      2.10.1", "2      2.10.1", "3      2.10.1", 
	"4      2.10.1", "5      2.10.1", "", "Object Information:", 
	"LumiBatch (storageMode: lockedEnvironment)", "assayData: 90 features, 2 samples ", 
	"  element names: beadNum, detection, exprs, se.exprs ", "protocolData: none", 
	"phenoData", "  sampleNames: 5356583020_A 5356583020_B", "  varLabels: sampleID", 
	"  varMetadata: labelDescription", "featureData", "  featureNames: ILMN_1704173 ILMN_1734742 ... ILMN_1755897 (90 total)", 
	"  fvarLabels: ProbeID TargetID", "  fvarMetadata: labelDescription", 
	"experimentData: use 'experimentData(object)'", "Annotation:  ", 
	"Control Data: Available", "QC information: Please run summary(x, 'QC') for details!"
	)
	expect_equivalent(res, expected.res)
	
	################################################################################
	res <- head(controlData(x))
	expected.res <- structure(list(controlType = c("biotin", "biotin", "cy3_hyb", 
	"cy3_hyb", "cy3_hyb", "cy3_hyb"), ProbeID = c("ILMN_1343048", 
	"ILMN_1343049", "ILMN_2038769", "ILMN_1343052", "ILMN_2038771", 
	"ILMN_2038770"), `5356583020_A` = c(16902.39, 20562.41, 17585.26, 
	410.2783, 3419.561, 15378.33), `5356583020_B` = c(16668.27, 20507.38, 
	15846.65, 407.3946, 3164.278, 15831.14)), .Names = c("controlType", 
	"ProbeID", "5356583020_A", "5356583020_B"), row.names = c("1", 
	"2", "3", "4", "5", "6"), class = "data.frame")
	expect_equivalent(res, expected.res)
	
	################################################################################
	
})

test_that("iDAT zipfile", {
	library(lumi)
	options(width=80)
	
	path <- system.file("extdata", package="lumidat")
	files <- c("5356583020_A_Grn.idat", "5356583020_B_Grn.idat")
	manifestfile <- system.file("extdata", "HumanHT-12_V3_0_R1_99999999.txt", package="lumidat")
	zipfile <- tempfile(fileext=".zip")
	zip(zipfile, file.path(path, files), flags = "-r9Xq")
	x <- lumiR.idat(zipfile=zipfile, manifestfile=manifestfile, probeID="ProbeID", verbose=FALSE)
	expect_that(x, is_a("LumiBatch"))
	expect_equivalent(dim(exprs(x)), c(90,2))
	
	res <- capture.output(x)[-c(7,12:19)]
	expected.res <- c("Summary of data information:", "\t Data File Information:", 
	"\t\tIllumina Inc. GenomeStudio version 1.8.0", "\t\tNormalization = none", 
	"\t\tArray Content = HumanHT-12_V3_0_R1_99999999.bgx.xml", "\t\tError Model = none", 
	"\t\tLocal Settings = en-US", "\t\t", "", "Major Operation History:", 
	"2         convertNuID = FALSE, lib.mapping = NULL, dec = \".\", parseColumnName = parseColumnName, ", 
	"3               checkDupId = checkDupId, QC = QC, columnNameGrepPattern = columnNameGrepPattern, ", 
	"4                        inputAnnotation = FALSE, annotationColumn = \"SYMBOL\", verbose = verbose)", 
	"5                            lumiQ(x.lumi = x.lumi, detectionTh = detectionTh, verbose = verbose)", 
	"  lumiVersion", "1      2.10.1", "2      2.10.1", "3      2.10.1", 
	"4      2.10.1", "5      2.10.1", "", "Object Information:", 
	"LumiBatch (storageMode: lockedEnvironment)", "assayData: 90 features, 2 samples ", 
	"  element names: beadNum, detection, exprs, se.exprs ", "protocolData: none", 
	"phenoData", "  sampleNames: 5356583020_A 5356583020_B", "  varLabels: sampleID", 
	"  varMetadata: labelDescription", "featureData", "  featureNames: ILMN_1704173 ILMN_1734742 ... ILMN_1755897 (90 total)", 
	"  fvarLabels: ProbeID TargetID", "  fvarMetadata: labelDescription", 
	"experimentData: use 'experimentData(object)'", "Annotation:  ", 
	"Control Data: Available", "QC information: Please run summary(x, 'QC') for details!")
	expect_equivalent(res, expected.res)
	
	################################################################################
	res <- head(controlData(x))
	expected.res <- structure(list(controlType = c("biotin", "biotin", "cy3_hyb", 
	"cy3_hyb", "cy3_hyb", "cy3_hyb"), ProbeID = c("ILMN_1343048", 
	"ILMN_1343049", "ILMN_2038769", "ILMN_1343052", "ILMN_2038771", 
	"ILMN_2038770"), `5356583020_A` = c(16902.39, 20562.41, 17585.26, 
	410.2783, 3419.561, 15378.33), `5356583020_B` = c(16668.27, 20507.38, 
	15846.65, 407.3946, 3164.278, 15831.14)), .Names = c("controlType", 
	"ProbeID", "5356583020_A", "5356583020_B"), row.names = c("1", 
	"2", "3", "4", "5", "6"), class = "data.frame")
	expect_equivalent(res, expected.res)
	
	################################################################################
	
})

test_that("idat files + CLM", {
	library(lumi)
	options(width=80)
	
	path <- system.file("extdata", package="lumidat")
	files <- c("5356583020_A_Grn.idat", "5356583020_B_Grn.idat")
	manifestfile <- system.file("extdata", "HumanHT-12_V3_0_R1_99999999.txt", package="lumidat")
	clmfile <- system.file("extdata", "5356583020.clm", package="lumidat")
	x <- lumiR.idat(files, path, manifestfile=manifestfile, clmfile=clmfile, probeID="ProbeID")
	expect_that(x, is_a("LumiBatch"))
	expect_equivalent(dim(exprs(x)), c(90,2))
	expect_equivalent(sampleNames(x), c("A", "B"))
	
	res <- grep("sampleNames", capture.output(x), value=TRUE)
	expected.res <- "  sampleNames: A B"
	expect_equivalent(res, expected.res)
	
	################################################################################
	res <- colnames(controlData(x))
	expected.res <- c("controlType", "ProbeID", "A", "B")
	expect_equivalent(res, expected.res)
})

################################################################################
################################################################################

# Detection pvalue slot
test_that("detection p-values are stable", {
	library(lumi)
	options(width=80)
	
	path <- system.file("extdata", package="lumidat")
	files <- c("5356583020_A_Grn.idat", "5356583020_B_Grn.idat")
	manifestfile <- system.file("extdata", "HumanHT-12_V3_0_R1_99999999.txt", package="lumidat")
	x <- lumiR.idat(files, path, manifestfile=manifestfile, probeID="ProbeID")
	expect_that(x, is_a("LumiBatch"))
	
	res <- detection(x)[1:15,]
	expected.res <- structure(c(0.0133333, 0, 0, 0, 0, 0, 0.9866667, 0.4666666, 0.0133333, 
	0, 0, 0.48, 0, 0, 0.0133333, 0.1333333, 0, 0, 0, 0, 0, 0.7066667, 
	0.6533333, 0, 0, 0, 0.8266667, 0, 0, 0), .Dim = c(15L, 2L), .Dimnames = list(
	    c("ILMN_1704173", "ILMN_1734742", "ILMN_1670130", "ILMN_2309245", 
	    "ILMN_1809344", "ILMN_1793729", "ILMN_2385662", "ILMN_2205695", 
	    "ILMN_1779257", "ILMN_1711453", "ILMN_1762281", "ILMN_1668162", 
	    "ILMN_1804174", "ILMN_1763663", "ILMN_1912287"), c("5356583020_A", 
	    "5356583020_B")))
	expect_equivalent(res, expected.res)
})


################################################################################
################################################################################

# probeID=c("ProbeID", "ArrayAddressID", "Sequence", "NuID"), 

test_that("probeID=ProbeID", {
	library(lumi)
	options(width=80)
	
	path <- system.file("extdata", package="lumidat")
	files <- c("5356583020_A_Grn.idat", "5356583020_B_Grn.idat")
	manifestfile <- system.file("extdata", "HumanHT-12_V3_0_R1_99999999.txt", package="lumidat")
	x <- lumiR.idat(files, path, manifestfile=manifestfile, probeID="ProbeID")
	expect_that(x, is_a("LumiBatch"))
	expect_equivalent(dim(exprs(x)), c(90,2))
	
	################################################################################
	res <- featureNames(x)[1:10]
	expected.res <- c("ILMN_1704173", "ILMN_1734742", "ILMN_1670130", "ILMN_2309245", 
	"ILMN_1809344", "ILMN_1793729", "ILMN_2385662", "ILMN_2205695", 
	"ILMN_1779257", "ILMN_1711453")
	expect_equivalent(res, expected.res)

	################################################################################
	res <- controlData(x)$ProbeID[1:5]
	expected.res <- c("ILMN_1343048", "ILMN_1343049", "ILMN_2038769", "ILMN_1343052", "ILMN_2038771")
	expect_equivalent(res, expected.res)
})

test_that("probeID=ArrayAddressID", {
	library(lumi)
	options(width=80)
	
	path <- system.file("extdata", package="lumidat")
	files <- c("5356583020_A_Grn.idat", "5356583020_B_Grn.idat")
	manifestfile <- system.file("extdata", "HumanHT-12_V3_0_R1_99999999.txt", package="lumidat")
	x <- lumiR.idat(files, path, manifestfile=manifestfile, probeID="ArrayAddressID")
	expect_that(x, is_a("LumiBatch"))
	expect_equivalent(dim(exprs(x)), c(90,2))
	
	################################################################################
	res <- featureNames(x)[1:10]
	expected.res <- c("4260551", "4760255", "6510553", "650553", "4250544", "870110", 
	"2100176", "1740731", "6420520", "6020544")
	expect_equivalent(res, expected.res)

	################################################################################
	res <- controlData(x)$ProbeID[1:5]
	expected.res <- c("5090180", "6510136", "1110170", "1450438", "2510500")
	expect_equivalent(res, expected.res)
})

test_that("probeID=Sequence", {
	library(lumi)
	options(width=80)
	
	path <- system.file("extdata", package="lumidat")
	files <- c("5356583020_A_Grn.idat", "5356583020_B_Grn.idat")
	manifestfile <- system.file("extdata", "HumanHT-12_V3_0_R1_99999999.txt", package="lumidat")
	x <- lumiR.idat(files, path, manifestfile=manifestfile, probeID="Sequence")
	expect_that(x, is_a("LumiBatch"))
	expect_equivalent(dim(exprs(x)), c(90,2))
	
	################################################################################
	res <- featureNames(x)[1:10]
	expected.res <- c("AGCTCACATGCAGTAGACTTGGGCAGGCAAAGGGGGCACCAAGGGCACAG", "CAGGCCTTCCCTGACCCAGCCAGGAACAAACAAGGGACCAAGTGCACACA", 
	"CCACACACTCACCACTCCCAGCTTCTCGTGTCCAGTGAAACCCCTGAACC", "AGGGGCGTTCTCCCAAAGATTAGGTCGTTTTCCAAAGAGCCGCGTCCCGG", 
	"TCACACAGCTAGGTCTTTTCAGAAGTGGTGGAAATTGGCAGCTGGGGTAC", "CTTGCCTAGAGAACACACATGGGCTTTGGAGCCCGACAGACCTGGGCTTG", 
	"GCTGAGAGAGAAAGAAGCAGCTCTTGAAGAAATGCGTAAGAAGATGCACC", "CCAGGCTGGTCTCAGAATTCTGGTGAGTGATCCTCCCACAGTGGACTTCC", 
	"CATGGATGCCAACCGGTCACCCAGGAGGATGGCAAAGAGAGTCGCATCTC", "GCCCACAAGGAGGACTTACAATGTGGAACTTAGTCTCTTTCCCTCACTCC")
	expect_equivalent(res, expected.res)

	################################################################################
	res <- controlData(x)$ProbeID[1:5]
	expected.res <- c(
		"GAATAAAGAACAATCTGCTGATGATCCCTCCGTGGATCTGATTCGTGTAA", "CCATGTGATACGAGGGCGCGTAGTTTGCATTATCGTTTTTATCGTTTCAA", 
		"AATTAAAACGATGCACACAGGGTTTAGCGCGTACACGTATTGCATTATGC", "TCTGTCACTGTCAGGAAAGTGGTAAAACTGCAACTCAATTACTGCAATGC", 
		"GCCCCGTATTCAGTGTCGCTGATTTGTATTGTCTGAAGTTGTTTTTACGT")
	expect_equivalent(res, expected.res)
})

test_that("probeID=NuID", {
	library(lumi)
	options(width=80)
	
	path <- system.file("extdata", package="lumidat")
	files <- c("5356583020_A_Grn.idat", "5356583020_B_Grn.idat")
	manifestfile <- system.file("extdata", "HumanHT-12_V3_0_R1_99999999.txt", package="lumidat")
	x <- lumiR.idat(files, path, manifestfile=manifestfile, probeID="NuID")
	expect_that(x, is_a("LumiBatch"))
	expect_equivalent(dim(exprs(x)), c(90,2))
	
	################################################################################
	res <- featureNames(x)[1:10]
	expected.res <- c("lJ0Tksh_pKQKqRQqRI", "NSl9XhUlKBAQqFC5EQ", "iURHRR1Sfdu1LgFXgU", 
	"QKpvdUCPK2.1AiWbVo", "Q0RJyt.0guugPpJ6rE", "0flyIEROp.olYSF6n4", 
	"ZniIgIJJ34IDmwgjkU", "9Up63SD3ri411RLofU", "QTo5QWtFSijpAiLZN0", 
	"NlRCih8Q7oHy3f1dHU")
	expect_equivalent(res, expected.res)

	################################################################################
	res <- controlData(x)$ProbeID[1:5]
	expected.res <- c("0gwIEN5441dbo3j27A", "BU7jGKmbL_Tzb.zb9A", "9DwBjkRKvyZsRs_Tzk", 
	"o3tHtKAusAeQdDx5Dk", "ulWz0u2eP7Pt4L7.xs")
	expect_equivalent(res, expected.res)
})

test_that("probeID=ERROR", {
	library(lumi)
	options(width=80)
	
	path <- system.file("extdata", package="lumidat")
	files <- c("5356583020_A_Grn.idat", "5356583020_B_Grn.idat")
	manifestfile <- system.file("extdata", "HumanHT-12_V3_0_R1_99999999.txt", package="lumidat")
	expect_error(x <- lumiR.idat(files, path, manifestfile=manifestfile, probeID="ERROR"))
})

################################################################################
################################################################################

# collapseMode=c("none", "max", "median", "mean"), 
test_that("collapseMode=none aka backgroundCorrect=FALSE", {
	library(lumi)
	options(width=80)
	
	path <- system.file("extdata", package="lumidat")
	files <- c("5356583020_A_Grn.idat", "5356583020_B_Grn.idat")
	manifestfile <- system.file("extdata", "HumanHT-12_V3_0_R1_99999999.txt", package="lumidat")
	x <- lumiR.idat(files, path, manifestfile=manifestfile, probeID="ProbeID", collapseMode="none")
	expect_that(x, is_a("LumiBatch"))
	expect_equivalent(dim(exprs(x)), c(90,2))
	
	# these probes are from genes with >1 probes
	probes <- c("ILMN_1744897", "ILMN_1801608", "ILMN_1711283", "ILMN_2296644", 
	"ILMN_1652609", "ILMN_1807467", "ILMN_1761911", "ILMN_2121316", 
	"ILMN_1810533", "ILMN_1659166", "ILMN_1679194", "ILMN_1755897")
	x <- x[probes, ]
	
	# > data.frame(as(featureData(x), "data.frame"),exprs(x), check.names=FALSE)
	#                   ProbeID TargetID 5356583020_A 5356583020_B
	# ILMN_1744897 ILMN_1744897    KCNN3     52.41120     60.19562
	# ILMN_1801608 ILMN_1801608    KCNN3     69.06544     73.10770
	# ILMN_1711283 ILMN_1711283  PCDHGA9     67.22208     51.06295
	# ILMN_2296644 ILMN_2296644  PCDHGA9     58.62053     65.66985
	# ILMN_1652609 ILMN_1652609    RGNEF    112.96380     82.18444
	# ILMN_1807467 ILMN_1807467    RGNEF     57.53606     73.94498
	# ILMN_1761911 ILMN_1761911 SCYL1BP1    187.19250    227.17070
	# ILMN_2121316 ILMN_2121316 SCYL1BP1    142.36550    138.79320
	# ILMN_1810533 ILMN_1810533  SLC6A15     76.14130     66.54485
	# ILMN_1659166 ILMN_1659166  SLC6A15     58.84699     55.11340
	# ILMN_1679194 ILMN_1679194   UGT2B7     58.90470     66.43938
	# ILMN_1755897 ILMN_1755897   UGT2B7     62.35626     67.29924
	
	################################################################################
	res <- exprs(x)[1:6,]
	expected.res <- structure(c(52.4112, 69.06544, 67.22208, 58.62053, 112.9638, 
	57.53606, 60.19562, 73.1077, 51.06295, 65.66985, 82.18444, 73.94498
	), .Dim = c(6L, 2L), .Dimnames = list(c("ILMN_1744897", "ILMN_1801608", 
	"ILMN_1711283", "ILMN_2296644", "ILMN_1652609", "ILMN_1807467"
	), c("5356583020_A", "5356583020_B")))
	expect_equivalent(res, expected.res)

	################################################################################
	res <- controlData(x)[1:6,]
	expected.res <- structure(list(controlType = c("biotin", "biotin", "cy3_hyb", 
	"cy3_hyb", "cy3_hyb", "cy3_hyb"), ProbeID = c("ILMN_1343048", 
	"ILMN_1343049", "ILMN_2038769", "ILMN_1343052", "ILMN_2038771", 
	"ILMN_2038770"), `5356583020_A` = c(16902.39, 20562.41, 17585.26, 
	410.2783, 3419.561, 15378.33), `5356583020_B` = c(16668.27, 20507.38, 
	15846.65, 407.3946, 3164.278, 15831.14)), .Names = c("controlType", 
	"ProbeID", "5356583020_A", "5356583020_B"), row.names = c("1", 
	"2", "3", "4", "5", "6"), class = "data.frame")
	expect_equivalent(res, expected.res)
})


test_that("collapseMode=max", {
	library(lumi)
	options(width=80)
	
	path <- system.file("extdata", package="lumidat")
	files <- c("5356583020_A_Grn.idat", "5356583020_B_Grn.idat")
	manifestfile <- system.file("extdata", "HumanHT-12_V3_0_R1_99999999.txt", package="lumidat")
	expect_message(
		x <- lumiR.idat(files, path, manifestfile=manifestfile, probeID="ProbeID", collapseMode="max"), 
		"addControlData2lumi can't handle data that's been collapsed. forcing controls=FALSE"
	)
	expect_that(x, is_a("LumiBatch"))
	expect_equivalent(dim(exprs(x)), c(82,2))
	
	# subset to the genes which have >1 probes
	x <- x[c("KCNN3", "NOMO1", "PCDHGA9", "RGNEF", "SCYL1BP1", "SLC6A15", "UGT2B7"), ]
	# here's the data from collapseMode=none; i'll choose the right probe & double check below
	# > data.frame(as(featureData(x), "data.frame"),exprs(x), check.names=FALSE)
	#                   ProbeID TargetID 5356583020_A 
	# ILMN_1744897 ILMN_1744897    KCNN3     52.41120 
	# ILMN_1801608 ILMN_1801608    KCNN3     69.06544 <<<<
	# ILMN_1711283 ILMN_1711283  PCDHGA9     67.22208 <<<< 
	# ILMN_2296644 ILMN_2296644  PCDHGA9     58.62053 
	# ILMN_1652609 ILMN_1652609    RGNEF    112.96380 <<<< 
	# ILMN_1807467 ILMN_1807467    RGNEF     57.53606 
	# ILMN_1761911 ILMN_1761911 SCYL1BP1    187.19250 <<<< 
	# ILMN_2121316 ILMN_2121316 SCYL1BP1    142.36550 
	# ILMN_1810533 ILMN_1810533  SLC6A15     76.14130 <<<< 
	# ILMN_1659166 ILMN_1659166  SLC6A15     58.84699 
	# ILMN_1679194 ILMN_1679194   UGT2B7     58.90470 
	# ILMN_1755897 ILMN_1755897   UGT2B7     62.35626 <<<< 
	
	################################################################################
	res <- exprs(x)
	expected.res <- structure(c(69.06544, 402.4396, 67.22208, 112.9638, 187.1925, 
	76.1413, 62.35626, 73.1077, 251.4016, 65.66985, 82.18444, 227.1707, 
	66.54485, 67.29924), .Dim = c(7L, 2L), .Dimnames = list(c("KCNN3", 
	"NOMO1", "PCDHGA9", "RGNEF", "SCYL1BP1", "SLC6A15", "UGT2B7"), 
	c("5356583020_A", "5356583020_B")))
	#          5356583020_A 5356583020_B
	# KCNN3        69.06544     73.10770
	# NOMO1       402.43960    251.40160
	# PCDHGA9      67.22208     65.66985
	# RGNEF       112.96380     82.18444
	# SCYL1BP1    187.19250    227.17070
	# SLC6A15      76.14130     66.54485
	# UGT2B7       62.35626     67.29924
	expect_equivalent(res, expected.res)
	# ^^^^ agrees with my manual assessment
	################################################################################
	res <- detection(x)
	expected.res <- structure(c(0.0933333, 0, 0.1066667, 0, 0, 0.0266666, 0.1866667, 
	0.0266666, 0, 0.0933333, 0, 0, 0.0933333, 0.0933333), .Dim = c(7L, 
	2L), .Dimnames = list(c("KCNN3", "NOMO1", "PCDHGA9", "RGNEF", 
	"SCYL1BP1", "SLC6A15", "UGT2B7"), c("5356583020_A", "5356583020_B"
	)))
	#          5356583020_A 5356583020_B
	# KCNN3       0.0933333    0.0266666
	# NOMO1       0.0000000    0.0000000
	# PCDHGA9     0.1066667    0.0933333
	# RGNEF       0.0000000    0.0000000
	# SCYL1BP1    0.0000000    0.0000000
	# SLC6A15     0.0266666    0.0933333
	# UGT2B7      0.1866667    0.0933333
	expect_equivalent(res, expected.res)
	
})

test_that("collapseMode=median", {
	library(lumi)
	options(width=80)
	
	path <- system.file("extdata", package="lumidat")
	files <- c("5356583020_A_Grn.idat", "5356583020_B_Grn.idat")
	manifestfile <- system.file("extdata", "HumanHT-12_V3_0_R1_99999999.txt", package="lumidat")
	expect_message(
		x <- lumiR.idat(files, path, manifestfile=manifestfile, probeID="ProbeID", collapseMode="median"), 
		"addControlData2lumi can't handle data that's been collapsed. forcing controls=FALSE"
	)
	expect_that(x, is_a("LumiBatch"))
	expect_equivalent(dim(exprs(x)), c(82,2))
	
	# subset to the genes which have >1 probes
	x <- x[c("KCNN3", "NOMO1", "PCDHGA9", "RGNEF", "SCYL1BP1", "SLC6A15", "UGT2B7"), ]
	# here's the data from collapseMode=none; i'll choose the right probe & double check below
	# > data.frame(as(featureData(x), "data.frame"),exprs(x), check.names=FALSE)
	#                   ProbeID TargetID 5356583020_A 
	# ILMN_1744897 ILMN_1744897    KCNN3     52.41120
	# ILMN_1801608 ILMN_1801608    KCNN3     69.06544
	# ILMN_1711283 ILMN_1711283  PCDHGA9     67.22208 
	# ILMN_2296644 ILMN_2296644  PCDHGA9     58.62053
	# ILMN_1652609 ILMN_1652609    RGNEF    112.96380 
	# ILMN_1807467 ILMN_1807467    RGNEF     57.53606
	# ILMN_1761911 ILMN_1761911 SCYL1BP1    187.19250 
	# ILMN_2121316 ILMN_2121316 SCYL1BP1    142.36550
	# ILMN_1810533 ILMN_1810533  SLC6A15     76.14130 
	# ILMN_1659166 ILMN_1659166  SLC6A15     58.84699
	# ILMN_1679194 ILMN_1679194   UGT2B7     58.90470
	# ILMN_1755897 ILMN_1755897   UGT2B7     62.35626 
	
	################################################################################
	res <- exprs(x)
	expected.res <- structure(c(60.73832, 381.8736, 62.9213, 85.24995, 164.779, 67.49415, 
	60.63048, 66.65166, 224.6189, 58.3664, 78.06471, 182.982, 60.82912, 
	66.86931), .Dim = c(7L, 2L), .Dimnames = list(c("KCNN3", "NOMO1", 
	"PCDHGA9", "RGNEF", "SCYL1BP1", "SLC6A15", "UGT2B7"), c("5356583020_A", 
	"5356583020_B")))
	#          5356583020_A 5356583020_B
	# KCNN3        60.73832     66.65166
	# NOMO1       381.87360    224.61890
	# PCDHGA9      62.92130     58.36640
	# RGNEF        85.24995     78.06471
	# SCYL1BP1    164.77900    182.98200
	# SLC6A15      67.49415     60.82912
	# UGT2B7       60.63048     66.86931
	expect_equivalent(res, expected.res)
	## ah ok, collapseMode=median reverts to the 'arithmetic mean' of the middle 2 elements.
	
	################################################################################
	## same as for collapseMode=max... is that right??
	res <- detection(x)
	expected.res <- structure(c(0.0933333, 0, 0.1066667, 0, 0, 0.0266666, 0.1866667, 
	0.0266666, 0, 0.0933333, 0, 0, 0.0933333, 0.0933333), .Dim = c(7L, 
	2L), .Dimnames = list(c("KCNN3", "NOMO1", "PCDHGA9", "RGNEF", 
	"SCYL1BP1", "SLC6A15", "UGT2B7"), c("5356583020_A", "5356583020_B"
	)))
	#          5356583020_A 5356583020_B
	# KCNN3       0.0933333    0.0266666
	# NOMO1       0.0000000    0.0000000
	# PCDHGA9     0.1066667    0.0933333
	# RGNEF       0.0000000    0.0000000
	# SCYL1BP1    0.0000000    0.0000000
	# SLC6A15     0.0266666    0.0933333
	# UGT2B7      0.1866667    0.0933333
	expect_equivalent(res, expected.res)
	
})

test_that("collapseMode=mean", {
	library(lumi)
	options(width=80)
	
	path <- system.file("extdata", package="lumidat")
	files <- c("5356583020_A_Grn.idat", "5356583020_B_Grn.idat")
	manifestfile <- system.file("extdata", "HumanHT-12_V3_0_R1_99999999.txt", package="lumidat")
	expect_message(
		x <- lumiR.idat(files, path, manifestfile=manifestfile, probeID="ProbeID", collapseMode="mean"), 
		"addControlData2lumi can't handle data that's been collapsed. forcing controls=FALSE"
	)
	expect_that(x, is_a("LumiBatch"))
	expect_equivalent(dim(exprs(x)), c(82,2))
	
	# subset to the genes which have >1 probes
	x <- x[c("KCNN3", "NOMO1", "PCDHGA9", "RGNEF", "SCYL1BP1", "SLC6A15", "UGT2B7"), ]
	
	################################################################################
	res <- exprs(x)
	expected.res <- structure(c(60.73832, 286.5653, 62.9213, 85.24995, 164.779, 67.49415, 
	60.63048, 66.65166, 182.4367, 58.3664, 78.06471, 182.982, 60.82912, 
	66.86931), .Dim = c(7L, 2L), .Dimnames = list(c("KCNN3", "NOMO1", 
	"PCDHGA9", "RGNEF", "SCYL1BP1", "SLC6A15", "UGT2B7"), c("5356583020_A", 
	"5356583020_B")))
	#          5356583020_A 5356583020_B
	# KCNN3        60.73832     66.65166
	# NOMO1       286.56530    182.43670
	# PCDHGA9      62.92130     58.36640
	# RGNEF        85.24995     78.06471
	# SCYL1BP1    164.77900    182.98200
	# SLC6A15      67.49415     60.82912
	# UGT2B7       60.63048     66.86931
	expect_equivalent(res, expected.res)
	
	################################################################################
	## same as for collapseMode=max... is that right??
	res <- detection(x)
	expected.res <- structure(c(0.4333333, 0.00888889999999998, 0.2666667, 0.2333333, 
	0, 0.2133333, 0.2933333, 0.12, 0.00888889999999998, 0.42, 0.0133333, 
	0, 0.3, 0.0933333), .Dim = c(7L, 2L), .Dimnames = list(c("KCNN3", 
	"NOMO1", "PCDHGA9", "RGNEF", "SCYL1BP1", "SLC6A15", "UGT2B7"), 
	    c("5356583020_A", "5356583020_B")))
	#          5356583020_A 5356583020_B
	# KCNN3       0.4333333    0.1200000
	# NOMO1       0.0088889    0.0088889
	# PCDHGA9     0.2666667    0.4200000
	# RGNEF       0.2333333    0.0133333
	# SCYL1BP1    0.0000000    0.0000000
	# SLC6A15     0.2133333    0.3000000
	# UGT2B7      0.2933333    0.0933333
	expect_equivalent(res, expected.res)
	# this **is** different to collapseMode=mean, great.
})


################################################################################
################################################################################

test_that("backgroundCorrect=TRUE", {
	library(lumi)
	options(width=80)
	
	path <- system.file("extdata", package="lumidat")
	files <- c("5356583020_A_Grn.idat", "5356583020_B_Grn.idat")
	manifestfile <- system.file("extdata", "HumanHT-12_V3_0_R1_99999999.txt", package="lumidat")
	x <- lumiR.idat(files, path, manifestfile=manifestfile, probeID="ProbeID", backgroundCorrect=TRUE)
	expect_that(x, is_a("LumiBatch"))
	expect_equivalent(dim(exprs(x)), c(90,2))
	
	################################################################################
	res <- exprs(x)[1:6,]
	expected.res <- structure(c(30.52452, 1113.839, 1028.59, 724.6708, 394.9726, 
	1315.341, 6.253441, 802.3579, 902.5901, 507.5767, 466.2357, 1211.94
	), .Dim = c(6L, 2L), .Dimnames = list(c("ILMN_1704173", "ILMN_1734742", 
	"ILMN_1670130", "ILMN_2309245", "ILMN_1809344", "ILMN_1793729"
	), c("5356583020_A", "5356583020_B")))
	expect_equivalent(res, expected.res)

	################################################################################
	res <- controlData(x)[1:6,]
	expected.res <- structure(list(controlType = c("biotin", "biotin", "cy3_hyb", 
	"cy3_hyb", "cy3_hyb", "cy3_hyb"), ProbeID = c("ILMN_1343048", 
	"ILMN_1343049", "ILMN_2038769", "ILMN_1343052", "ILMN_2038771", 
	"ILMN_2038770"), `5356583020_A` = c(16902.39, 20562.41, 17585.26, 
	410.2783, 3419.561, 15378.33), `5356583020_B` = c(16668.27, 20507.38, 
	15846.65, 407.3946, 3164.278, 15831.14)), .Names = c("controlType", 
	"ProbeID", "5356583020_A", "5356583020_B"), row.names = c("1", 
	"2", "3", "4", "5", "6"), class = "data.frame")
	expect_equivalent(res, expected.res)
	
})

################################################################################
################################################################################

# forbidden params: c("convertNuID", "lib.mapping", "dec", "inputAnnotation", "annotationColumn"),
test_that("forbidden parameters throw errors", {
	library(lumi)
	options(width=80)
	
	path <- system.file("extdata", package="lumidat")
	files <- c("5356583020_A_Grn.idat", "5356583020_B_Grn.idat")
	manifestfile <- system.file("extdata", "HumanHT-12_V3_0_R1_99999999.txt", package="lumidat")
	
	################################################################################
	expect_error(x <- lumiR.idat(files, path, manifestfile=manifestfile, convertNuID=TRUE))
	expect_error(x <- lumiR.idat(files, path, manifestfile=manifestfile, lib.mapping="garbage"))
	expect_error(x <- lumiR.idat(files, path, manifestfile=manifestfile, dec="."))
	expect_error(x <- lumiR.idat(files, path, manifestfile=manifestfile, inputAnnotation="garbage"))
	expect_error(x <- lumiR.idat(files, path, manifestfile=manifestfile, annotationColumn="garbage"))
})

################################################################################
################################################################################

# # lumi::lumiR arguments:
# detectionTh=0.01, na.rm=TRUE, parseColumnName=FALSE, checkDupId=TRUE, QC=TRUE, 
# columnNameGrepPattern=list(exprs='AVG_SIGNAL', se.exprs='BEAD_STD', detection='DETECTION', beadNum='Avg_NBEADS'), 

## this seem to match the A or B from within the array barcodes, eg 5356583020_A_AVG_SIGNAL
test_that("parseColumnName=TRUE", {
	library(lumi)
	options(width=80)
	
	path <- system.file("extdata", package="lumidat")
	files <- c("5356583020_A_Grn.idat", "5356583020_B_Grn.idat")
	manifestfile <- system.file("extdata", "HumanHT-12_V3_0_R1_99999999.txt", package="lumidat")
	
	x <- lumiR.idat(files, path, manifestfile=manifestfile, parseColumnName=TRUE)
	res <- sampleNames(x)
	expected.res <- c("A", "B")
	expect_equivalent(res, expected.res)
	
})

# does nothing with these input files
test_that("checkDupId=TRUE", {
	library(lumi)
	options(width=80)
	
	path <- system.file("extdata", package="lumidat")
	files <- c("5356583020_A_Grn.idat", "5356583020_B_Grn.idat")
	manifestfile <- system.file("extdata", "HumanHT-12_V3_0_R1_99999999.txt", package="lumidat")
	
	x <- lumiR.idat(files, path, manifestfile=manifestfile, checkDupId=TRUE)
	expect_that(x, is_a("LumiBatch"))
})


# boring too...
test_that("QC=FALSE", {
	library(lumi)
	options(width=80)
	
	path <- system.file("extdata", package="lumidat")
	files <- c("5356583020_A_Grn.idat", "5356583020_B_Grn.idat")
	manifestfile <- system.file("extdata", "HumanHT-12_V3_0_R1_99999999.txt", package="lumidat")
	
	x <- lumiR.idat(files, path, manifestfile=manifestfile, QC=FALSE)
	expect_that(x, is_a("LumiBatch"))
})

# columnNameGrepPattern=list(exprs='AVG_SIGNAL', se.exprs='BEAD_STD', detection='DETECTION', beadNum='Avg_NBEADS')
test_that("columnNameGrepPattern=MIN_SIGNAL", {
	library(lumi)
	options(width=80)
	
	path <- system.file("extdata", package="lumidat")
	files <- c("5356583020_A_Grn.idat", "5356583020_B_Grn.idat")
	manifestfile <- system.file("extdata", "HumanHT-12_V3_0_R1_99999999.txt", package="lumidat")
	
	x <- lumiR.idat(files, path, manifestfile=manifestfile, columnNameGrepPattern=list(exprs='MIN_SIGNAL', se.exprs='BEAD_STD', detection='DETECTION', beadNum='Avg_NBEADS'))
	expect_that(x, is_a("LumiBatch"))
	
	res <- exprs(x)[1:6,]
	expected.res <- structure(c(88.54737, 1171.861, 1086.613, 782.6937, 452.9955, 
	1373.364, 62.09662, 858.2011, 958.4333, 563.4199, 522.0789, 1267.783
	), .Dim = c(6L, 2L), .Dimnames = list(c("ILMN_1704173", "ILMN_1734742", 
	"ILMN_1670130", "ILMN_2309245", "ILMN_1809344", "ILMN_1793729"
	), c("5356583020_A", "5356583020_B")))
	expect_equivalent(res, expected.res)
	
	# @TODO deliberate error to make MJC consider:
	# why are MIN_SIGNAL,MAX_Signal,AVG_SIGNAL all identical? is it just this test data??
	expect_equivalent(TRUE, FALSE)
	
})

test_that("columnNameGrepPattern=MAX_Signal", {
	library(lumi)
	options(width=80)
	
	path <- system.file("extdata", package="lumidat")
	files <- c("5356583020_A_Grn.idat", "5356583020_B_Grn.idat")
	manifestfile <- system.file("extdata", "HumanHT-12_V3_0_R1_99999999.txt", package="lumidat")
	
	x <- lumiR.idat(files, path, manifestfile=manifestfile, columnNameGrepPattern=list(exprs='MAX_Signal', se.exprs='BEAD_STD', detection='DETECTION', beadNum='Avg_NBEADS'))
	expect_that(x, is_a("LumiBatch"))
	
	res <- exprs(x)[1:6,]
	expected.res <- structure(c(88.54737, 1171.861, 1086.613, 782.6937, 452.9955, 
	1373.364, 62.09662, 858.2011, 958.4333, 563.4199, 522.0789, 1267.783
	), .Dim = c(6L, 2L), .Dimnames = list(c("ILMN_1704173", "ILMN_1734742", 
	"ILMN_1670130", "ILMN_2309245", "ILMN_1809344", "ILMN_1793729"
	), c("5356583020_A", "5356583020_B")))
	expect_equivalent(res, expected.res)
	
	# @TODO deliberate error to make MJC consider:
	# why are MIN_SIGNAL,MAX_Signal,AVG_SIGNAL all identical? is it just this test data??
	expect_equivalent(TRUE, FALSE)
})

# # dropping the final terms in columnNameGrepPattern don't have an effet, at least on detection(), se.exprs(),
# # so default values must be enforced
# test_that("columnNameGrepPattern simple, exprs only", {
# 	library(lumi)
# 	options(width=80)
# 	
# 	path <- system.file("extdata", package="lumidat")
# 	files <- c("5356583020_A_Grn.idat", "5356583020_B_Grn.idat")
# 	manifestfile <- system.file("extdata", "HumanHT-12_V3_0_R1_99999999.txt", package="lumidat")
# 	
# 	x <- lumiR.idat(files, path, manifestfile=manifestfile, columnNameGrepPattern=list(exprs='AVG_Signal'))
# 	expect_that(x, is_a("LumiBatch"))
# })

