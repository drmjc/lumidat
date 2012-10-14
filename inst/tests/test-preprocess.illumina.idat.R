context("preprocess.illumina.idat suite")

test_that("iDAT files+path as input", {
	path <- system.file("extdata", package="lumidat")
	outdir <- tempdir()
	idat.files <- c("5356583020_A_Grn.idat", "5356583020_B_Grn.idat")
	manifestfile <- system.file("extdata", "HumanHT-12_V3_0_R1_99999999.txt", package="lumidat")
	files <- preprocess.illumina.idat(idat.files, path, manifestfile=manifestfile, probeID="ProbeID", outdir=outdir)

	expect_equal(length(files), 2)
	expect_equal(length(grep("Sample Probe Profile.txt", files)), 1)
	expect_equal(length(grep("Control Probe Profile.txt", files)), 1)
	expect_equal(all(file.exists(files)), TRUE)
	
	df <- read.delim(files[1], stringsAsFactors=FALSE, skip=6)
	res <- df[1:6, ]
	expected.res <- structure(list(TargetID = c("AIFM3", "ARHGDIA", "ARID3A", "BIN1", 
	"BTBD10", "C15ORF39"), ProbeID = c("ILMN_1704173", "ILMN_1734742", 
	"ILMN_1670130", "ILMN_2309245", "ILMN_1809344", "ILMN_1793729"
	), MIN_Signal.5356583020_A = c(88.54737, 1171.861, 1086.613, 
	782.6937, 452.9955, 1373.364), AVG_Signal.5356583020_A = c(88.54737, 
	1171.861, 1086.613, 782.6937, 452.9955, 1373.364), MAX_Signal.5356583020_A = c(88.54737, 
	1171.861, 1086.613, 782.6937, 452.9955, 1373.364), NARRAYS.5356583020_A = c(1L, 
	1L, 1L, 1L, 1L, 1L), ARRAY_STDEV.5356583020_A = c(NaN, NaN, NaN, 
	NaN, NaN, NaN), BEAD_STDEV.5356583020_A = c(36.75683, 209.5618, 
	171.5067, 183.8675, 86.37054, 466.5978), Avg_NBEADS.5356583020_A = c(14L, 
	21L, 15L, 17L, 30L, 15L), Detection.5356583020_A = c(0.9866667, 
	1, 1, 1, 1, 1), MIN_Signal.5356583020_B = c(62.09662, 858.2011, 
	958.4333, 563.4199, 522.0789, 1267.783), AVG_Signal.5356583020_B = c(62.09662, 
	858.2011, 958.4333, 563.4199, 522.0789, 1267.783), MAX_Signal.5356583020_B = c(62.09662, 
	858.2011, 958.4333, 563.4199, 522.0789, 1267.783), NARRAYS.5356583020_B = c(1L, 
	1L, 1L, 1L, 1L, 1L), ARRAY_STDEV.5356583020_B = c(NaN, NaN, NaN, 
	NaN, NaN, NaN), BEAD_STDEV.5356583020_B = c(18.05332, 230.0427, 
	121.8764, 156.202, 130.108, 374.7991), Avg_NBEADS.5356583020_B = c(14L, 
	19L, 21L, 19L, 24L, 20L), Detection.5356583020_B = c(0.8666667, 
	1, 1, 1, 1, 1)), .Names = c("TargetID", "ProbeID", "MIN_Signal.5356583020_A", 
	"AVG_Signal.5356583020_A", "MAX_Signal.5356583020_A", "NARRAYS.5356583020_A", 
	"ARRAY_STDEV.5356583020_A", "BEAD_STDEV.5356583020_A", "Avg_NBEADS.5356583020_A", 
	"Detection.5356583020_A", "MIN_Signal.5356583020_B", "AVG_Signal.5356583020_B", 
	"MAX_Signal.5356583020_B", "NARRAYS.5356583020_B", "ARRAY_STDEV.5356583020_B", 
	"BEAD_STDEV.5356583020_B", "Avg_NBEADS.5356583020_B", "Detection.5356583020_B"
	), row.names = c(NA, 6L), class = "data.frame")
	expect_identical(res, expected.res)

	df <- read.delim(files[2], stringsAsFactors=FALSE, check.names=FALSE)
	res <- df[1:6, ]
	expected.res <- structure(list(TargetID = c("biotin", "biotin", "cy3_hyb", "cy3_hyb", 
	"cy3_hyb", "cy3_hyb"), ProbeID = c("ILMN_1343048", "ILMN_1343049", 
	"ILMN_2038769", "ILMN_1343052", "ILMN_2038771", "ILMN_2038770"
	), `5356583020_A.AVG_Signal` = c(16902.39, 20562.41, 17585.26, 
	410.2783, 3419.561, 15378.33), `5356583020_A.BEAD_STDERR` = c(5062.077, 
	4918.893, 4229.713, 101.1036, 735.4537, 4070.494), `5356583020_A.Avg_NBEADS` = c(85L, 
	67L, 74L, 70L, 64L, 81L), `5356583020_A.Detection Pval` = c(1, 
	1, 1, 1, 1, 1), `5356583020_B.AVG_Signal` = c(16668.27, 20507.38, 
	15846.65, 407.3946, 3164.278, 15831.14), `5356583020_B.BEAD_STDERR` = c(6532.465, 
	5919.903, 3202.289, 115.9396, 647.5594, 4153.448), `5356583020_B.Avg_NBEADS` = c(60L, 
	84L, 62L, 51L, 77L, 63L), `5356583020_B.Detection Pval` = c(1, 
	1, 1, 1, 1, 1)), .Names = c("TargetID", "ProbeID", "5356583020_A.AVG_Signal", 
	"5356583020_A.BEAD_STDERR", "5356583020_A.Avg_NBEADS", "5356583020_A.Detection Pval", 
	"5356583020_B.AVG_Signal", "5356583020_B.BEAD_STDERR", "5356583020_B.Avg_NBEADS", 
	"5356583020_B.Detection Pval"), row.names = c(NA, 6L), class = "data.frame")
	expect_identical(res, expected.res)
})

test_that("zipfile as input", {
	path <- system.file("extdata", package="lumidat")
	outdir <- tempdir()
	idat.files <- c("5356583020_A_Grn.idat", "5356583020_B_Grn.idat")
	zipfile <- tempfile(fileext=".zip")
	zip(zipfile, file.path(path, idat.files), flags="-r9Xq")
	manifestfile <- system.file("extdata", "HumanHT-12_V3_0_R1_99999999.txt", package="lumidat")
	files <- preprocess.illumina.idat(zipfile=zipfile, manifestfile=manifestfile, probeID="ProbeID", outdir=outdir)

	expect_equal(length(files), 2)
	expect_equal(length(grep("Sample Probe Profile.txt", files)), 1)
	expect_equal(length(grep("Control Probe Profile.txt", files)), 1)
	expect_equal(all(file.exists(files)), TRUE)
	
	df <- read.delim(files[1], stringsAsFactors=FALSE, skip=6)
	res <- df[1:6, ]
	expected.res <- structure(list(TargetID = c("AIFM3", "ARHGDIA", "ARID3A", "BIN1", 
	"BTBD10", "C15ORF39"), ProbeID = c("ILMN_1704173", "ILMN_1734742", 
	"ILMN_1670130", "ILMN_2309245", "ILMN_1809344", "ILMN_1793729"
	), MIN_Signal.5356583020_A = c(88.54737, 1171.861, 1086.613, 
	782.6937, 452.9955, 1373.364), AVG_Signal.5356583020_A = c(88.54737, 
	1171.861, 1086.613, 782.6937, 452.9955, 1373.364), MAX_Signal.5356583020_A = c(88.54737, 
	1171.861, 1086.613, 782.6937, 452.9955, 1373.364), NARRAYS.5356583020_A = c(1L, 
	1L, 1L, 1L, 1L, 1L), ARRAY_STDEV.5356583020_A = c(NaN, NaN, NaN, 
	NaN, NaN, NaN), BEAD_STDEV.5356583020_A = c(36.75683, 209.5618, 
	171.5067, 183.8675, 86.37054, 466.5978), Avg_NBEADS.5356583020_A = c(14L, 
	21L, 15L, 17L, 30L, 15L), Detection.5356583020_A = c(0.9866667, 
	1, 1, 1, 1, 1), MIN_Signal.5356583020_B = c(62.09662, 858.2011, 
	958.4333, 563.4199, 522.0789, 1267.783), AVG_Signal.5356583020_B = c(62.09662, 
	858.2011, 958.4333, 563.4199, 522.0789, 1267.783), MAX_Signal.5356583020_B = c(62.09662, 
	858.2011, 958.4333, 563.4199, 522.0789, 1267.783), NARRAYS.5356583020_B = c(1L, 
	1L, 1L, 1L, 1L, 1L), ARRAY_STDEV.5356583020_B = c(NaN, NaN, NaN, 
	NaN, NaN, NaN), BEAD_STDEV.5356583020_B = c(18.05332, 230.0427, 
	121.8764, 156.202, 130.108, 374.7991), Avg_NBEADS.5356583020_B = c(14L, 
	19L, 21L, 19L, 24L, 20L), Detection.5356583020_B = c(0.8666667, 
	1, 1, 1, 1, 1)), .Names = c("TargetID", "ProbeID", "MIN_Signal.5356583020_A", 
	"AVG_Signal.5356583020_A", "MAX_Signal.5356583020_A", "NARRAYS.5356583020_A", 
	"ARRAY_STDEV.5356583020_A", "BEAD_STDEV.5356583020_A", "Avg_NBEADS.5356583020_A", 
	"Detection.5356583020_A", "MIN_Signal.5356583020_B", "AVG_Signal.5356583020_B", 
	"MAX_Signal.5356583020_B", "NARRAYS.5356583020_B", "ARRAY_STDEV.5356583020_B", 
	"BEAD_STDEV.5356583020_B", "Avg_NBEADS.5356583020_B", "Detection.5356583020_B"
	), row.names = c(NA, 6L), class = "data.frame")
	expect_identical(res, expected.res)

	df <- read.delim(files[2], stringsAsFactors=FALSE, check.names=FALSE)
	res <- df[1:6, ]
	expected.res <- structure(list(TargetID = c("biotin", "biotin", "cy3_hyb", "cy3_hyb", 
	"cy3_hyb", "cy3_hyb"), ProbeID = c("ILMN_1343048", "ILMN_1343049", 
	"ILMN_2038769", "ILMN_1343052", "ILMN_2038771", "ILMN_2038770"
	), `5356583020_A.AVG_Signal` = c(16902.39, 20562.41, 17585.26, 
	410.2783, 3419.561, 15378.33), `5356583020_A.BEAD_STDERR` = c(5062.077, 
	4918.893, 4229.713, 101.1036, 735.4537, 4070.494), `5356583020_A.Avg_NBEADS` = c(85L, 
	67L, 74L, 70L, 64L, 81L), `5356583020_A.Detection Pval` = c(1, 
	1, 1, 1, 1, 1), `5356583020_B.AVG_Signal` = c(16668.27, 20507.38, 
	15846.65, 407.3946, 3164.278, 15831.14), `5356583020_B.BEAD_STDERR` = c(6532.465, 
	5919.903, 3202.289, 115.9396, 647.5594, 4153.448), `5356583020_B.Avg_NBEADS` = c(60L, 
	84L, 62L, 51L, 77L, 63L), `5356583020_B.Detection Pval` = c(1, 
	1, 1, 1, 1, 1)), .Names = c("TargetID", "ProbeID", "5356583020_A.AVG_Signal", 
	"5356583020_A.BEAD_STDERR", "5356583020_A.Avg_NBEADS", "5356583020_A.Detection Pval", 
	"5356583020_B.AVG_Signal", "5356583020_B.BEAD_STDERR", "5356583020_B.Avg_NBEADS", 
	"5356583020_B.Detection Pval"), row.names = c(NA, 6L), class = "data.frame")
	expect_identical(res, expected.res)
})


test_that("specifying CLM file", {
	path <- system.file("extdata", package="lumidat")
	outdir <- tempdir()
	idat.files <- c("5356583020_A_Grn.idat", "5356583020_B_Grn.idat")
	manifestfile <- system.file("extdata", "HumanHT-12_V3_0_R1_99999999.txt", package="lumidat")
	clmfile <- system.file("extdata", "5356583020.clm", package="lumidat")
	files <- preprocess.illumina.idat(idat.files, path, manifestfile=manifestfile, clmfile=clmfile, probeID="ProbeID", outdir=outdir)

	expect_equal(length(files), 2)
	expect_equal(length(grep("Sample Probe Profile.txt", files)), 1)
	expect_equal(length(grep("Control Probe Profile.txt", files)), 1)
	expect_equal(all(file.exists(files)), TRUE)
	
	################################################################################
	df <- read.delim(files[1], stringsAsFactors=FALSE, skip=6)
	res <- df[1:6, ]
	expected.res <- structure(list(TargetID = c("AIFM3", "ARHGDIA", "ARID3A", "BIN1", 
	"BTBD10", "C15ORF39"), ProbeID = c("ILMN_1704173", "ILMN_1734742", 
	"ILMN_1670130", "ILMN_2309245", "ILMN_1809344", "ILMN_1793729"
	), MIN_Signal.A = c(88.54737, 1171.861, 1086.613, 782.6937, 452.9955, 
	1373.364), AVG_Signal.A = c(88.54737, 1171.861, 1086.613, 782.6937, 
	452.9955, 1373.364), MAX_Signal.A = c(88.54737, 1171.861, 1086.613, 
	782.6937, 452.9955, 1373.364), NARRAYS.A = c(1L, 1L, 1L, 1L, 
	1L, 1L), ARRAY_STDEV.A = c(NaN, NaN, NaN, NaN, NaN, NaN), BEAD_STDEV.A = c(36.75683, 
	209.5618, 171.5067, 183.8675, 86.37054, 466.5978), Avg_NBEADS.A = c(14L, 
	21L, 15L, 17L, 30L, 15L), Detection.A = c(0.9866667, 1, 1, 1, 
	1, 1), MIN_Signal.B = c(62.09662, 858.2011, 958.4333, 563.4199, 
	522.0789, 1267.783), AVG_Signal.B = c(62.09662, 858.2011, 958.4333, 
	563.4199, 522.0789, 1267.783), MAX_Signal.B = c(62.09662, 858.2011, 
	958.4333, 563.4199, 522.0789, 1267.783), NARRAYS.B = c(1L, 1L, 
	1L, 1L, 1L, 1L), ARRAY_STDEV.B = c(NaN, NaN, NaN, NaN, NaN, NaN
	), BEAD_STDEV.B = c(18.05332, 230.0427, 121.8764, 156.202, 130.108, 
	374.7991), Avg_NBEADS.B = c(14L, 19L, 21L, 19L, 24L, 20L), Detection.B = c(0.8666667, 
	1, 1, 1, 1, 1)), .Names = c("TargetID", "ProbeID", "MIN_Signal.A", 
	"AVG_Signal.A", "MAX_Signal.A", "NARRAYS.A", "ARRAY_STDEV.A", 
	"BEAD_STDEV.A", "Avg_NBEADS.A", "Detection.A", "MIN_Signal.B", 
	"AVG_Signal.B", "MAX_Signal.B", "NARRAYS.B", "ARRAY_STDEV.B", 
	"BEAD_STDEV.B", "Avg_NBEADS.B", "Detection.B"), row.names = c(NA, 
	6L), class = "data.frame")
	expect_identical(res, expected.res)
	
	################################################################################
	df <- read.delim(files[2], stringsAsFactors=FALSE, check.names=FALSE)
	res <- df[1:6, ]
	expected.res <- structure(list(TargetID = c("biotin", "biotin", "cy3_hyb", "cy3_hyb", 
	"cy3_hyb", "cy3_hyb"), ProbeID = c("ILMN_1343048", "ILMN_1343049", 
	"ILMN_2038769", "ILMN_1343052", "ILMN_2038771", "ILMN_2038770"
	), A.AVG_Signal = c(16902.39, 20562.41, 17585.26, 410.2783, 3419.561, 
	15378.33), A.BEAD_STDERR = c(5062.077, 4918.893, 4229.713, 101.1036, 
	735.4537, 4070.494), A.Avg_NBEADS = c(85L, 67L, 74L, 70L, 64L, 
	81L), `A.Detection Pval` = c(1, 1, 1, 1, 1, 1), B.AVG_Signal = c(16668.27, 
	20507.38, 15846.65, 407.3946, 3164.278, 15831.14), B.BEAD_STDERR = c(6532.465, 
	5919.903, 3202.289, 115.9396, 647.5594, 4153.448), B.Avg_NBEADS = c(60L, 
	84L, 62L, 51L, 77L, 63L), `B.Detection Pval` = c(1, 1, 1, 1, 
	1, 1)), .Names = c("TargetID", "ProbeID", "A.AVG_Signal", "A.BEAD_STDERR", 
	"A.Avg_NBEADS", "A.Detection Pval", "B.AVG_Signal", "B.BEAD_STDERR", 
	"B.Avg_NBEADS", "B.Detection Pval"), row.names = c(NA, 6L), class = "data.frame")
	
	expect_identical(res, expected.res)
	
})


# probeID=c("ProbeID", "ArrayAddressID", "Sequence", "NuID"), 

test_that("probeID=ArrayAddressID works", {
	path <- system.file("extdata", package="lumidat")
	outdir <- tempdir()
	idat.files <- c("5356583020_A_Grn.idat", "5356583020_B_Grn.idat")
	manifestfile <- system.file("extdata", "HumanHT-12_V3_0_R1_99999999.txt", package="lumidat")
	files <- preprocess.illumina.idat(idat.files, path, manifestfile=manifestfile, probeID="ArrayAddressID", outdir=outdir)

	expect_equal(length(files), 2)
	expect_equal(length(grep("Sample Probe Profile.txt", files)), 1)
	expect_equal(length(grep("Control Probe Profile.txt", files)), 1)
	expect_equal(all(file.exists(files)), TRUE)
	
	################################################################################
	df <- read.delim(files[1], stringsAsFactors=FALSE, skip=6)
	res <- df$ProbeID
	expected.res <- c(4260551L, 4760255L, 6510553L, 650553L, 4250544L, 870110L, 2100176L, 
	1740731L, 6420520L, 6020544L, 7320187L, 6020725L, 2480717L, 5080746L, 
	2710020L, 7570154L, 2480672L, 1850044L, 1340112L, 4780731L, 1470398L, 
	6400477L, 5130397L, 3710372L, 3780097L, 4850010L, 6520037L, 2470450L, 
	4810470L, 5960372L, 290020L, 2260747L, 2360681L, 7400692L, 3130152L, 
	7320189L, 4780195L, 2570411L, 5900364L, 3360129L, 7160367L, 6650692L, 
	1570184L, 5570068L, 1710221L, 6110274L, 3870215L, 6180280L, 5860608L, 
	3360368L, 1030609L, 6960706L, 4150010L, 2640750L, 6180112L, 3310608L, 
	4480719L, 2570463L, 540717L, 6940338L, 5960196L, 4760674L, 1430450L, 
	2650040L, 3310056L, 2480719L, 6660333L, 4180259L, 7510243L, 6650189L, 
	520431L, 4560520L, 650102L, 4230026L, 4920193L, 4730168L, 830441L, 
	1400546L, 290685L, 6270445L, 50475L, 1500307L, 7150059L, 5340180L, 
	4230136L, 4900209L, 1300239L, 2630356L, 1190064L, 5420450L)
	expect_identical(res, expected.res)
	
})


test_that("probeID=Sequence works", {
	path <- system.file("extdata", package="lumidat")
	outdir <- tempdir()
	idat.files <- c("5356583020_A_Grn.idat", "5356583020_B_Grn.idat")
	manifestfile <- system.file("extdata", "HumanHT-12_V3_0_R1_99999999.txt", package="lumidat")
	files <- preprocess.illumina.idat(idat.files, path, manifestfile=manifestfile, probeID="Sequence", outdir=outdir)

	expect_equal(length(files), 2)
	expect_equal(length(grep("Sample Probe Profile.txt", files)), 1)
	expect_equal(length(grep("Control Probe Profile.txt", files)), 1)
	expect_equal(all(file.exists(files)), TRUE)
	
	################################################################################
	df <- read.delim(files[1], stringsAsFactors=FALSE, skip=6)
	res <- df$ProbeID[1:10]
	expected.res <- c("AGCTCACATGCAGTAGACTTGGGCAGGCAAAGGGGGCACCAAGGGCACAG", "CAGGCCTTCCCTGACCCAGCCAGGAACAAACAAGGGACCAAGTGCACACA", 
	"CCACACACTCACCACTCCCAGCTTCTCGTGTCCAGTGAAACCCCTGAACC", "AGGGGCGTTCTCCCAAAGATTAGGTCGTTTTCCAAAGAGCCGCGTCCCGG", 
	"TCACACAGCTAGGTCTTTTCAGAAGTGGTGGAAATTGGCAGCTGGGGTAC", "CTTGCCTAGAGAACACACATGGGCTTTGGAGCCCGACAGACCTGGGCTTG", 
	"GCTGAGAGAGAAAGAAGCAGCTCTTGAAGAAATGCGTAAGAAGATGCACC", "CCAGGCTGGTCTCAGAATTCTGGTGAGTGATCCTCCCACAGTGGACTTCC", 
	"CATGGATGCCAACCGGTCACCCAGGAGGATGGCAAAGAGAGTCGCATCTC", "GCCCACAAGGAGGACTTACAATGTGGAACTTAGTCTCTTTCCCTCACTCC"
	)
	expect_identical(res, expected.res)
})


test_that("probeID=NuID works", {
	path <- system.file("extdata", package="lumidat")
	outdir <- tempdir()
	idat.files <- c("5356583020_A_Grn.idat", "5356583020_B_Grn.idat")
	manifestfile <- system.file("extdata", "HumanHT-12_V3_0_R1_99999999.txt", package="lumidat")
	files <- preprocess.illumina.idat(idat.files, path, manifestfile=manifestfile, probeID="NuID", outdir=outdir)

	expect_equal(length(files), 2)
	expect_equal(length(grep("Sample Probe Profile.txt", files)), 1)
	expect_equal(length(grep("Control Probe Profile.txt", files)), 1)
	expect_equal(all(file.exists(files)), TRUE)
	
	################################################################################
	df <- read.delim(files[1], stringsAsFactors=FALSE, skip=6)
	res <- df$ProbeID
	expected.res <- c("lJ0Tksh_pKQKqRQqRI", "NSl9XhUlKBAQqFC5EQ", "iURHRR1Sfdu1LgFXgU", 
	"QKpvdUCPK2.1AiWbVo", "Q0RJyt.0guugPpJ6rE", "0flyIEROp.olYSF6n4", 
	"ZniIgIJJ34IDmwgjkU", "9Up63SD3ri411RLofU", "QTo5QWtFSijpAiLZN0", 
	"NlRCih8Q7oHy3f1dHU", "0UJD25LqOJx.lJyKWU", "rtCnUep15THUpc_0e4", 
	"KyqQynMZxJcruyylEU", "WjiN7uk.cmd3A87d6Q", "EuUnlPkeXRP9fyO.iQ", 
	"0CtCopYKLu31Un__KU", "TSCUT9C74ElEbHLpC4", "K30srb.UieSsA59VHk", 
	"B59AqJL9z9a5UaI3Yo", "cfJqib6pX1.zKrN3lc", "ow0oVJA.VYJDqtWup0", 
	"c6CC_6g8ZY6VgNCkik", "cIsqOaCSYutKLublAI", "xpTSCGKSc_i8ZJSJQo", 
	"uAtIy0uIFJyDU43vuo", "HSKdOrfFOrjQoKnrB0", "NlIi4iu8.S8qy0qCXc", 
	"lFkiX_F2sXQzgj0q30", "6feJMxuxv7lTlMg.iU", "ooSh0qCxhyHoDgZqdU", 
	"fegCQD_j_69DUU38dI", "rgIj.rKIqeKeEmqI5Q", "f9eA.E.jOSMf0nVgOE", 
	"fgteSVG3cKSR7UMEt0", "Q3VLorv0mjnmt8V_yU", "uponqPQOQk96dSiEZ0", 
	"TeqCiiqCUtSKeUT3.s", "E_qksA7Eh1qgohSgoo", "NEX0oqCV8.er4HVfU4", 
	"EJxLDQh7u8_sioyETI", "igl7ntl1PCsUjs.p4o", "HnnVoDkDiHXoB6uoPo", 
	"f8rl3H6iagiiOEoJhk", "xhBTlRSFedTRXDlFRA", "rp13_p1x6D80lNLk3c", 
	"xUWKuCKQulcJAIblLs", "QtVBXBWhekTEIT0kjo", "QBeq3pMkJNKJiVdDeI", 
	"3KvUoTLeKQjoqziqX0", "i6n5BujgiOCAPTQvNY", "0gI6ol2ty6OYERXx54", 
	"rIobfPN67en54.DBcg", "rhNHa3uH2uQ3X0qCWo", "Q8gegTQvuHg_eS3TAE", 
	"ZVWqYlqkmeR2uZFpJs", "HN1aK.CnSSjSeteN_0", "9ckqJrioiaej9_ajeQ", 
	"6jh6c.oAXq5x5Qers0", "l9NFipeiuR9ugRxaiM", "ueuFSIRPPhXdcl1P10", 
	"9gVhUokj4ovE0hczoo", "uRDOULMF6.byK7cS1I", "9Vu9RtOo91tFqKt_Bo", 
	"Z_KKCqEvkt0fqEpRJQ", "E9VUedIhpiBzvufeh0", "c5KhFJJcXXFUnt74iQ", 
	"BFOlu73Rep17Zdefug", "Z7q7I6r2J.XEhc90os", "BOkEIILakIKIILBOik", 
	"fqiMRwVUgOHnCXfl_0", "iJUrcDsOvr8_9zBVJU", "rrDGQSF.KPTVCUAOqg", 
	"3ZfhT0oSUT4jUSnr3U", "uvcl_OX6GWEuCqCpPc", "rTUh4oIVoB8qlG4lGk", 
	"oUuoJKG3qfxr3KhKII", "N3KrO.Rfj.pWvf0u9Y", "xUBB8BSTR.7fQ.l9S4", 
	"KS7gsL3TUE.delNT3U", "i4U3n76eXq06.n6c40", "upWX307QohTkj44igk", 
	"3kE_IgB5VTdCQg6b68", "97CRbQTSEXHToTuvlo", "KSXe7i6OPSvlIhIFdI", 
	"0lQpzJV19DPhcGqiAI", "uklAkodfpQUYyILXjo", "xXl7eXuF7sbPEp.KFI", 
	"Zh7iRpOlJtN1Kfh9Do", "WwEX2r5JUYXRevUsUc", "WicBF9q_SVGF0Xr1LE")
	
	expect_identical(res, expected.res)
})

################################################################################

# collapseMode=c("none", "max", "median", "mean"), 

test_that("collapseMode=none works", {
	path <- system.file("extdata", package="lumidat")
	outdir <- tempdir()
	idat.files <- c("5356583020_A_Grn.idat", "5356583020_B_Grn.idat")
	manifestfile <- system.file("extdata", "HumanHT-12_V3_0_R1_99999999.txt", package="lumidat")
	files <- preprocess.illumina.idat(idat.files, path, manifestfile=manifestfile, probeID="Probe", collapseMode="none", outdir=outdir)

	expect_equal(length(files), 2)
	expect_equal(length(grep("Sample Probe Profile.txt", files)), 1)
	expect_equal(length(grep("Control Probe Profile.txt", files)), 1)
	expect_equal(all(file.exists(files)), TRUE)
	
	################################################################################
	df <- read.delim(files[1], stringsAsFactors=FALSE, skip=6)
	expect_equal(nrow(df), 90)
	
	res <- df$AVG_Signal.5356583020_A
	expected.res <- c(88.54737, 1171.861, 1086.613, 782.6937, 452.9955, 1373.364, 
	45.57533, 57.69991, 97.63488, 227.8882, 1042.684, 57.50295, 114.5324, 
	265.9068, 88.38239, 56.64154, 82.54376, 48.48536, 61.01905, 78.97118, 
	50.4348, 60.58284, 83.12913, 59.57077, 55.27691, 41.89008, 60.99137, 
	54.97004, 59.08458, 62.21026, 56.95772, 64.65753, 47.86235, 56.09657, 
	72.27145, 70.67132, 59.31025, 67.65513, 55.99724, 57.69477, 58.80888, 
	157.8183, 52.4112, 69.06544, 67.61997, 51.81778, 58.11007, 58.27464, 
	4004.936, 52.5428, 55.66942, 59.79181, 68.53111, 46.03054, 63.22699, 
	53.61557, 73.92967, 503.2489, 282.5659, 62.29892, 53.42324, 166.1604, 
	75.38274, 381.8736, 402.4396, 247.954, 55.99047, 67.22208, 58.62053, 
	678.9595, 4426.523, 64.53059, 112.9638, 57.53606, 12576.61, 205.1812, 
	187.1925, 142.3655, 76.1413, 58.84699, 659.4638, 309.2188, 564.6039, 
	94.59225, 110.703, 321.0327, 709.85, 438.7442, 58.9047, 62.35626
	)
	expect_equivalent(res, expected.res)
	
})

# TODO, there's 90 rows in df, but 82 unique TargetID's. ie just like collapseMode=none;
# surely there should be just 82 rows in res ??
# test_that("collapseMode=max works", {
# 	path <- system.file("extdata", package="lumidat")
# 	outdir <- tempdir()
# 	idat.files <- c("5356583020_A_Grn.idat", "5356583020_B_Grn.idat")
# 	manifestfile <- system.file("extdata", "HumanHT-12_V3_0_R1_99999999.txt", package="lumidat")
# 	files <- preprocess.illumina.idat(idat.files, path, manifestfile=manifestfile, probeID="Probe", collapseMode="max", outdir=outdir)
# 
# 	expect_equal(length(files), 2)
# 	expect_equal(length(grep("Sample Probe Profile.txt", files)), 1)
# 	expect_equal(length(grep("Control Probe Profile.txt", files)), 1)
# 	expect_equal(all(file.exists(files)), TRUE)
# 	
# 	################################################################################
# 	df <- read.delim(files[1], stringsAsFactors=FALSE, skip=6)
# 	expect_equal(nrow(df), 90)
# 	
# 	res <- df$AVG_Signal.5356583020_A
# 	expected.res <- c(88.54737, 1171.861, 1086.613, 782.6937, 452.9955, 1373.364, 
# 	45.57533, 57.69991, 97.63488, 227.8882, 1042.684, 57.50295, 114.5324, 
# 	265.9068, 88.38239, 56.64154, 82.54376, 48.48536, 61.01905, 78.97118, 
# 	50.4348, 60.58284, 83.12913, 59.57077, 55.27691, 41.89008, 60.99137, 
# 	54.97004, 59.08458, 62.21026, 56.95772, 64.65753, 47.86235, 56.09657, 
# 	72.27145, 70.67132, 59.31025, 67.65513, 55.99724, 57.69477, 58.80888, 
# 	157.8183, 52.4112, 69.06544, 67.61997, 51.81778, 58.11007, 58.27464, 
# 	4004.936, 52.5428, 55.66942, 59.79181, 68.53111, 46.03054, 63.22699, 
# 	53.61557, 73.92967, 503.2489, 282.5659, 62.29892, 53.42324, 166.1604, 
# 	75.38274, 381.8736, 402.4396, 247.954, 55.99047, 67.22208, 58.62053, 
# 	678.9595, 4426.523, 64.53059, 112.9638, 57.53606, 12576.61, 205.1812, 
# 	187.1925, 142.3655, 76.1413, 58.84699, 659.4638, 309.2188, 564.6039, 
# 	94.59225, 110.703, 321.0327, 709.85, 438.7442, 58.9047, 62.35626
# 	)
# 	expect_equivalent(res, expected.res)
# 	
# })
# 

################################################################################
################################################################################

test_that("backgroundCorrect=FALSE", {
	path <- system.file("extdata", package="lumidat")
	outdir <- tempdir()
	idat.files <- c("5356583020_A_Grn.idat", "5356583020_B_Grn.idat")
	manifestfile <- system.file("extdata", "HumanHT-12_V3_0_R1_99999999.txt", package="lumidat")
	files <- preprocess.illumina.idat(idat.files, path, manifestfile=manifestfile, probeID="ProbeID", backgroundCorrect=FALSE, outdir=outdir)

	expect_equal(length(files), 2)
	expect_equal(length(grep("Sample Probe Profile.txt", files)), 1)
	expect_equal(length(grep("Control Probe Profile.txt", files)), 1)
	expect_equal(all(file.exists(files)), TRUE)
	
	df <- read.delim(files[1], stringsAsFactors=FALSE, skip=6)
	res <- df$AVG_Signal.5356583020_A
	expected.res <- c(88.54737, 1171.861, 1086.613, 782.6937, 452.9955, 1373.364, 
	45.57533, 57.69991, 97.63488, 227.8882, 1042.684, 57.50295, 114.5324, 
	265.9068, 88.38239, 56.64154, 82.54376, 48.48536, 61.01905, 78.97118, 
	50.4348, 60.58284, 83.12913, 59.57077, 55.27691, 41.89008, 60.99137, 
	54.97004, 59.08458, 62.21026, 56.95772, 64.65753, 47.86235, 56.09657, 
	72.27145, 70.67132, 59.31025, 67.65513, 55.99724, 57.69477, 58.80888, 
	157.8183, 52.4112, 69.06544, 67.61997, 51.81778, 58.11007, 58.27464, 
	4004.936, 52.5428, 55.66942, 59.79181, 68.53111, 46.03054, 63.22699, 
	53.61557, 73.92967, 503.2489, 282.5659, 62.29892, 53.42324, 166.1604, 
	75.38274, 381.8736, 402.4396, 247.954, 55.99047, 67.22208, 58.62053, 
	678.9595, 4426.523, 64.53059, 112.9638, 57.53606, 12576.61, 205.1812, 
	187.1925, 142.3655, 76.1413, 58.84699, 659.4638, 309.2188, 564.6039, 
	94.59225, 110.703, 321.0327, 709.85, 438.7442, 58.9047, 62.35626)
	expect_equivalent(res, expected.res)
	
	# controls shouldn't change:
	df <- read.delim(files[2], stringsAsFactors=FALSE, check.names=FALSE)
	res <- df[,"5356583020_A.AVG_Signal"]
	expected.res <- c(16902.39, 20562.41, 17585.26, 410.2783, 3419.561, 15378.33, 
	2753.392, 410.9207, 25780.96, 34023.34, 610.7968, 4953.948, 2697.165, 
	26310.26, 7641.937, 55.80587, 50.49826, 5045.01, 3166.717, 17585.26, 
	92.17317, 3419.561, 226.2774, 15378.33, 2753.392, 54.64553, 73.54171, 
	50.45399, 47.25956, 50.85957, 58.00498, 51.28456, 58.34096, 56.62717, 
	61.71053, 51.18339, 46.08123, 57.51138, 69.44334, 52.17262, 65.5533, 
	53.98127, 55.92076, 54.97611, 46.28537, 54.5341, 56.16732, 48.66447, 
	69.75181, 61.06767, 62.44195, 63.55336, 57.39806, 54.35235, 60.54974, 
	111.3958, 53.59858, 60.24976, 50.15711, 58.53523, 71.07463, 49.90788, 
	77.31194, 60.30117, 54.28918, 62.43774, 50.06254, 60.85621, 56.52922, 
	59.71457, 64.6293, 52.41349, 57.03858, 67.88881, 63.83892, 49.25184, 
	55.90054, 71.01463, 52.44993, 55.49572, 58.83201, 58.76225, 62.23101, 
	50.59364, 59.4303, 59.00262, 61.29153, 54.15204, 59.29384, 57.19605, 
	60.40684, 55.44492, 59.80412, 56.6338, 48.49327, 59.32286, 46.62666, 
	41.80057, 59.03254, 52.70005)
	expect_equivalent(res, expected.res)
	
})

test_that("backgroundCorrect=TRUE", {
	path <- system.file("extdata", package="lumidat")
	outdir <- tempdir()
	idat.files <- c("5356583020_A_Grn.idat", "5356583020_B_Grn.idat")
	manifestfile <- system.file("extdata", "HumanHT-12_V3_0_R1_99999999.txt", package="lumidat")
	files <- preprocess.illumina.idat(idat.files, path, manifestfile=manifestfile, probeID="ProbeID", backgroundCorrect=TRUE, outdir=outdir)

	expect_equal(length(files), 2)
	expect_equal(length(grep("Sample Probe Profile.txt", files)), 1)
	expect_equal(length(grep("Control Probe Profile.txt", files)), 1)
	expect_equal(all(file.exists(files)), TRUE)
	
	df <- read.delim(files[1], stringsAsFactors=FALSE, skip=6)
	res <- df$AVG_Signal.5356583020_A
	expected.res <- c(30.52452, 1113.839, 1028.59, 724.6708, 394.9726, 1315.341, 
	-12.4475, -0.32294, 39.61203, 169.8654, 984.6613, -0.5199, 56.50955, 
	207.884, 30.35954, -1.38131, 24.52092, -9.53748, 2.996201, 20.94833, 
	-7.58804, 2.55999, 25.10628, 1.547928, -2.74593, -16.1328, 2.968529, 
	-3.0528, 1.061729, 4.187412, -1.06513, 6.634686, -10.1605, -1.92628, 
	14.24861, 12.64847, 1.287407, 9.632282, -2.02561, -0.32808, 0.786034, 
	99.7955, -5.61165, 11.0426, 9.597126, -6.20507, 0.087227, 0.251797, 
	3946.913, -5.48005, -2.35342, 1.768963, 10.50827, -11.9923, 5.204147, 
	-4.40728, 15.90683, 445.2261, 224.543, 4.276077, -4.59961, 108.1375, 
	17.35989, 323.8508, 344.4168, 189.9312, -2.03238, 9.199238, 0.597679, 
	620.9367, 4368.5, 6.50774, 54.94099, -0.48678, 12518.59, 147.1583, 
	129.1697, 84.34264, 18.11846, 0.824142, 601.441, 251.1959, 506.5811, 
	36.56941, 52.6802, 263.0098, 651.8271, 380.7214, 0.881851, 4.333416)
	expect_equivalent(res, expected.res)

	# controls shouldn't change:
	df <- read.delim(files[2], stringsAsFactors=FALSE, check.names=FALSE)
	res <- df[,"5356583020_A.AVG_Signal"]
	expected.res <- c(16902.39, 20562.41, 17585.26, 410.2783, 3419.561, 15378.33, 
	2753.392, 410.9207, 25780.96, 34023.34, 610.7968, 4953.948, 2697.165, 
	26310.26, 7641.937, 55.80587, 50.49826, 5045.01, 3166.717, 17585.26, 
	92.17317, 3419.561, 226.2774, 15378.33, 2753.392, 54.64553, 73.54171, 
	50.45399, 47.25956, 50.85957, 58.00498, 51.28456, 58.34096, 56.62717, 
	61.71053, 51.18339, 46.08123, 57.51138, 69.44334, 52.17262, 65.5533, 
	53.98127, 55.92076, 54.97611, 46.28537, 54.5341, 56.16732, 48.66447, 
	69.75181, 61.06767, 62.44195, 63.55336, 57.39806, 54.35235, 60.54974, 
	111.3958, 53.59858, 60.24976, 50.15711, 58.53523, 71.07463, 49.90788, 
	77.31194, 60.30117, 54.28918, 62.43774, 50.06254, 60.85621, 56.52922, 
	59.71457, 64.6293, 52.41349, 57.03858, 67.88881, 63.83892, 49.25184, 
	55.90054, 71.01463, 52.44993, 55.49572, 58.83201, 58.76225, 62.23101, 
	50.59364, 59.4303, 59.00262, 61.29153, 54.15204, 59.29384, 57.19605, 
	60.40684, 55.44492, 59.80412, 56.6338, 48.49327, 59.32286, 46.62666, 
	41.80057, 59.03254, 52.70005)
	expect_equivalent(res, expected.res)
	
})


test_that("prefix works", {
	path <- system.file("extdata", package="lumidat")
	outdir <- tempdir()
	idat.files <- c("5356583020_A_Grn.idat", "5356583020_B_Grn.idat")
	manifestfile <- system.file("extdata", "HumanHT-12_V3_0_R1_99999999.txt", package="lumidat")
	files <- preprocess.illumina.idat(idat.files, path, manifestfile=manifestfile, probeID="ProbeID", prefix="MYPREFIX", outdir=outdir)

	expect_equal(length(files), 2)
	expect_equal(length(grep("Sample Probe Profile.txt", files)), 1)
	expect_equal(length(grep("Control Probe Profile.txt", files)), 1)
	expect_equal(all(file.exists(files)), TRUE)
	expect_equal(files[1], file.path(outdir, "MYPREFIX Sample Probe Profile.txt"))
	expect_equal(files[2], file.path(outdir, "MYPREFIX Control Probe Profile.txt"))
	expect_equal(all(file.exists(files)), TRUE)
})
