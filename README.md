lumidat
=======

An R package for processing Illumina gene expression idat files

Description
===========
This package enables the decryption, and
preprocessing of Illumina gene expression iDAT files
(aka version 1 iDAT files). Previously, the only option
that Illumina gene expression customers had was to rely
upon a microarray core facility with access to an
Illumina Scanner and a copy of Illumina BeadStudio, or
GenomeStudio to pre-process the array data. Our
intention is to allow users to pre-process their
probe-level Illumina data using this software, thereby
enabling them to choose from the collection of
sophisticated normalisation and pre-processing
procedures, which have recently been demonstrated by
Shi et al, 2010 to simultaneously improve noise and
bias, via the *lumi*, or *limma* pipelines. We
have made every effort to reproduce the GenomeStudio
output, down to the Detection P-Value calculation, the
order of the rows, the background correction procedure
(which is only applied to gene-probes), the file formats and names.

Installation
============
    library(devtools)
    install_github("lumidat", "drmjc")

License
=======
Importantly, note that this code is provided under a GenePattern License agreement. See LICENSE file for more info.

Usage
=====
There's an extensive package overview in the package-wide helpfile:

    library(lumidat)
    ?lumidat

Vignette
========
    #
    # import & initial QC
    #
    idat.files <- dir("/path/to/idat/files", full.path=TRUE)
    x.raw <- lumiR.idat(idat.files, probeID="ProbeID", manifestfile="/Volumes/GRIW/ICGCPancreas/icgc_data/icgc_gex/HumanHT-12_V4_0_R2_15002873_B.txt", convertNuID=FALSE,inputAnnotation=FALSE, QC=TRUE, memory="-Xmx4080m")
    dir.create("QC")
    doMA <- doPairs <- FALSE
    plot.lumi.QC.all(x.raw, "QC/01.unnorm/", "raw", "Unnormalised", doMA, doPairs)
    
    #
    # optional, import a probe-level annotation object
    #
    load("Rmisc/IlluminaHT4.anno-2.0.RDa.gz")
    fData(x.raw) <- IlluminaHT4.anno
    fvarMetadata(x.raw) <- structure(list(
    	labelDescription = c(
    		"The Illumina microarray probe identifier", 
    		"Array Address code to identify the probe at the bead-level", 
    		"Lumi's nuID (universal naming scheme for oligos)", 
    		"Quality grade assigned to the probe", 
    		"Coding status of the target sequence: intergenic / intronic / transcriptomic? / transcriptomic", 
    		"50mer probe sequence", 
    		"Genomic coordinates of second best matches between the probe and the genome", 
    		"Genomic coordinates of sequences as alignable with the probe as its main target", 
    		"Overlapping RepeatMasked sequences", 
    		"Overlapping annotated SNPs", 
    		"Entrez Gene ID's after remapping probes", 
    		"Probe's genomic coordinates (hg19, mm9 or rn4)", 
    		"Gene symbol derived by reannotation", 
    		"A more descriptive identifier for controls", 
    		"An identifier for probes designated as controls by Illumina", 
    		"Gene title (derived from EntrezReannotated)",
    		"Probe-level Description"
    	)), 
    	.Names = "labelDescription", 
    	row.names = c(
    		"ProbeID", 
    		"ArrayAddressID", "NuID", "ProbeQuality", "CodingZone", "ProbeSequence", 
    		"SecondMatches", "OtherGenomicMatches", "RepeatMask", "ContainsSNP", 
    		"EntrezReannotated", "GenomicLocation", "SymbolReannotated", 
    		"ReporterGroupName", "ReporterGroupID", "GeneNameReannotated", "Description"
    	),
    	class = "data.frame"
    )
    
    #
    # potentially exclude 'bad' arrays
    #
    x.passedqc <- x.raw
    plot.lumi.QC.all(x.passedqc, "QC/02.passedqc/", "raw", "Unnormalised", doMA, doPairs)
    
    #
    # VST transform
    #
    x.transformed   <- lumiT(x.passedqc, method="vst")
    plot.lumi.QC.all(x.transformed, "QC/03.transformed/", "vst", "Transformed", MA=FALSE, pairs=FALSE)
    
    #
    # normalise. one of rsn, ssn, 
    #
    x.norm <- lumiN(x.transformed, "rsn")
    plot.lumi.QC.all(x.norm, "QC/04.norm/", "rsn", "RSN Normalised", MA=FALSE, pairs=FALSE)
    
    #
    # average replicates
    #
    x.averaged <- average.replicates(x.norm, sub("\\.[12]$", "", sampleNames(x.norm)))
    
    #
    # filter on probe quality:
    #
    Rkeys(illuminaHumanv4PROBEQUALITY)
    # [1] "No match"    "Bad"         "Perfect***"  "Perfect"     "Perfect****"
    # [6] "Good****"    "Good"        "Good***"    
    # keep the good, and perfect probes:
    idx <- grep("Good|Perfect", fData(x.averaged)$ProbeQuality)
    x.averaged.hiqual <- x.averaged[idx,]     
    
    #
    # collapse to 1 row per gene
    # (using 2 of MJC's libraries)
    library(microarrays)
    library(metaGSEA)
    x.averaged.genes <- collapse(x.averaged, T, F, FUN=var, "SymbolReannotated")
    x.averaged.hiqual.genes <- collapse(x.averaged.hiqual, T, F, FUN=var, "SymbolReannotated")
    

    
