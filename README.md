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

Usage
=====
There's an extensive package overview in the package-wide helpfile:

    library(lumidat)
    ?lumidat

    
