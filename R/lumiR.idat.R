#' @include read.illumina.idat.R
roxygen()

#' Read Illumina gene expression iDAT files.
#' This function can decrypt Illumina gene expression iDAT files (aka version 1 iDAT files). It will
#' temporarily create GenomeStudio-compatible output files, and then run \code{\link[lumi]{lumiR}}.
#'
#' See \code{\link{read.illumina.idat}} for more details on Illumina arrays, manifest files, and probeID naming options.
#' See \code{\link[lumi]{lumiR}} for more details on all paramaters from \eqn{detectionTh} onwards.
#'
#' @param files A character vector of at least one file name
#' @param path The path to the directories where the files are. This defaults to the current 
#' working directory.
#' @param probeID This controls which value to identify each probe by. Allowable values are 
#' "ArrayAddressID", "ProbeID", "Sequence", "NuID". See \code{\link{read.illumina.idat}}
#' for more details.
#' @param manifestfile The full path to the Array manifest file in TXT format.
#' @param detectionTh the p-value threshold of determining detectability of the expression.
#'  See more details in \code{\link[lumi]{lumiQ}}.
#' @param na.rm determine whether to remove NA
#' @param parseColumnName determine whether to parse the column names and retrieve the sample 
#' information (Assume the sample information is separated by "\_".)
#' @param checkDupId determine whether to check duplicated TargetIDs or ProbeIds. The 
#' duplicated ones will be averaged.
#' @param QC determine whether to do quality control assessment after read in the data.
#' @param columnNameGrepPattern the string grep patterns used to determine the slot 
#' corresponding columns.
#' @param verbose a boolean to decide whether to print out some messages
#' @param ... other parameters used by \code{\link[utils]{read.table}} function
#'
#' @return return a LumiBatch object
#' @author Mark Cowley, with contributions from Mark Pinese, David Eby.
# @TODO represent manifest files the way that Affymetrix CDF's are (ie in data packages).
# @TODO Import probe-level annotation from the manifest file into the resulting LumiBatch object.
#' @seealso \code{\link{read.illumina.idat}}
#' \code{\link[lumi]{lumiR}}
#' @export
#' @examples
#' library(lumi)
#' path <- system.file("extdata", package="lumidat")
#' files <- c("5356583020_A_Grn.idat", "5356583020_B_Grn.idat")
#' manifestfile <- system.file("extdata", "HumanHT-12_V3_0_R1_99999999.txt", package="lumidat")
#' res <- lumiR.idat(files, path, manifestfile, probeID="NuID")
#'
lumiR.idat <- function(files=NULL, path=NULL, probeID=c("ArrayAddressID", "ProbeID", "Sequence", "NuID"), manifestfile=NULL,
detectionTh=0.01, na.rm=TRUE, parseColumnName=FALSE, checkDupId=TRUE, QC=TRUE, columnNameGrepPattern=list(exprs='AVG_SIGNAL', se.exprs='BEAD_STD', detection='DETECTION', beadNum='Avg_NBEADS'), verbose=TRUE, ...) {
	require(lumi)
	
	outdir <-  tempdir()
	files <- read.illumina.idat(files=files, path=path, probeID=probeID, manifestfile=manifestfile, outdir=outdir, verbose=verbose, collapseMode="none", prefix=NULL, backgroundCorrect=FALSE)
	# files[1] = "<dir>/Sample Probe Profile.txt", files[2] = "<dir>/Control Probe Profile.txt"
	on.exit(unlink(files))	

	if(!file.exists(files[1])) stop("Error occurred during read.illumina.idat.")
	
	# execute the normal lumiR function from the 'lumi' package
	result <- lumiR(files[1], detectionTh=detectionTh, na.rm=na.rm, convertNuID=FALSE, lib.mapping=NULL, dec='.', parseColumnName=parseColumnName, checkDupId=checkDupId, QC=QC, columnNameGrepPattern=columnNameGrepPattern, inputAnnotation=FALSE, annotationColumn="SYMBOL", verbose=verbose, ...)
	
	return( result )
}
