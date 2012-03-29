#' @include preprocess.illumina.idat.R
NULL

#' Read Illumina gene expression iDAT files.
#' 
#' This function can decrypt Illumina gene expression iDAT files (aka version 1 iDAT files). It will
#' temporarily create GenomeStudio-compatible output files, and then run \code{\link[lumi]{lumiR}}.
#'
#' See \code{\link{preprocess.illumina.idat}} for more details on Illumina arrays, manifest files, and probeID naming options.
#' See \code{\link[lumi]{lumiR}} for more details on all paramaters from \code{detectionTh} onwards.
#' 
#' @inheritParams preprocess.illumina.idat
#' @param detectionTh the p-value threshold of determining detectability of the expression.
#'  See more details in \code{\link[lumi]{lumiQ}}.
#' @param na.rm logical: remove \code{NA}?
#' @param parseColumnName logical: parse the column names and retrieve the sample 
#' information? (Assume the sample information is separated by \dQuote{\_}.)
#' @param checkDupId logical: check duplicated TargetIDs or ProbeIds? The 
#' duplicated ones will be averaged.
#' @param QC logical:do quality control assessment after read in the data?
#' @param columnNameGrepPattern list of named character(1)'s: the grep patterns used to 
#' determine which column name corresponds to which slot. Eg the column named \dQuote{AVG_SIGNAL}
#' will be put into the \code{exprs} slot.
#' @param \dots other parameters used by \code{\link[utils]{read.table}} function
#'
#' @return return a \code{\link[lumi]{LumiBatch-class}} object
#' @author Mark Cowley, with contributions from Mark Pinese, David Eby.
#' @seealso \code{\link{preprocess.illumina.idat}} \code{\link[lumi]{lumiR}}
#' @export
#' @importFrom lumi lumiR
#' @importClassesFrom lumi LumiBatch
#' @examples
#' library(lumi)
#' path <- system.file("extdata", package="lumidat")
#' files <- c("5356583020_A_Grn.idat", "5356583020_B_Grn.idat")
#' manifestfile <- system.file("extdata", "HumanHT-12_V3_0_R1_99999999.txt", package="lumidat")
#' res <- lumiR.idat(files, path, manifestfile, probeID="NuID")
#' res
#' 
lumiR.idat <- function(files=NULL, path=NULL, probeID=c("ArrayAddressID", "ProbeID", "Sequence", "NuID"), manifestfile=NULL,
  detectionTh=0.01, na.rm=TRUE, parseColumnName=FALSE, checkDupId=TRUE, QC=TRUE, 
  columnNameGrepPattern=list(exprs='AVG_SIGNAL', se.exprs='BEAD_STD', detection='DETECTION', beadNum='Avg_NBEADS'), 
  verbose=TRUE, memory="-Xmx1024m", ...) {
	
	outdir <-  tempdir()
	files <- preprocess.illumina.idat(files=files, path=path, probeID=probeID, manifestfile=manifestfile, outdir=outdir, verbose=verbose, collapseMode="none", prefix=NULL, backgroundCorrect=FALSE, memory=memory)
	# files[1] = "<dir>/Sample Probe Profile.txt", files[2] = "<dir>/Control Probe Profile.txt"
	on.exit(unlink(files))	

	if(!file.exists(files[1])) stop("Error occurred during preprocess.illumina.idat.")
	
	# execute the normal lumiR function from the 'lumi' package
	result <- lumiR(files[1], detectionTh=detectionTh, na.rm=na.rm, convertNuID=FALSE, lib.mapping=NULL, dec='.', parseColumnName=parseColumnName, checkDupId=checkDupId, QC=QC, columnNameGrepPattern=columnNameGrepPattern, inputAnnotation=FALSE, annotationColumn="SYMBOL", verbose=verbose, ...)
	
	return( result )
}
