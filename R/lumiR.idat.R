#' @include preprocess.illumina.idat.R
NULL

#' Read Illumina gene expression iDAT files.
#' 
#' This function can decrypt Illumina gene expression iDAT files (aka version 1 iDAT files). It will
#' temporarily create GenomeStudio-compatible output files, and then run \code{\link[lumi]{lumiR}}.
#'
#' See \code{\link{preprocess.illumina.idat}} for more details on Illumina arrays, manifest files, 
#' probeID naming options, collapsing probes to genes, and background correcting. Note that lumi
#' provides \code{\link[lumi]{lumiB}} for doing background correction (so we usually set 
#' \code{backgroundCorrect=FALSE}). Also, we usually prefer to
#' collapse from probes to genes after: normalizing, and filtering on probe quality, and expression level,
#' thus we usually select \code{collapseMode="none"}.
#' 
#' See \code{\link[lumi]{lumiR}} for more details on all parameters from \code{detectionTh} onwards.
#' 
#' @inheritParams preprocess.illumina.idat
#' @param controls logical: if \code{TRUE}, the controlData slot will be filled, via 
#'  \code{lumi::\link[lumi]{addControlData2lumi}}. We recommend this to be \code{TRUE}.
#' @param detectionTh the p-value threshold of determining detectability of the expression.
#'  See more details in \code{\link[lumi]{lumiQ}}.
#' @param na.rm logical: remove \code{NA}?
#' @param parseColumnName logical: parse the column names and retrieve the sample 
#' information? (Assume the sample information is separated by \dQuote{\_} from the 
#' column type, eg \dQuote{5356583020_A_AVG_Signal}.)
#' @param checkDupId logical: check duplicated TargetIDs or ProbeIds? The 
#' duplicated ones will be averaged.
#' @param QC logical: do quality control assessment after reading in the data?
#' @param columnNameGrepPattern list of named character(1)'s: the grep patterns used to 
#' determine which column name corresponds to which slot. Eg the column named \dQuote{AVG_SIGNAL}
#' will be put into the \code{exprs} slot.
#' @param \dots other parameters used by \code{\link[utils]{read.table}} function
#'
#' @return return a \code{\link[lumi]{LumiBatch-class}} object
#' @author Mark Cowley, with contributions from Mark Pinese, David Eby.
#' @seealso \code{\link{preprocess.illumina.idat}} \code{\link[lumi]{lumiR}}
#' @export
#' @importFrom lumi lumiR addControlData2lumi
#' @importClassesFrom lumi LumiBatch
#' @examples
#' # iDAT files+path as input
#' library(lumi)
#' path <- system.file("extdata", package="lumidat")
#' files <- c("5356583020_A_Grn.idat", "5356583020_B_Grn.idat")
#' manifestfile <- system.file("extdata", "HumanHT-12_V3_0_R1_99999999.txt", package="lumidat")
#' res <- lumiR.idat(files, path, manifestfile=manifestfile, probeID="ProbeID")
#' res
#' 
#' # zipfile as input
#' path <- system.file("extdata", package="lumidat")
#' files <- c("5356583020_A_Grn.idat", "5356583020_B_Grn.idat")
#' manifestfile <- system.file("extdata", "HumanHT-12_V3_0_R1_99999999.txt", package="lumidat")
#' zipfile <- tempfile(fileext=".zip")
#' zip(zipfile, file.path(path, files), flags="-r9Xq")
#' res <- lumiR.idat(zipfile=zipfile, manifestfile=manifestfile, probeID="ProbeID", verbose=FALSE)
#' res
#'
#' # clm file
#' clmfile <- system.file("extdata", "5356583020.clm", package="lumidat")
#' res <- lumiR.idat(zipfile=zipfile, manifestfile=manifestfile, probeID="ProbeID", verbose=FALSE, clmfile=clmfile)
#' sampleNames(res)
#' # [1] "A" "B"
#' 
#' \dontrun{
#' # Get the Human HT12 V4 manifest file:
#' manifestfile <- download_illumina_manifest_file("HumanHT-12_V4_0_R2_15002873_B", "txt")
#' }
#' 
lumiR.idat <- function(
	# preprocess.illumina.idat arguments:
	files=NULL, path=NULL, zipfile=NULL, clmfile=NULL,
	probeID=c("ProbeID", "ArrayAddressID", "Sequence", "NuID"), manifestfile=NULL,
	collapseMode=c("none", "max", "median", "mean"), 
	backgroundCorrect=FALSE,
	controls=TRUE,
	# lumi::lumiR arguments:
	detectionTh=0.01, na.rm=TRUE, parseColumnName=FALSE, checkDupId=TRUE, QC=TRUE, 
	columnNameGrepPattern=list(exprs='AVG_SIGNAL', se.exprs='BEAD_STD', detection='DETECTION', beadNum='Avg_NBEADS'), 
	# preprocess.illumina.idat arguments:
	verbose=FALSE, memory="-Xmx1024m", ...) {

	collapseMode <- match.arg(collapseMode)
	probeID <- match.arg(probeID)
	
	# since I force some settings when calling lumi::lumiR below, these parameters are not allowed:
	forbidden_params <- intersect(
		c("convertNuID", "lib.mapping", "dec", "inputAnnotation", "annotationColumn"),
		names(list(...))
	)
	if( length(forbidden_params) ) stop(sprintf("You can't specify these parameters: %s", paste(forbidden_params, collapse=", ")))
	
	outdir <-  tempdir()
	files <- preprocess.illumina.idat(
		files=files, path=path, zipfile=zipfile, manifestfile=manifestfile, clmfile=clmfile, 
		probeID=probeID, collapseMode=collapseMode, backgroundCorrect=backgroundCorrect, 
		outdir=outdir, prefix=NULL, verbose=verbose, memory=memory
	)
	# files[1] = "<dir>/Sample Probe Profile.txt", files[2] = "<dir>/Control Probe Profile.txt"
	on.exit(unlink(files))	

	if(!file.exists(files[1])) stop("Error occurred during preprocess.illumina.idat.")
	
	# execute the normal lumiR function from the 'lumi' package
	result <- lumiR(files[1], detectionTh=detectionTh, na.rm=na.rm, convertNuID=FALSE, lib.mapping=NULL, dec='.', parseColumnName=parseColumnName, checkDupId=checkDupId, QC=QC, columnNameGrepPattern=columnNameGrepPattern, inputAnnotation=FALSE, annotationColumn="SYMBOL", verbose=verbose, ...)
	
	if( backgroundCorrect && controls ) {
		message("careful setting both backgroundCorrect=TRUE and controls=TRUE; you should avoid further background correcting using lumi::lumiB for example.")
	}
	
	if( controls ) {
		if( collapseMode == "none" ) {
			result <- addControlData2lumi(files[2], result)
		}
		else {
			message("addControlData2lumi can't handle data that's been collapsed. forcing controls=FALSE")
		}
	}
	
	return( result )
}
