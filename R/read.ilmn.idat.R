#' @include preprocess.illumina.idat.R
NULL

#' Read Illumina gene expression iDAT files.
#' 
#' This function can decrypt Illumina gene expression iDAT files (aka version 1 iDAT files). 
#' It will temporarily create GenomeStudio-compatible output files using \code{\link{preprocess.illumina.idat}},
#' and then run \code{\link[limma]{read.ilmn}}.
#' 
#' \code{limma::read.ilmn} allows you control over importing the control data via the \code{ctrlpath} parameter.
#' In read.ilmn.idat, you use the \code{controls} parameter: if \code{controls=TRUE}, then both the gene-probes, 
#' and the control-probes will be imported into an \code{\link[limma]{EListRaw-class}} object. 
#'
#' See \code{\link{preprocess.illumina.idat}} for more details on Illumina arrays, manifest files, 
#' probeID naming options, collapsing probes to genes, and background correcting. Note that limma
#' provides \code{\link[limma]{neqc}} for doing background correction (so we usually set 
#' \code{backgroundCorrect=FALSE}). Also, we usually prefer to
#' collapse from probes to genes after: normalizing, and filtering on probe quality, and expression level,
#' thus we usually select \code{collapseMode="none"}.
#' See \code{\link[limma]{read.ilmn}} for more details on \code{expr}, \code{other.columns} parameters.
#' We do not allow access to the annotation parameter, as lumidat doesn't add any additional columns
#' which could be used for annotating probes (ie, not the \dQuote{SYMBOL} parameter as listed in \code{?read.ilmn})
#'
#' @inheritParams preprocess.illumina.idat
#' @param controls logical: if \code{TRUE} (the default), gene-probes and the control-probes will be
#'  imported; if \code{FALSE}, only the gene-probes will be imported. See 
#' \code{\link[limma]{read.ilmn}} for more details.
#' @param expr character string giving the keyword in the names of the expression intensity columns.
#'  default=\code{AVG_Signal}; other options are \dQuote{MIN_Signal}, \dQuote{MAX_Signal}.
#'  See \code{\link[limma]{read.ilmn}} for more details.
#' @param other.columns character vector giving the keywords in the names of extra columns 
#' required, such as \dQuote{Detection}, \dQuote{Avg_NBEADS}, \dQuote{BEAD_STDEV} etc. Each keyword corresponds 
#' to one type of column. See \code{\link[limma]{read.ilmn}} for more details.
#' @param \dots any other parameters are passed on to \code{\link[limma]{read.columns}}.
#' 
#' @return An \code{\link[limma]{EListRaw-class}} object with the following components:
#' \item{E}{numeric matrix of raw intensities.}
#' \item{genes}{data.frame of probe annotation.}
#' \item{targets}{data.frame of sample information.}
#' \item{other}{list of other column data.}
#' 
#' @seealso \code{\link[limma]{read.ilmn}}, \code{\link{preprocess.illumina.idat}}
#' @export
#' @importFrom limma read.ilmn
#' @importClassesFrom limma EListRaw
#' @author Mark Cowley \email{m.cowley@@garvan.org.au}
#' @examples
#' # iDAT files+path as input
#' path <- system.file("extdata", package="lumidat")
#' files <- c("5356583020_A_Grn.idat", "5356583020_B_Grn.idat")
#' manifestfile <- system.file("extdata", "HumanHT-12_V3_0_R1_99999999.txt", package="lumidat")
#' res <- read.ilmn.idat(files, path, manifestfile=manifestfile, probeID="NuID", controls=TRUE)
#' res
#' 
#' # zipfile as input
#' path <- system.file("extdata", package="lumidat")
#' files <- c("5356583020_A_Grn.idat", "5356583020_B_Grn.idat")
#' manifestfile <- system.file("extdata", "HumanHT-12_V3_0_R1_99999999.txt", package="lumidat")
#' zipfile <- tempfile(fileext=".zip")
#' zip(zipfile, file.path(path, files), flags="-r9Xq")
#' res <- read.ilmn.idat(zipfile=zipfile, manifestfile=manifestfile, probeID="ProbeID", controls=FALSE)
#' res
#' 
#' # clm file
#' clmfile <- system.file("extdata", "5356583020.clm", package="lumidat")
#' res <- read.ilmn.idat(zipfile=zipfile, manifestfile=manifestfile, probeID="ProbeID", controls=FALSE, clmfile=clmfile)
#' res$targets
#' #    SampleNames
#' # 1            A
#' # 2            B
#' 
#' \dontrun{
#' # Get the Human HT12 V4 manifest file:
#' manifestfile <- download_illumina_manifest_file("HumanHT-12_V4_0_R2_15002873_B", "txt")
#' }
#' 
read.ilmn.idat <- function (
	# preprocess.illumina.idat arguments:
	files=NULL, path=NULL, zipfile=NULL, clmfile=NULL,
	probeID=c("ProbeID", "ArrayAddressID", "Sequence", "NuID"), manifestfile=NULL,
	collapseMode=c("none", "max", "median", "mean"), 
	backgroundCorrect=FALSE,
	controls=TRUE,
	# limma::read.ilmn arguments
	expr="AVG_Signal", other.columns="Detection", 
	# preprocess.illumina.idat arguments:
	verbose=FALSE, memory="-Xmx1024m", ...) {

	collapseMode <- match.arg(collapseMode)
	probeID <- match.arg(probeID)

	# since I force some settings when calling limma::read.ilmn below, these parameters are not allowed:
	forbidden_params <- intersect(
		c("probeid", "annotation", "sep"),
		names(list(...))
	)
	if( length(forbidden_params) ) stop(sprintf("You can't specify these parameters: %s", paste(forbidden_params, collapse=", ")))
	
	outdir <-  tempdir()
	
	files <- preprocess.illumina.idat(
		files=files, path=path, zipfile=zipfile, manifestfile=manifestfile, clmfile=clmfile, 
		probeID=probeID, collapseMode=collapseMode, backgroundCorrect=backgroundCorrect, 
		outdir=outdir, prefix=NULL, verbose=verbose, memory=memory
	)
	# files[1] = "<dir>/Sample Probe Profile.txt"
	# files[2] = "<dir>/Control Probe Profile.txt"
	on.exit(unlink(files))	

	all(file.exists(files)) || stop("Error occurred during preprocess.illumina.idat.")
	
	if( backgroundCorrect && controls ) {
		message("careful setting both backgroundCorrect=TRUE and controls=TRUE; there's the real potential that you'll end up background correcting again within limma::neqc for example.")
	}
	
	if( controls )
		result <- read.ilmn(files=files[1], ctrlfiles=files[2], probeid="ProbeID", annotation="TargetID", expr = expr, other.columns = other.columns, sep="\t", verbose=verbose, ...)
	else
		result <- read.ilmn(files=files[1],                     probeid="ProbeID", annotation="TargetID", expr = expr, other.columns = other.columns, sep="\t", verbose=verbose, ...)
	
	return( result )
}