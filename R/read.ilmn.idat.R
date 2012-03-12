#' @include read.illumina.idat.R
NULL

#' Read Illumina gene expression iDAT files.
#' 
#' This function can decrypt Illumina gene expression iDAT files (aka version 1 iDAT files). 
#' It will temporarily create GenomeStudio-compatible output files, and then run 
#' \code{\link[limma]{read.ilmn}}.
#' 
#' If \eqn{controls=TRUE}, then both the gene-probes, and the control-probes will be imported 
#' into an \code{\link[limma]{EListRaw-class}} object. 
#'
#' See \code{\link{read.illumina.idat}} for more details on Illumina arrays, manifest files, and probeID naming options.
#' See \code{\link[limma]{read.ilmn}} for more details on \sQuote{other.columns} parameter.
#'
#' @param files character vector of >= 1 iDAT file name
#' @param path character string giving the directory containing the iDAT files. The default 
#' is the current working directory.
#' @param probeID This controls which value to identify each probe by. Allowable values are 
#' \dQuote{ArrayAddressID}, \dQuote{ProbeID}, \dQuote{Sequence}, \dQuote{NuID}. See \code{\link{read.illumina.idat}}
#'  for more details.
#' @param manifestfile The full path to the Array manifest file in TXT format.
#' @param controls logical: if TRUE (the default), gene-probes and the control-probes will be
#'  imported; if FALSE, only the gene-probes will be imported. See 
#' \code{\link[limma]{read.ilmn}} for more details.
#' @param other.columns character vector giving the keywords in the names of extra columns 
#' required, such as "Detection", "Avg_NBEADS", "BEAD_STDEV" etc. Each keyword corresponds 
#' to one type of columns. See \code{\link[limma]{read.ilmn}} for more details.
#' @param verbose logical, \eqn{TRUE} to report names of profile files being read.
#' @param memory the maximum amount of memory to use. see \code{\link{read.illumina.idat}}.
#' @param ... any other parameters are passed on to \code{\link[limma]{read.columns}}.
#' @return An \code{\link[limma]{EListRaw-class}} object with the following components:
#' \tabular{ll}{
#' \code{E} \tab numeric matrix of raw intensities.\cr
#' \code{genes} \tab data.frame of probe annotation.\cr
#' \code{targets} \tab data.frame of sample information.\cr
#' \code{other} \tab list of other column data.\cr
#' }
#' @seealso \code{\link[limma]{read.ilmn}}, \code{\link{read.illumina.idat}}
#' @export
#' @author Mark Cowley \email{m.cowley@@garvan.org.au}
#' @examples
#' path <- system.file("extdata", package="lumidat")
#' files <- c("5356583020_A_Grn.idat", "5356583020_B_Grn.idat")
#' manifestfile <- system.file("extdata", "HumanHT-12_V3_0_R1_99999999.txt", package="lumidat")
#' res <- read.ilmn.idat(files, path, manifestfile, probeID="NuID", controls=TRUE)
#'
read.ilmn.idat <- function (files=NULL, path=NULL, probeID=c("ArrayAddressID", "ProbeID", "Sequence", "NuID"), manifestfile=NULL, controls=TRUE, other.columns=NULL, verbose=TRUE, memory="-Xmx1024m", ...) {
	require(limma)
	
	outdir <-  tempdir()
	files <- read.illumina.idat(files=files, path=path, probeID=probeID, manifestfile=manifestfile, outdir=outdir, verbose=verbose, collapseMode="none", prefix=NULL, backgroundCorrect=FALSE, memory=memory)
	# files[1] = "<dir>/Sample Probe Profile.txt", files[2] = "<dir>/Control Probe Profile.txt"
	on.exit(unlink(files))	

	if( !all(file.exists(files)) ) stop("Error occurred during read.illumina.idat.")
	
	if( !controls ) files[2] <- NULL
	
	result <- read.ilmn(files=files[1], ctrlfiles=files[2], probeid="ProbeID", annotation="TargetID", expr = "AVG_Signal", sep="\t", verbose=verbose, ...)
	
	return( result )
}