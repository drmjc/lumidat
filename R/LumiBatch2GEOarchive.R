#' Convert a LumiBatch into a GEOarchive formatted file
#'
#' To upload Illumina gene expression data to GEO requires a specific 
#' format: \url{http://www.ncbi.nlm.nih.gov/geo/info/geo_illu.html}\cr
#' Briefly, you need both the unnormalised, and normalised data, both
#' of which must include the expression level, and detection p-val for 
#' each probe, using ILMN_ -style probe ID's.
#' \code{LumiBatch} objects conform to this style already.
#' 
#' @section Probe naming:
#' Due to there being a number of different naming conventions for Illumina 
#' probes (see ?\code{\link{preprocess.illumina.idat}}) the probe
#' ID style can differ depending on the \code{probeID} type selected either by
#' Illumina's GenomeStudio, or one of the \code{lumi} methods provided by
#' this package: \code{\link{preprocess.illumina.idat}}, \code{\link{lumiR.idat}}.
#' If your probe ID's differ from the \dQuote{ILMN_} style, then you will
#' either need to convert them to ILMN_ style, using an appropriate
#' probe-level annotation package, eg \code{illuminaHumanv4.db}, or find
#' a GEO Platform which uses your identifiers.
#' 
#' This function pulls the expression levels and detection pvalues from
#' a \code{LumiBatch} object and writes an xls file.
#' 
#' @section ids argument:
#' If you'd like to reorder the output rows, set the \code{ids}
#' argument to be a character vector of probe ID's. You can also use \code{ids} to
#' specify a subset of probes to report. In addition, \code{ids} may contain probes
#' that are not present in \code{featureNames(x)}; in this case these probes
#' will be in effect added to \code{x}, and exported as \code{\dQuote{null}}'s.
#' Yes, this is a bit of a hack, because it's hard to add missing probes to
#' a LumiBatch object retrospectively.
#' 
#' @note This shold work for \code{ExpressionSet} objects, but that's currently
#' untested.
#' 
#' @param x a \code{LumiBatch} object
#' @param file the path to the output file. it should end in tsv or txt
#' @param ids an optional vector of probe ID's. Default=NULL to use the
#'   featureNames within \code{x}. See Details.
#' @param round.digits The number of digits to round the data to. No point
#'   exporting anything much more precise than 4-5 decimal places. Set to
#'   \code{NULL} to skip this.
#' @return nothing. it writes a tab-delimited file.
#' @author Mark Cowley, 2012-03-29
#' @export
#' @importClassesFrom lumi LumiBatch
#' @importFrom Biobase exprs featureNames sampleNames
#' @importFrom lumi detection
#' @examples
#' \dontrun{
#' LumiBatch2GEOarchive(x.raw, "geosub-unnorm.xls")
#' LumiBatch2GEOarchive(x.norm, "geosub-norm.xls")
#' }
LumiBatch2GEOarchive <- function(x, file, ids=NULL, round.digits=5) {
	inherits(x, "LumiBatch") || stop("x must be a LumiBatch object")
	!missing(file) && is.character(file) && length(file) == 1 || stop("file must be the path to an output file")
	
	if( is.null(ids) ) ids <- featureNames(x)
	else if( any(duplicated(ids)) ) stop("ids can't contain duplicates")
	
	dat <- exprs(x)
	pval <- detection(x)
	# I think LumiBatch requires the following to be true, but check it anyway.
	nrow(dat) == nrow(pval) || stop("number of features in exprs and detection slots must be the same")
	ncol(dat) == ncol(pval) || stop("number of samples in exprs and detection slots must be the same")
	
	# do we need to add any rows due to ids being different to featureNames(x) ??
	missing.ids <- setdiff(ids, featureNames(x))
	if( length(missing.ids) > 0 ) {
		blanks <- matrix(NA, length(missing.ids), ncol(dat))
		dimnames(blanks) <- list(missing.ids, sampleNames(x))
		dat  <- rbind(dat, blanks)
		pval <- rbind(pval, blanks)
	}
	
	# limit & reorder dat, pval to 'ids'
	dat <- dat[match(ids, rownames(dat)), ]
	pval <- pval[match(ids, rownames(pval)), ]
	
	# round
	if( !is.null(round.digits) && round.digits > 0 ) {
		dat <- round(dat, digits=round.digits)
		pval <- round(pval, digits=round.digits)
	}
	
	# collate these data, sample 1, pval, sample 2, pval, ..., sample N, pval
	res <- data.frame(ID_REF=ids, matrix(NA, nrow=nrow(dat), ncol=ncol(dat)*2))
	res[seq(2,ncol(res),by=2)] <- dat
	res[seq(3,ncol(res),by=2)] <- pval

	# fix colnames
	cn <- c("ID_REF", rep(sampleNames(x),each=2))
	cn[seq(3,ncol(res),by=2)] <- "Detection Pval"
	colnames(res) <- cn
	
	# write table
	write.table(res, file, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE, na="null")
}
