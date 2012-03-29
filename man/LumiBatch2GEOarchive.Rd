\name{LumiBatch2GEOarchive}
\alias{LumiBatch2GEOarchive}
\title{Convert a LumiBatch into a GEOarchive formatted file}
\usage{
  LumiBatch2GEOarchive(x, file, ids = NULL,
    round.digits = 5)
}
\arguments{
  \item{x}{a \code{LumiBatch} object}

  \item{file}{the path to the output file. it should end in
  tsv or txt}

  \item{ids}{an optional vector of probe ID's. Default=NULL
  to use the featureNames within \code{x}. See Details.}

  \item{round.digits}{The number of digits to round the
  data to. No point exporting anything much more precise
  than 4-5 decimal places. Set to \code{NULL} to skip
  this.}
}
\value{
  nothing. it writes a tab-delimited file.
}
\description{
  To upload Illumina gene expression data to GEO requires a
  specific format:
  \url{http://www.ncbi.nlm.nih.gov/geo/info/geo_illu.html}\cr
  Briefly, you need both the unnormalised, and normalised
  data, both of which must include the expression level,
  and detection p-val for each probe, using ILMN_ -style
  probe ID's. \code{LumiBatch} objects conform to thus
  style already.
}
\note{
  This shold work for \code{ExpressionSet} objects, but
  that's currently untested.
}
\section{Probe naming}{
  Due to there being a number of different naming
  conventions for Illumina probes (see
  ?\code{\link{preprocess.illumina.idat}}) the probe ID
  style can differ depending on the \code{probeID} type
  selected either by Illumina's GenomeStudio, or one of the
  \code{lumi} methods provided by this package:
  \code{\link{preprocess.illumina.idat}},
  \code{\link{lumiR.idat}}. If your probe ID's differ from
  the \dQuote{ILMN_} style, then you will either need to
  convert them to ILMN_ style, using an appropriate
  probe-level annotation package, eg
  \code{illuminaHumanv4.db}, or find a GEO Platform which
  uses your identifiers.

  This function pulls the expression levels and detection
  pvalues from a \code{LumiBatch} object and writes an xls
  file.
}

\section{ids argument}{
  If you'd like to reorder the output rows, set the
  \code{ids} argument to be a character vector of probe
  ID's. You can also use \code{ids} to specify a subset of
  probes to report. In addition, \code{ids} may contain
  probes that are not present in \code{featureNames(x)}; in
  this case these probes will be in effect added to 'x',
  and exported as \code{\dQuote{null}}'s. Yes, this is a
  bit of a hack, because it's hard to add missing probes to
  a LumiBatch object retrospectively.
}
\examples{
\dontrun{
LumiBatch2GEOarchive(x.raw, "geosub-unnorm.xls")
LumiBatch2GEOarchive(x.norm, "geosub-norm.xls")
}
}
\author{
  Mark Cowley, 2012-03-29
}

