\name{read.ilmn.idat}
\alias{read.ilmn.idat}
\title{Read Illumina gene expression iDAT files.}
\usage{
  read.ilmn.idat(files = NULL, path = NULL, zipfile = NULL,
    clmfile = NULL,
    probeID = c("ProbeID", "ArrayAddressID", "Sequence", "NuID"),
    manifestfile = NULL,
    collapseMode = c("none", "max", "median", "mean"),
    backgroundCorrect = FALSE, controls = TRUE,
    expr = "AVG_Signal", other.columns = "Detection",
    verbose = TRUE, memory = "-Xmx1024m", ...)
}
\arguments{
  \item{controls}{logical: if \code{TRUE} (the default),
  gene-probes and the control-probes will be imported; if
  \code{FALSE}, only the gene-probes will be imported. See
  \code{\link[limma]{read.ilmn}} for more details.}

  \item{expr}{character string giving the keyword in the
  names of the expression intensity columns.
  default=\code{AVG_Signal}; other options are
  \dQuote{MIN_Signal}, \dQuote{MAX_Signal}.  See
  \code{\link[limma]{read.ilmn}} for more details.}

  \item{other.columns}{character vector giving the keywords
  in the names of extra columns required, such as
  \dQuote{Detection}, \dQuote{Avg_NBEADS},
  \dQuote{BEAD_STDEV} etc. Each keyword corresponds to one
  type of column. See \code{\link[limma]{read.ilmn}} for
  more details.}

  \item{\dots}{any other parameters are passed on to
  \code{\link[limma]{read.columns}}.}

  \item{files}{A character vector of at least one file
  name}

  \item{path}{The path to the directories where the files
  are. This defaults to the current working directory.}

  \item{zipfile}{the path to a zip file containing idat
  files. This over-rides the \code{files} argument.}

  \item{clmfile}{\code{NULL}, or the path to a
  GenePattern-CLM file [5]. This represents a mechanism for
  renaming, and reordering the samples in the resulting
  object. Column 1 is the idat name (not path), Column 2 is
  the desired sample name, column 3 is the biological class
  of each sample. The order of the rows specifies the order
  of the samples in the result.}

  \item{probeID}{This controls which value to identify each
  probe by. Allowable values are \dQuote{ArrayAddressID},
  \dQuote{ProbeID}, \dQuote{Sequence}, \dQuote{NuID}.}

  \item{manifestfile}{The full path to the Array manifest
  file in TXT format. See
  \code{\link{list_illumina_manifest_files}} and then
  \code{\link{download_illumina_manifest_file}} to download
  manifest files directly from Illumina.}

  \item{collapseMode}{Collapse probes to genes using the
  specified mode. Valid values are \dQuote{none} (the
  default), \dQuote{max}, \dQuote{median} and \dQuote{mean}
  (the GenomeStudio default).}

  \item{backgroundCorrect}{logical, if TRUE then peform
  background correction on only the gene-level probes; if
  FALSE, do no correction.}

  \item{verbose}{logical, if \code{TRUE}, print informative
  messages}

  \item{memory}{The maximum Java memory heapsize. Default
  -Xmx1024m, which allows 1GB of RAM.}
}
\value{
  An \code{\link[limma]{EListRaw-class}} object with the
  following components: \item{E}{numeric matrix of raw
  intensities.} \item{genes}{data.frame of probe
  annotation.} \item{targets}{data.frame of sample
  information.} \item{other}{list of other column data.}
}
\description{
  This function can decrypt Illumina gene expression iDAT
  files (aka version 1 iDAT files). It will temporarily
  create GenomeStudio-compatible output files using
  \code{\link{preprocess.illumina.idat.R}}, and then run
  \code{\link[limma]{read.ilmn}}.
}
\details{
  \code{limma::read.ilmn} allows you control over importing
  the control data via the \code{ctrlpath} parameter. In
  read.ilmn.idat, you use the \code{controls} parameter: if
  \code{controls=TRUE}, then both the gene-probes, and the
  control-probes will be imported into an
  \code{\link[limma]{EListRaw-class}} object.

  See \code{\link{preprocess.illumina.idat}} for more
  details on Illumina arrays, manifest files, probeID
  naming options, collapsing probes to genes, and
  background correcting. Note that limma provides
  \code{\link[limma]{neqc}} for doing background correction
  (so we usually set \code{backgroundCorrect=FALSE}). Also,
  we usually prefer to collapse from probes to genes after:
  normalizing, and filtering on probe quality, and
  expression level, thus we usually select
  \code{collapseMode="none"}. See
  \code{\link[limma]{read.ilmn}} for more details on
  \code{expr}, \code{other.columns} parameters. We do not
  allow access to the annotation parameter, as lumidat
  doesn't add any additional columns which could be used
  for annotating probes (ie, not the \dQuote{SYMBOL}
  parameter as listed in \code{?read.ilmn})
}
\examples{
# iDAT files+path as input
path <- system.file("extdata", package="lumidat")
files <- c("5356583020_A_Grn.idat", "5356583020_B_Grn.idat")
manifestfile <- system.file("extdata", "HumanHT-12_V3_0_R1_99999999.txt", package="lumidat")
res <- read.ilmn.idat(files, path, manifestfile=manifestfile, probeID="NuID", controls=TRUE)
res

# zipfile as input
path <- system.file("extdata", package="lumidat")
files <- c("5356583020_A_Grn.idat", "5356583020_B_Grn.idat")
manifestfile <- system.file("extdata", "HumanHT-12_V3_0_R1_99999999.txt", package="lumidat")
zipfile <- tempfile(fileext=".zip")
zip(zipfile, file.path(path, files))
res <- read.ilmn.idat(zipfile=zipfile, manifestfile=manifestfile, probeID="ProbeID", controls=FALSE)
res

# clm file
clmfile <- system.file("extdata", "5356583020.clm", package="lumidat")
res <- read.ilmn.idat(zipfile=zipfile, manifestfile=manifestfile, probeID="ProbeID", controls=FALSE, clmfile=clmfile)
sampleNames(res)
# [1] "A" "B"
}
\author{
  Mark Cowley \email{m.cowley@garvan.org.au}
}
\seealso{
  \code{\link[limma]{read.ilmn}},
  \code{\link{preprocess.illumina.idat}}
}

