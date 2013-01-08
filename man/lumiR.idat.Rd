\name{lumiR.idat}
\alias{lumiR.idat}
\title{Read Illumina gene expression iDAT files.}
\usage{
  lumiR.idat(files = NULL, path = NULL, zipfile = NULL,
    clmfile = NULL,
    probeID = c("ProbeID", "ArrayAddressID", "Sequence", "NuID"),
    manifestfile = NULL,
    collapseMode = c("none", "max", "median", "mean"),
    backgroundCorrect = FALSE, controls = TRUE,
    detectionTh = 0.01, na.rm = TRUE,
    parseColumnName = FALSE, checkDupId = TRUE, QC = TRUE,
    columnNameGrepPattern = list(exprs = "AVG_SIGNAL", se.exprs = "BEAD_STD", detection = "DETECTION", beadNum = "Avg_NBEADS"),
    verbose = FALSE, memory = "-Xmx1024m", ...)
}
\arguments{
  \item{controls}{logical: if \code{TRUE}, the controlData
  slot will be filled, via
  \code{lumi::\link[lumi]{addControlData2lumi}}. We
  recommend this to be \code{TRUE}.}

  \item{detectionTh}{the p-value threshold of determining
  detectability of the expression.  See more details in
  \code{\link[lumi]{lumiQ}}.}

  \item{na.rm}{logical: remove \code{NA}?}

  \item{parseColumnName}{logical: parse the column names
  and retrieve the sample information? (Assume the sample
  information is separated by \dQuote{\_} from the column
  type, eg \dQuote{5356583020_A_AVG_Signal}.)}

  \item{checkDupId}{logical: check duplicated TargetIDs or
  ProbeIds? The duplicated ones will be averaged.}

  \item{QC}{logical: do quality control assessment after
  reading in the data?}

  \item{columnNameGrepPattern}{list of named
  character(1)'s: the grep patterns used to determine which
  column name corresponds to which slot. Eg the column
  named \dQuote{AVG_SIGNAL} will be put into the
  \code{exprs} slot.}

  \item{\dots}{other parameters used by
  \code{\link[utils]{read.table}} function}

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

  \item{memory}{The maximum Java memory heapsize. eg:
  "-Xmx1024m", or "-Xmx4g", which reserves 1GB or 4G of
  RAM, respectively. default="-Xmx1024m"}
}
\value{
  return a \code{\link[lumi]{LumiBatch-class}} object
}
\description{
  This function can decrypt Illumina gene expression iDAT
  files (aka version 1 iDAT files). It will temporarily
  create GenomeStudio-compatible output files, and then run
  \code{\link[lumi]{lumiR}}.
}
\details{
  See \code{\link{preprocess.illumina.idat}} for more
  details on Illumina arrays, manifest files, probeID
  naming options, collapsing probes to genes, and
  background correcting. Note that lumi provides
  \code{\link[lumi]{lumiB}} for doing background correction
  (so we usually set \code{backgroundCorrect=FALSE}). Also,
  we usually prefer to collapse from probes to genes after:
  normalizing, and filtering on probe quality, and
  expression level, thus we usually select
  \code{collapseMode="none"}.

  See \code{\link[lumi]{lumiR}} for more details on all
  parameters from \code{detectionTh} onwards.
}
\examples{
# iDAT files+path as input
library(lumi)
path <- system.file("extdata", package="lumidat")
files <- c("5356583020_A_Grn.idat", "5356583020_B_Grn.idat")
manifestfile <- system.file("extdata", "HumanHT-12_V3_0_R1_99999999.txt", package="lumidat")
res <- lumiR.idat(files, path, manifestfile=manifestfile, probeID="ProbeID")
res

# zipfile as input
path <- system.file("extdata", package="lumidat")
files <- c("5356583020_A_Grn.idat", "5356583020_B_Grn.idat")
manifestfile <- system.file("extdata", "HumanHT-12_V3_0_R1_99999999.txt", package="lumidat")
zipfile <- tempfile(fileext=".zip")
zip(zipfile, file.path(path, files), flags="-r9Xq")
res <- lumiR.idat(zipfile=zipfile, manifestfile=manifestfile, probeID="ProbeID", verbose=FALSE)
res

# clm file
clmfile <- system.file("extdata", "5356583020.clm", package="lumidat")
res <- lumiR.idat(zipfile=zipfile, manifestfile=manifestfile, probeID="ProbeID", verbose=FALSE, clmfile=clmfile)
sampleNames(res)
# [1] "A" "B"

\dontrun{
# Get the Human HT12 V4 manifest file:
manifestfile <- download_illumina_manifest_file("HumanHT-12_V4_0_R2_15002873_B", "txt")
}
}
\author{
  Mark Cowley, with contributions from Mark Pinese, David
  Eby.
}
\seealso{
  \code{\link{preprocess.illumina.idat}}
  \code{\link[lumi]{lumiR}}
}

