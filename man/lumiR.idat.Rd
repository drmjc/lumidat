\name{lumiR.idat}
\alias{lumiR.idat}
\title{Read Illumina gene expression iDAT files.}
\usage{
  lumiR.idat(files = NULL, path = NULL,
    probeID = c("ArrayAddressID", "ProbeID", "Sequence", "NuID"),
    manifestfile = NULL, detectionTh = 0.01, na.rm = TRUE,
    parseColumnName = FALSE, checkDupId = TRUE, QC = TRUE,
    columnNameGrepPattern = list(exprs = "AVG_SIGNAL", se.exprs = "BEAD_STD", detection = "DETECTION", beadNum = "Avg_NBEADS"),
    verbose = TRUE, memory = "-Xmx1024m", ...)
}
\arguments{
  \item{detectionTh}{the p-value threshold of determining
  detectability of the expression.  See more details in
  \code{\link[lumi]{lumiQ}}.}

  \item{na.rm}{logical: remove \code{NA}?}

  \item{parseColumnName}{logical: parse the column names
  and retrieve the sample information? (Assume the sample
  information is separated by \dQuote{\_}.)}

  \item{checkDupId}{logical: check duplicated TargetIDs or
  ProbeIds? The duplicated ones will be averaged.}

  \item{QC}{logical:do quality control assessment after
  read in the data?}

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

  \item{probeID}{This controls which value to identify each
  probe by. Allowable values are \dQuote{ArrayAddressID},
  \dQuote{ProbeID}, \dQuote{Sequence}, \dQuote{NuID}.}

  \item{manifestfile}{The full path to the Array manifest
  file in TXT format.}

  \item{verbose}{logical, if TRUE, print informative
  messages}

  \item{memory}{The maximum Java memory heapsize. Default
  -Xmx1024m, which allows 1GB of RAM.}
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
  details on Illumina arrays, manifest files, and probeID
  naming options. See \code{\link[lumi]{lumiR}} for more
  details on all paramaters from \code{detectionTh}
  onwards.
}
\examples{
library(lumi)
path <- system.file("extdata", package="lumidat")
files <- c("5356583020_A_Grn.idat", "5356583020_B_Grn.idat")
manifestfile <- system.file("extdata", "HumanHT-12_V3_0_R1_99999999.txt", package="lumidat")
res <- lumiR.idat(files, path, manifestfile, probeID="NuID")
res
}
\author{
  Mark Cowley, with contributions from Mark Pinese, David
  Eby.
}
\seealso{
  \code{\link{preprocess.illumina.idat}}
  \code{\link[lumi]{lumiR}}
}

