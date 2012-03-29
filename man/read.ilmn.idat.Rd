\name{read.ilmn.idat}
\alias{read.ilmn.idat}
\title{Read Illumina gene expression iDAT files.}
\usage{
  read.ilmn.idat(files = NULL, path = NULL,
    probeID = c("ArrayAddressID", "ProbeID", "Sequence", "NuID"),
    manifestfile = NULL, controls = TRUE,
    other.columns = NULL, verbose = TRUE,
    memory = "-Xmx1024m", ...)
}
\arguments{
  \item{controls}{logical: if TRUE (the default),
  gene-probes and the control-probes will be imported; if
  FALSE, only the gene-probes will be imported. See
  \code{\link[limma]{read.ilmn}} for more details.}

  \item{other.columns}{character vector giving the keywords
  in the names of extra columns required, such as
  "Detection", "Avg_NBEADS", "BEAD_STDEV" etc. Each keyword
  corresponds to one type of column. See
  \code{\link[limma]{read.ilmn}} for more details.}

  \item{\dots}{any other parameters are passed on to
  \code{\link[limma]{read.columns}}.}

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
  An \code{\link[limma]{EListRaw-class}} object with the
  following components: \item{E}{numeric matrix of raw
  intensities.} \item{genes}{data.frame of probe
  annotation.} \item{targets}{data.frame of sample
  information.} \item{other}{list of other column data.}
}
\description{
  This function can decrypt Illumina gene expression iDAT
  files (aka version 1 iDAT files). It will temporarily
  create GenomeStudio-compatible output files, and then run
  \code{\link[limma]{read.ilmn}}.
}
\details{
  If \code{controls=TRUE}, then both the gene-probes, and
  the control-probes will be imported into an
  \code{\link[limma]{EListRaw-class}} object.

  See \code{\link{preprocess.illumina.idat}} for more
  details on Illumina arrays, manifest files, and probeID
  naming options. See \code{\link[limma]{read.ilmn}} for
  more details on \sQuote{other.columns} parameter.
}
\examples{
path <- system.file("extdata", package="lumidat")
files <- c("5356583020_A_Grn.idat", "5356583020_B_Grn.idat")
manifestfile <- system.file("extdata", "HumanHT-12_V3_0_R1_99999999.txt", package="lumidat")
res <- read.ilmn.idat(files, path, manifestfile, probeID="NuID", controls=TRUE)
}
\author{
  Mark Cowley \email{m.cowley@garvan.org.au}
}
\seealso{
  \code{\link[limma]{read.ilmn}},
  \code{\link{preprocess.illumina.idat}}
}

