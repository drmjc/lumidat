\name{read.illumina.manifest}
\alias{read.illumina.manifest}
\title{Import an Illumina manifest file...}
\usage{read.illumina.manifest(f, verbose=FALSE)
}
\description{Import an Illumina manifest file}
\details{Import an Illumina manifest file. Only manifest files in TXT format are supported, ie, BGX files can't be read.
These files can be downloaded from Illumina [1,2]}
\value{a list with 2 elements, probes and controls. Both elements are \dQuote{data.frame}'s with as
many columns as found in the manifest file}
\seealso{[1] \url{http://www.switchtoi.com/annotationfiles.ilmn}

[2] \url{http://www.switchtoi.com/annotationprevfiles.ilmn}}
\author{Mark Cowley, 2011-06-02}
\arguments{\item{f}{the path to the manifest file name}
\item{verbose}{logical if \dQuote{TRUE} print informative messages. Default \dQuote{FALSE}.}
}
\examples{manifestfile <- system.file("extdata", "HumanHT-12_V3_0_R1_99999999.txt", package="lumidat")
manifest <- read.illumina.manifest(manifestfile, TRUE)
head(manifest$probes)
head(manifest$controls)}

