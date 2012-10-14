\name{download_illumina_manifest_file}
\alias{download_illumina_manifest_file}
\title{Download an Illumina array manifest}
\usage{
  download_illumina_manifest_file(array.name,
    type = c("txt", "bgx"), dir = ".", verbose = FALSE)
}
\arguments{
  \item{array.name}{the \emph{exact} array name, including
  capitalisation, as it appears in [1,2] or
  \code{\link{list_illumina_manifest_files}}.}

  \item{type}{one of \dQuote{txt} or \dQuote{bgx}}

  \item{dir}{the directory to download to. defaults to
  current working dir.}

  \item{verbose}{logical}
}
\value{
  a character vector(1) of the path to the downloaded
  manifest file
}
\description{
  Illumina array manifest files (current and archived) are
  hosted in [1,2], and exist in 2 formats: bgx (compressed)
  and txt (most commonly used). Given the exact name of the
  array, this will download either the bgx or txt formatted
  manifest file to dir.
}
\section{manifest files}{
  These are the current (2012-10-14) array manifest names:
  \tabular{l}{ HumanHT-12_V3_0_R2_11283641_A \cr
  HumanHT-12_V3_0_R3_11283641_A \cr
  HumanHT-12_V4_0_R1_15002873_B \cr
  HumanHT-12_V4_0_R2_15002873_B \cr
  HumanHT-12_V4_0_R2_15002873_B_WGDASL \cr
  HumanMI_V1_R2_XS0000122-MAP \cr
  HumanMI_V2_R0_XS0000124-MAP \cr
  HumanRef-8_V2_0_R4_11223162_A \cr
  HUMANREF-8_V3_0_R1_11282963_A_WGDASL \cr
  HumanRef-8_V3_0_R2_11282963_A \cr
  HumanRef-8_V3_0_R3_11282963_A \cr
  HumanWG-6_V2_0_R4_11223189_A \cr
  HumanWG-6_V3_0_R2_11282955_A \cr
  HumanWG-6_V3_0_R3_11282955_A \cr
  MouseMI_V1_R2_XS0000127-MAP \cr
  MouseMI_V2_R0_XS0000129-MAP \cr
  MouseRef-8_V1_1_R4_11234312_A \cr
  MouseRef-8_V2_0_R2_11278551_A \cr
  MouseRef-8_V2_0_R3_11278551_A \cr
  MouseWG-6_V1_1_R4_11234304_A \cr
  MouseWG-6_V2_0_R2_11278593_A \cr
  MouseWG-6_V2_0_R3_11278593_A \cr
  RatRef-12_V1_0_R5_11222119_A \cr } For the archived ones,
  see \code{\link{list_illumina_manifest_files}}.
}
\examples{
\dontrun{
download_illumina_manifest_file("HumanHT-12_V4_0_R2_15002873_B", "txt")
}
}
\author{
  Mark Cowley, 2012-10-14
}
\seealso{
  \code{\link{list_illumina_manifest_files}}
}

