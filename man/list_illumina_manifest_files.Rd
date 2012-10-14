\name{list_illumina_manifest_files}
\alias{list_illumina_manifest_files}
\title{List the available Illumina manifest files.}
\usage{
  list_illumina_manifest_files(all = FALSE,
    verbose = FALSE)
}
\arguments{
  \item{all}{logical: if \code{TRUE}, then report all
  current and archived manifest names;
  default=\code{FALSE}, where only the current manifest
  names are reported}

  \item{verbose}{logical: if \code{TRUE} prints info when
  downloading the HTML pages from Illumina.}
}
\value{
  a character vector of manifest names
}
\description{
  Note, you must have an internet connection for this to
  work.
}
\examples{
\dontrun{
#
# output as of 2012-10-14
#
list_illumina_manifest_files()
#  [1] "HumanHT-12_V4_0_R2_15002873_B_WGDASL"
#  [2] "HumanHT-12_V4_0_R2_15002873_B"
#  [3] "HumanHT-12_V4_0_R1_15002873_B"
#  [4] "HumanHT-12_V3_0_R2_11283641_A"
#  [5] "HumanHT-12_V3_0_R3_11283641_A"
#  [6] "HumanRef-8_V2_0_R4_11223162_A"
#  [7] "HUMANREF-8_V3_0_R1_11282963_A_WGDASL"
#  [8] "HumanRef-8_V3_0_R2_11282963_A"
#  [9] "HumanRef-8_V3_0_R3_11282963_A"
# [10] "HumanWG-6_V2_0_R4_11223189_A"
# [11] "HumanWG-6_V3_0_R2_11282955_A"
# [12] "HumanWG-6_V3_0_R3_11282955_A"
# [13] "MouseRef-8_V1_1_R4_11234312_A"
# [14] "MouseRef-8_V2_0_R2_11278551_A"
# [15] "MouseRef-8_V2_0_R3_11278551_A"
# [16] "MouseWG-6_V1_1_R4_11234304_A"
# [17] "MouseWG-6_V2_0_R2_11278593_A"
# [18] "MouseWG-6_V2_0_R3_11278593_A"
# [19] "RatRef-12_V1_0_R5_11222119_A"
# [20] "HumanMI_V1_R2_XS0000122-MAP"
# [21] "HumanMI_V2_R0_XS0000124-MAP"
# [22] "MouseMI_V1_R2_XS0000127-MAP"
# [23] "MouseMI_V2_R0_XS0000129-MAP"

list_illumina_manifest_files(all=TRUE)
# another 20 arrays
}
}
\author{
  Mark Cowley, 2012-10-14
}
\references{
  [1] \url{http://www.switchtoi.com/annotationfiles.ilmn}
  \cr [2]
  \url{http://www.switchtoi.com/annotationprevfiles.ilmn}
  \cr
}
\seealso{
  \code{\link{download_illumina_manifest_file}}
}

