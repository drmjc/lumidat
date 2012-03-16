#' Preprocess Illumina Gene Expression iDAT files.
#' 
#' This package enables the decryption, and preprocessing of lllumina gene expression iDAT files 
#' (aka version 1 iDAT files). Previously, the only option that Illumina gene expression users 
#' had was to rely upon a miroarray core facility with access to an Illumina Scanner and a copy of 
#' Illumina BeadStudio, or GenomeStudio to pre-process the array data. Our intention is to allow 
#' users to pre-process their probe-level Illumina data using this software, thereby enabling them 
#' to choose from the collection of sophisticated normalisation and pre-processing procedures, which 
#' have recently been demonstrated by [1] to simultaneously improve noise and bias, 
#' via the \eqn{lumi}, or \eqn{limma} pipelines. 
#' We have made every effort to reproduce the GenomeStudio output [2], down to the Detection P-Value
#' calculation, the order of the rows, the background correction procedure, the file headers and
#' ProbeID naming, and sample naming.
#' 
#' The bulk of the work is carried out by a Java program, which aims to reproduce the behaviour of 
#' Illumina GenomeStudio, thereby producing \dQuote{Sample Probe Profile.txt} and \dQuote{Control Probe 
#' Profile.txt} files. 
#' This Java program is initially based on the \emph{IlluminaExpressionFileCreator} GenePattern module [3], 
#' created by the Broad Institute (most of the coding done by David Eby).
#' 
#' @section Key functions:
#' \code{\link{read.illumina.idat}} provides a number of options for reading iDAT files, 
#' including background correction, summarisation at the gene-level. Unlike GenomeStudio, it
#' provides greater control over how the gene-level summarisation is performed,
#' and provides a number of different options for naming probes, including NuID [4]. 
#'
#' \code{\link{lumiR.idat}} provides a 
#' replacement for \code{\link[lumi]{lumiR}}, and allows users to use the \emph{lumi} pipeline starting 
#' with iDAT files. \code{\link{read.ilmn.idat}} provides a replacement for 
#'
#' \code{\link[limma]{read.ilmn}}, and allows users to use the \emph{limma} pipeline starting with iDAT 
#' files.
#' 
#' This package will NOT read iDAT files from Illumina Infinium SNP arrays (aka Version 3 iDAT files); 
#' however the \eqn{crlmm} package can read them.
#' 
#' Array manifest files are required, and can be downloaded from Illumina [5], [6]. These files 
#'  describe each of the probes on the array, including the probe sequence, as well as probe-level
#'  annotations. See \code{\link{read.illumina.idat}} for more info.
#'
#' \tabular{ll}{
#' Package: \tab lumidat\cr
#' Type: \tab Package\cr
#' Version: \tab 1.0.3\cr
#' Date: \tab 2012-03-16\cr
#' License: \tab GenePattern license\cr
#' LazyLoad: \tab yes\cr
#' }
#'
#' @seealso \code{\link{read.illumina.idat}}, R\code{\link{lumiR.idat}}, \code{\link[limma]{read.ilmn}}
#'
#' @name lumidat-package
#' @aliases lumidat
#' @docType package
#' @author Mark Cowley <m.cowley@@garvan.org.au>, with contributions from Mark Pinese and David Eby
#' @references
#' [1] Shi, W., Oshlack, A., & Smyth, G. K. (2010). Optimizing the noise versus bias trade-off for 
#' Illumina whole genome expression BeadChips Nucleic acids research, 38(22), e204. (doi:10.1093/nar/gkq871) \cr
#' [2] Genome Studio Gene Expression Module User Guide, version 1.0, Illumina \cr
#' [3] \url{http://genepattern.broadinstitute.org} \cr
#' [4] Du, P., Kibbe, W. A., & Lin, S. M. (2007). nuID: a universal naming scheme of oligonucleotides
#' for illumina, affymetrix, and other microarrays Biol Direct, 2, 16. (doi:10.1186/1745-6150-2-16) \cr
#' [5] \url{http://www.switchtoi.com/annotationfiles.ilmn} \cr
#' [6] \url{http://www.switchtoi.com/annotationprevfiles.ilmn} \cr
#' @examples
#' \dontrun{example(GenomeStudioFileCreator)}
NULL
