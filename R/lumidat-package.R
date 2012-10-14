#' Preprocess Illumina Gene Expression iDAT files.
#' 
#' This package enables the decryption, and preprocessing of lllumina gene expression iDAT files 
#' (aka version 1 iDAT files). Previously, the only option that Illumina gene expression users 
#' had was to rely upon a miroarray core facility with access to an Illumina Scanner and a copy of 
#' Illumina BeadStudio, or GenomeStudio to pre-process the array data. Our intention is to allow 
#' users to pre-process their probe-level Illumina data using this software, thereby enabling them 
#' to choose from the collection of sophisticated normalisation and pre-processing procedures, which 
#' have recently been demonstrated by [1] to simultaneously improve noise and bias, 
#' via the \code{lumi}, or \code{limma} pipelines. 
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
#' \code{\link{preprocess.illumina.idat}} provides a number of options for reading iDAT files, 
#' including background correction, summarisation at the gene-level. Unlike GenomeStudio, it
#' provides greater control over how the gene-level summarisation is performed,
#' and provides a number of different options for naming probes, including NuID [4]. 
#'
#' \code{\link{lumiR.idat}} provides a 
#' replacement for \code{\link[lumi]{lumiR}}, and allows users to use the \emph{lumi} pipeline starting 
#' with iDAT files.
#' 
#' \code{\link{read.ilmn.idat}} provides a replacement for 
#' \code{\link[limma]{read.ilmn}}, and allows users to use the \emph{limma} pipeline starting with iDAT 
#' files.
#' 
#' This package will NOT read iDAT files from Illumina Infinium SNP arrays (aka Version 3 iDAT files); 
#' however the \code{crlmm} package can read them.
#' 
#' @section Array manifest files:
#' Array manifest files are required, and can be downloaded from Illumina [5], [6]. These files 
#'  describe each of the probes on the array, including the probe sequence, as well as probe-level
#'  annotations. See \code{\link{preprocess.illumina.idat}} for more info.
#' To make this easier, we have added \code{\link{list_illumina_manifest_files}}, and 
#' \code{\link{download_illumina_manifest_file}}.
#'
#' @section Analysis pipeline:
#' # \cr
#' # import & initial QC \cr
#' # \cr
#' idat.files <- dir("/path/to/idat/files", full.path=TRUE) \cr
#' x.raw <- lumiR.idat(idat.files, probeID="ProbeID", manifestfile="/Volumes/GRIW/ICGCPancreas/icgc_data/icgc_gex/HumanHT-12_V4_0_R2_15002873_B.txt", convertNuID=FALSE,inputAnnotation=FALSE, QC=TRUE, memory="-Xmx4080m") \cr
#' dir.create("QC") \cr
#' doMA <- doPairs <- FALSE \cr
#' plot.lumi.QC.all(x.raw, "QC/01.unnorm/", "raw", "Unnormalised", doMA, doPairs) \cr
#'  \cr
#' # \cr
#' # optional, import a probe-level annotation object \cr
#' # \cr
#' load("Rmisc/IlluminaHT4.anno-2.0.RDa.gz") \cr
#' fData(x.raw) <- IlluminaHT4.anno \cr
#' fvarMetadata(x.raw) <- structure(list( \cr
#' 	labelDescription = c( \cr
#' 		"The Illumina microarray probe identifier",  \cr
#' 		"Array Address code to identify the probe at the bead-level",  \cr
#' 		"Lumi's nuID (universal naming scheme for oligos)",  \cr
#' 		"Quality grade assigned to the probe",  \cr
#' 		"Coding status of the target sequence: intergenic / intronic / transcriptomic? / transcriptomic",  \cr
#' 		"50mer probe sequence",  \cr
#' 		"Genomic coordinates of second best matches between the probe and the genome",  \cr
#' 		"Genomic coordinates of sequences as alignable with the probe as its main target",  \cr
#' 		"Overlapping RepeatMasked sequences",  \cr
#' 		"Overlapping annotated SNPs",  \cr
#' 		"Entrez Gene ID's after remapping probes",  \cr
#' 		"Probe's genomic coordinates (hg19, mm9 or rn4)",  \cr
#' 		"Gene symbol derived by reannotation",  \cr
#' 		"A more descriptive identifier for controls",  \cr
#' 		"An identifier for probes designated as controls by Illumina",  \cr
#' 		"Gene title (derived from EntrezReannotated)", \cr
#' 		"Probe-level Description" \cr
#' 	)),  \cr
#' 	.Names = "labelDescription",  \cr
#' 	row.names = c( \cr
#' 		"ProbeID",  \cr
#' 		"ArrayAddressID", "NuID", "ProbeQuality", "CodingZone", "ProbeSequence",  \cr
#' 		"SecondMatches", "OtherGenomicMatches", "RepeatMask", "ContainsSNP",  \cr
#' 		"EntrezReannotated", "GenomicLocation", "SymbolReannotated",  \cr
#' 		"ReporterGroupName", "ReporterGroupID", "GeneNameReannotated", "Description" \cr
#' 	), \cr
#' 	class = "data.frame" \cr
#' ) \cr
#'  \cr
#' # \cr
#' # potentially exclude 'bad' arrays \cr
#' # \cr
#' x.passedqc <- x.raw \cr
#' plot.lumi.QC.all(x.passedqc, "QC/02.passedqc/", "raw", "Unnormalised", doMA, doPairs) \cr
#'  \cr
#' # \cr
#' # VST transform \cr
#' # \cr
#' x.transformed   <- lumiT(x.passedqc, method="vst") \cr
#' plot.lumi.QC.all(x.transformed, "QC/03.transformed/", "vst", "Transformed", MA=FALSE, pairs=FALSE) \cr
#'  \cr
#' # \cr
#' # normalise. one of rsn, ssn,  \cr
#' # \cr
#' x.norm <- lumiN(x.transformed, "rsn") \cr
#' plot.lumi.QC.all(x.norm, "QC/04.norm/", "rsn", "RSN Normalised", MA=FALSE, pairs=FALSE) \cr
#'  \cr
#' # \cr
#' # average replicates \cr
#' # \cr
#' x.averaged <- average.replicates(x.norm, sub("\\.[12]$", "", sampleNames(x.norm))) \cr
#'  \cr
#' # \cr
#' # filter on probe quality: \cr
#' # \cr
#' Rkeys(illuminaHumanv4PROBEQUALITY) \cr
#' # [1] "No match"    "Bad"         "Perfect***"  "Perfect"     "Perfect****" \cr
#' # [6] "Good****"    "Good"        "Good***"     \cr
#' # keep the good, and perfect probes: \cr
#' idx <- grep("Good|Perfect", fData(x.averaged)$ProbeQuality) \cr
#' x.averaged.hiqual <- x.averaged[idx,]      \cr
#'  \cr
#' # \cr
#' # collapse to 1 row per gene \cr
#' # (using 2 of MJC's libraries) \cr
#' library(microarrays) \cr
#' library(metaGSEA) \cr
#' x.averaged.genes <- collapse(x.averaged, T, F, FUN=var, "SymbolReannotated") \cr
#' x.averaged.hiqual.genes <- collapse(x.averaged.hiqual, T, F, FUN=var, "SymbolReannotated") \cr
#' 
#' \tabular{ll}{
#' Package: \tab lumidat\cr
#' Type: \tab Package\cr
#' Version: \tab 1.0.5\cr
#' Date: \tab 2012-07-31\cr
#' License: \tab GenePattern license\cr
#' LazyLoad: \tab yes\cr
#' }
#'
#' @seealso \code{\link{preprocess.illumina.idat}}, \code{\link{lumiR.idat}}, \code{\link[limma]{read.ilmn}}
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
#' \dontrun{example(preprocess.illumina.idat)}
#' \dontrun{example(lumiR.idat)}
#' \dontrun{example(read.ilmn.idat)}
NULL
