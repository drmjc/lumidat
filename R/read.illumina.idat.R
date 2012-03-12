#' Preprocess Illumina gene expression iDAT files.
#' 
#' This function can decrypt Illumina gene expression iDAT files (aka version 1 iDAT files). 
#' It will create \eqn{Sample Probe Profile.txt} and \eqn{Control Probe Profile.txt} files, 
#' similar to Illumina GenomeStudio version 1.8.0 [1]. We have made every effort to reproduce the 
#' GenomeStudio output, down to the Detection P-Value calculation, the order of the rows, the 
#' background correction procedure (which is only applied to gene-probes), the output file headers
#' sample and ProbeID naming, however we cannot guarantee that the output will be identical to
#' that produced by GenomeStudio.
#'
#' @section Array manifest files:
#' Array manifest files can be downloaded from Illumina [2,3]. It is important to
#'  use the TXT version, not the BGX version. You can use a newer release of the manifest file 
#' as long as you get the right array type. For example, HumanHT-12_V3 arrays can use the 
#' \eqn{HumanHT-12_V3_0_R2_11283641_A.txt} or \eqn{HumanHT-12_V3_0_R3_11283641_A.txt} 
#'  manifest files.
#' 
#' @section Probe naming:
#' There is some confusion about the naming of Illumina Probes. From the array manifest files, 
#' you can uniquely
#'  identify each probe by either the Illumina Probe ID (eg ILMN_1802252), the Array Address 
#' ID (eg 2640048), the Probe_Sequence, or the NuID [3].
#' Illumina's GenomeStudio uses the ArrayAddressID in the ProbeID column; 
#' I think the Probe_ID or NuID are better alternatives. 
#' If you choose \eqn{NuID}, note that unlike \code{\link[lumi]{addNuID2lumi}}, which uses
#' pre-built mapping tables, we calculate NuID's directly from the probe sequences within
#' the manifest file.
#' 
#' @section Probe collapsing:
#' Some genes are represented by multiple probes. You can either obtain values at the
#' probe-level, or the gene-level using the \sQuote{collapseMode} parameter:
#' \enumerate{
#' \item{\dQuote{none}}{The default: do no collapsing. Produces Sample Probe Profile.txt and 
#' Control Probe Profile.txt files}
#' \item{\dQuote{max}}{Choose the probe with the highest AVG_Intensity value}
#' \item{\dQuote{median}}{Choose the probe with the median AVG_Intensity value}
#' \item{\dQuote{mean}}{Determing the mean fluorescence across all probes}
#' }
#' Choosing either \dQuote{max}, \dQuote{median} or \dQuote{mean} creates 
#' \dQuote{Sample Gene Profile.txt} and \dQuote{Control Gene Profile.txt} files. 
#' \dQuote{mean} is the method that GenomeStudio uses, however we were unable to 
#' reproduce the mean values for \dQuote{numBeads}, \dQuote{Detection PValue} and 
#' \dQuote{Probe_STDERR} - we opted for taking the mean of these values. 
#' When using \dQuote{median} and there are an even number of probes for a gene, 
#' there is no single median probe, so from the 2 median probes, we take 
#' the conservative option of choosing the average AVG_Intensity, the smallest NumBeads, and the 
#' largest Probe_STDERR and Detection PVal, for those 2 median probes
#' 
#' @section Background correct:
#' Illumina GenomeStudio offers a background correction option. This estimates the background, and
#' subtracts it from only the gene-level probes; ie control probes are never background corrected.
#' The value is the mean AVG_Signal level of all the negative control probes on each array.
#' 
#' @section Memory usage:
#' At the heart of this function is a Java program, which requires that you specify an appropriate
#' maximum amount of memory. The default is -Xmx1024m, which reserves 1 GB RAM. We have analysed
#' 85 Human HT12 arrays using -Xmx2048m. If you get the following error:
#'     \sQuote{Exception in thread "main" java.lang.OutOfMemoryError: Java heap space}
#' Then you need to increase the amount of RAM, upto the maximum available in your system.
#' 
#' @section TODO:
#' represent manifest files the way that Affymetrix CDF's are (ie in data packages).
#' Import probe-level annotation from the manifest file into the resulting LumiBatch object.
#' Add options for background subtraction.
#' Add options for gene-level summarisation.
#' 
#' @param files A character vector of at least one file name
#' @param path The path to the directories where the files are. This defaults to the current
#'  working directory.
#' @param manifestfile The full path to the Array manifest file in TXT format.
#' @param probeID This controls which value to identify each probe by. Allowable values are 
#'    \dQuote{ArrayAddressID}, \dQuote{ProbeID}, \dQuote{Sequence}, \dQuote{NuID}.
#' @param collapseMode Collapse probes to genes using the specified mode. Valid values are 
#'    \dQuote{none} (the default), 
#'    \dQuote{max}, \dQuote{median} and \dQuote{mean} (the GenomeStudio default). 
#' @param backgroundCorrect logical, if TRUE then peform background correction on only the gene-level probes;
#'    if FALSE, do no correction.
#' @param outdir The full path to the output directory. If \eqn{NULL} the current working directory is used.
#' @param prefix An character[1] used as a file name prefix for the files that will be created.
#'    \dQuote{NONE} (the default) means no prefix will be used.
#' @param verbose logical, if TRUE, print informative messages
#' @param memory The maximum Java memory heapsize. Default -Xmx1024m, which allows 1GB of RAM.
#' @return If collapseMode=\dQuote{NONE}, then invisbly return a character[2] containing the 
#'   file paths of the \dQuote{Sample Probe Profile.txt} and \dQuote{Control Probe Profile.txt}
#'   files. Otherwise, invisbly return a character[2] containing the 
#'   file paths of the \dQuote{Sample Gene Profile.txt} and \dQuote{Control Gene Profile.txt}
#'   files.
#' @author Mark Cowley \email{m.cowley@@garvan.org.au}, with contributions from Mark Pinese, 
#'   David Eby.
#' @references
#' [1] Genome Studio Gene Expression Module User Guide, version 1.0, Illumina \cr
#' [2] \url{http://www.switchtoi.com/annotationfiles.ilmn} \cr
#' [3] \url{http://www.switchtoi.com/annotationprevfiles.ilmn} \cr
#' [4] Du, P., Kibbe, W. A., & Lin, S. M. (2007). nuID: a universal naming scheme of oligonucleotides
#' for illumina, affymetrix, and other microarrays Biol Direct, 2, 16. (doi:10.1186/1745-6150-2-16) \cr
#' @export
#' @examples
#' path <- system.file("extdata", package="lumidat")
#' outdir <- tempdir()
#' files <- c("5356583020_A_Grn.idat", "5356583020_B_Grn.idat")
#' manifestfile <- system.file("extdata", "HumanHT-12_V3_0_R1_99999999.txt", package="lumidat")
#' files <- read.illumina.idat(files, path, manifestfile, probeID="NuID", outdir=outdir)
#' files
#'
read.illumina.idat <- function(files=NULL, path=NULL, manifestfile=NULL, 
	probeID=c("ArrayAddressID", "ProbeID", "Sequence", "NuID"), 
	collapseMode=c("none", "max", "median", "mean"), 
	backgroundCorrect=FALSE,
	outdir=NULL, prefix=NULL, verbose=TRUE, memory="-Xmx1024m") {
	####################
	# check function arguments
	#
	
	# Check files + path
	files <- unique(files)
	length(files) > 0 || stop("Must supply at least 1 input file.")
	
	if( is.null(path) ) {
		if( !all(file.exists(files)) )
			stop("files must specify a vector of valid file names.")
	}
	else {
		files <- file.path(path, files)
		if( !all(file.exists(files)) )
			stop("files+path must specify a vector of valid file names.")
	}
	
	# check the manifestfile
	!is.null(manifestfile) || stop("Valid path to TXT version of an Illumina manifest file is required.")
	file.exists(manifestfile) || stop("Valid path to TXT version of an Illumina manifest file is required.")
	if( tolower(substring(manifestfile, nchar(manifestfile)-3, nchar(manifestfile))) == "bgx" )
		stop("You must use a TXT manifest file, not the BGX version.")
	
	# check probeID parameter.
	probeID <- match.arg(probeID[1], c("ArrayAddressID", "ProbeID", "Sequence", "NuID"), several.ok=FALSE)

	# check probeID parameter.
	collapseMode <- match.arg(collapseMode[1], c("none", "max", "median", "mean"), several.ok=FALSE)
	
	# check outdir
	if( is.null(outdir) ) outdir <- "." 
	
	# prefix
	if( is.null(prefix) ) prefixFlag <- ""
	else prefixFlag <- paste("-prefix", prefix)
	
	# backgroundCorrect
	backgroundCorrect <- ifelse(backgroundCorrect, "true", "false")
	
	# memory
	substring(memory, 1, 4) == "-Xmx" || stop("Invalid memory parameter. eg -Xmx2048m for 2GB RAM")
	####################
	
	####################
	# Setup the system call
	# Java stuff
	nzchar(Sys.which("java")) || stop("Can't find a valid Java")
	jar <- file.path(.path.package('lumidat'), 'bin', 'IlluminaGeneExpressionIdatReader-1.0.jar')
	file.exists(jar) || stop("Can't find jar.")
	
	# setup the command line
	quiet <- ifelse(verbose, "", "-quiet")
	flags <- sprintf("-outputDir %s -bg %s -collapse %s -manifestfile %s %s %s", 
		outdir, 
		backgroundCorrect,
		collapseMode,
		manifestfile,
		prefixFlag,
		quiet
	)
	files2 <- paste(shQuote(files), collapse=" ")
	cmd <- paste("java", memory, "-jar", jar, flags, files2)
	if(verbose) cat(cmd, "\n")
	####################
	
	####################
	# run the Java command
	system(cmd, intern=FALSE)
	####################

	####################
	# If it worked, then these files should have been created:
	filename.spp <- file.path(outdir, "Sample Probe Profile.txt")
	file.exists(filename.spp) || stop("Can't find the Sample Probe Profile.txt file created by running GenomeStudioFileCreator-1.0.jar")
	
	filename.cpp <- file.path(outdir, "Control Probe Profile.txt")
	file.exists(filename.cpp) || stop("Can't find the Control Probe Profile.txt file created by running GenomeStudioFileCreator-1.0.jar")
	####################

	res  <- c(filename.spp, filename.cpp)
	invisible(res)
}