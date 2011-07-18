#' Import an Illumina manifest file
#'
#' Import an Illumina manifest file. Only manifest files in TXT format are supported, ie, BGX files can't be read.
#' These files can be downloaded from Illumina [1,2]
#'
#' @param f the path to the manifest file name
#' @param verbose logical if \dQuote{TRUE} print informative messages. Default \dQuote{FALSE}.
#' @return a list with 2 elements, probes and controls. Both elements are \dQuote{data.frame}'s with as
#' many columns as found in the manifest file
#' @export
#' @seealso 
#' [1] \url{http://www.switchtoi.com/annotationfiles.ilmn}
#' 
#' [2] \url{http://www.switchtoi.com/annotationprevfiles.ilmn}
#' @author Mark Cowley, 2011-06-02
#' @examples
#' manifestfile <- system.file("extdata", "HumanHT-12_V3_0_R1_99999999.txt", package="lumidat")
#' manifest <- read.illumina.manifest(manifestfile, TRUE)
#' head(manifest$probes)
#' head(manifest$controls)
#'
read.illumina.manifest <- function(f, verbose=FALSE) {
	!missing(f) || stop("The path to a manifest file must be supplied")
	!grepl("bgx$", f, ignore.case=TRUE) || stop("BGX manifest files are not supported")
	file.exists(f) || stop("Can't find Illumina manifest file")
	if (verbose) cat(sprintf("Importing data from Illumina manifest file: %s\n", basename(f)))
	
	# Read the [Heading] section to determine the number of probes & controls to be imported.
	header <- readLines(f, 20)
	grepl("Illumina, Inc.", header[1]) || stop("File does not look like an Illumina manifest file")
	header[2] == "[Heading]" || stop("File does not look like an Illumina manifest file")
	
	numProbes <- grep("Number of Probes", header, value=TRUE)
	numProbes <- as.numeric(sub("Number of Probes[ \t]+", "", numProbes))
	numControls <- grep("Number of Controls", header, value=TRUE)
	numControls <- as.numeric(sub("Number of Controls[ \t]+", "", numControls))
	
	if( verbose ) cat(sprintf("Expecting %d probes, and %d control probes\n", numProbes, numControls))
	
	# import the probes.
	headerSkip <- which(header == "[Probes]")
	probes <- read.delim(f, skip=headerSkip, nrows=numProbes+1)
	if(nrow(probes) < numProbes) stop("too few probe rows found.")
	probes <- probes[1:numProbes, ]
	if( verbose ) cat(sprintf("Imported %d probes\n", numProbes))
	
	# import the controls.
	conSkip <- headerSkip + 1 + numProbes + 1 # first +1 for [Probes] column names; 2nd +1 for [Controls] line
	controls <- read.delim(f, skip=conSkip, nrows=numControls+1)
	if(nrow(controls) < numControls) stop("too few control probe rows found.")
	controls <- controls[1:numControls, ]
	if( verbose ) cat(sprintf("Imported %d control probes\n", numControls))
	
	# make result list.
	res <- list(probes=probes, controls=controls)
	res
}
