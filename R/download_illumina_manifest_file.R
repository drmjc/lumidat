#' Download an Illumina array manifest
#'
#' Illumina array manifest files (current and archived) are hosted in [1,2], and exist
#' in 2 formats: bgx (compressed) and txt (most commonly used). Given the
#' exact name of the array, this will download either the bgx or txt formatted manifest file
#' to dir.
#' 
#' @section manifest files:
#' These are the current (2012-10-14) array manifest names:
#' \tabular{l}{
#'  HumanHT-12_V3_0_R2_11283641_A \cr       
#'  HumanHT-12_V3_0_R3_11283641_A \cr       
#'  HumanHT-12_V4_0_R1_15002873_B \cr       
#'  HumanHT-12_V4_0_R2_15002873_B \cr       
#'  HumanHT-12_V4_0_R2_15002873_B_WGDASL \cr
#'  HumanMI_V1_R2_XS0000122-MAP \cr         
#'  HumanMI_V2_R0_XS0000124-MAP \cr         
#'  HumanRef-8_V2_0_R4_11223162_A \cr       
#'  HUMANREF-8_V3_0_R1_11282963_A_WGDASL \cr
#'  HumanRef-8_V3_0_R2_11282963_A \cr       
#'  HumanRef-8_V3_0_R3_11282963_A \cr       
#'  HumanWG-6_V2_0_R4_11223189_A \cr        
#'  HumanWG-6_V3_0_R2_11282955_A \cr        
#'  HumanWG-6_V3_0_R3_11282955_A \cr        
#'  MouseMI_V1_R2_XS0000127-MAP \cr         
#'  MouseMI_V2_R0_XS0000129-MAP \cr         
#'  MouseRef-8_V1_1_R4_11234312_A \cr       
#'  MouseRef-8_V2_0_R2_11278551_A \cr       
#'  MouseRef-8_V2_0_R3_11278551_A \cr       
#'  MouseWG-6_V1_1_R4_11234304_A \cr        
#'  MouseWG-6_V2_0_R2_11278593_A \cr        
#'  MouseWG-6_V2_0_R3_11278593_A \cr        
#'  RatRef-12_V1_0_R5_11222119_A \cr 
#' }
#' For the archived ones, see \code{\link{list_illumina_manifest_files}}.
#' 
#' @param array.name the \emph{exact} array name, including capitalisation, as it appears in [1,2] or \code{\link{list_illumina_manifest_files}}.
#' @param type one of \dQuote{txt} or \dQuote{bgx}
#' @param dir the directory to download to. defaults to current working dir.
#' @param verbose logical
#' 
#' @return a character vector(1) of the path to the downloaded manifest file
#' @author Mark Cowley, 2012-10-14
#' @export
#' @seealso \code{\link{list_illumina_manifest_files}}
#' 
#' @examples
#' \dontrun{
#' download_illumina_manifest_file("HumanHT-12_V4_0_R2_15002873_B", "txt")
#' }
download_illumina_manifest_file <- function(array.name, type=c("txt", "bgx"), dir=".", verbose=FALSE) {
	type <- match.arg(type)
	is.dir <- function(path) {
	    return( !is.na(file.info(path)$isdir) & file.info(path)$isdir )
	}
	is.dir(dir) || stop("dir must exist, and must be a directory name")
	
	manifest.table <- .parse_illumina_manifest_html(verbose)
	array.name %in% manifest.table$Name || stop(array.name, "not found in the manifest table. see ?list_illumina_manifest_files")
	manifest.table <- subset(manifest.table, manifest.table$Name == array.name, drop=FALSE)
	nrow(manifest.table) == 1 || stop("Multiple URL's found for", array.name)
	
	url <- if(type == "txt") manifest.table$TXT else manifest.table$BGX
	destfile.zip <- file.path(dir, basename(url))
	download.file(url, destfile.zip) == 0 || stop("download of manifest file failed:", url)
	destfile <- unzip(destfile.zip, exdir=dir)
	unlink(destfile.zip)
	
	destfile
}


# parse out the array name, BGX url, TXT url from current and archived URL's & create a data.frame
# .parse_illumina_manifest_html(TRUE)
.parse_illumina_manifest_html <- function(verbose=FALSE) {
	url.current <- "http://www.switchtoi.com/annotationfiles.ilmn"
	url.previous <- "http://www.switchtoi.com/annotationprevfiles.ilmn"
	
	res.current <- .parse_illumina_manifest_html_array_url(url.current, verbose=verbose)
	res.current$current <- TRUE
	res.prev <- .parse_illumina_manifest_html_array_url(url.previous, verbose=verbose)
	res.prev$current <- FALSE
	
	res <- rbind(res.current, res.prev)
	res <- res[order(res$Name), ]
	rownames(res) <- 1:nrow(res)
	
	res
}


# parse out the array name, BGX url, TXT url from a single URL & create a data.frame
# url.current <- "http://www.switchtoi.com/annotationfiles.ilmn"
# url.previous <- "http://www.switchtoi.com/annotationprevfiles.ilmn"
# .parse_illumina_manifest_html_array_url(url.current)
.parse_illumina_manifest_html_array_url <- function(url, verbose=FALSE) {
	capabilities("http/ftp") || stop("http capabilities are required")
	stopifnot(length(url) == 1)
	
	f <- tempfile(pattern=basename(url))
	download.file(url, f, quiet=!verbose) == 0 || stop("failed to download URL")
	a <- readLines(f)
	a <- grep("Text Version", a, value=T)

	names <- sub(".*>([^><]+).bgx</a>.*", "\\1", a)

	.parse <- function(x) {
		x <- sub(".*href=..", "", x)
		x <- sub(".zip.*", ".zip", x)
		x <- paste("http://www.switchtoi.com/", x, sep="")
		x
	}
	bgx.url <- .parse(sub("\\|.*", "", a))
	txt.url <- .parse(sub(".*\\|", "", a))
	
	res <- data.frame(Name=names, BGX=bgx.url, TXT=txt.url, stringsAsFactors=FALSE)
	res
}
