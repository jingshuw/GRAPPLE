#' Query NCBI's dbSNP for summary information on a set of SNPs. This function is modified from \code{rsnps} package to make it work with the NCBI data format change 
#'
#' @export
#' @param x A vector of SNPs (with or without 'rs' prefix)
#' @param key (character) NCBI Entrez API key. optional. 
#' See "NCBI Authenication" in [rsnps-package]
#' @param ... Curl options passed on to [crul::HttpClient]
#' @return data.frame with three columns: \code{snp_id}, \code{gene2} and \code{chrpos}. SNPs not found are omitted
#' @import httr
#' @import xml2
#' @importFrom crul HttpClient
#' @importFrom data.table data.table
ncbi_snp_summary <- function (x, key = NULL, ...) 
{
	httr::set_config(httr::config(http_version = 0))
	stopifnot(inherits(x, "character"))
    x <- gsub("^rs", "", x)
    key <- check_key(key %||% "")
    args <- rsnps_comp(list(db = "snp", retmode = "flt", rettype = "flt", 
        id = paste(x, collapse = ","), api_key = key))
    url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    cli <- HttpClient$new(url = url, opts = list(...))
    res <- cli$get(query = args)
    res$raise_for_status()
    tmp <- res$parse("UTF-8")
    xml <- xml2::read_xml(tmp)
#	return(xml)
	docsums <- xml2::xml_find_all(xml, "//DocumentSummary")
#	return(docsums)
    docsums <- xml2::xml_text(xml2::xml_find_all(xml, "//DOCSUM"))
	chrpos <- xml2::xml_text(xml2::xml_find_all(xml, "//CHRPOS"))
	snp_id <- xml2::xml_text(xml2::xml_find_all(xml, "//SNP_ID"))
	dats <- data.table(matrix(c(snp_id, chrpos, docsums), ncol = 3))#, stringsAsFactors = F) 
	colnames(dats) <- c("snp_id", "chrpos", "gene2")
	dats$gene2 <- sapply(dats$gene2, function(z) {
					   gn <- stract_(z, "GENE=[A-Za-z0-9-]+:[0-9]+,?([A-Za-z0-9-]+:[0-9]+)?")
					   gn <- sub("GENE=", "", gn)
					   z <- gn %||% NA_character_
					   z
	})
	return(dats)


#	print(dats)
#    rbl(dats)
}

#' Original check_key function from the rsnps R package
#' 
#' @keywords internal
#'
check_key <- function(key = "") {
  stopifnot(is.character(key))
  if (nzchar(key)) return(key)
  key <- Sys.getenv("ENTREZ_KEY")
  if (!nzchar(key)) NULL else key
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x


rsnps_comp <- function(x) Filter(Negate(is.null), x)

stract_ <- function(string, pattern) {
  regmatches(string, regexpr(pattern, string))
}

#rbl <- function(x) {
#  (xxxxx <- data.table::setDF(
#    data.table::rbindlist(x, use.names = TRUE, fill = TRUE)
#  ))
#}
#
#make_named_list <- function(x) {
#  as.list(stats::setNames(xml2::xml_text(x), 
#    tolower(xml2::xml_name(x))))
#}
#
