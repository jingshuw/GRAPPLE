#' Read GWAS summary statistics data
#'
#' @importFrom tools file_ext
#' @importFrom data.table fread
#'
read_gwas_summary <- function(file, message = TRUE) {

    if (message) {
        message(paste("Reading", file))
    }

    if (file_ext(file) %in% c("rda", ".rData")) {
        load(file)
    } else {
        dat <- fread(file)
    }
    dat
}

#' Get intersection of the SNPs in GWAS summary datasets
#'
#' @import data.table
#'
#' @details The first file is the selection dataset.
#'
get_snp_intersection <- function(files, message = TRUE) {

    sel_dat <- read_gwas_summary(files[1], message)
    snp_intersection <- unique(sel_dat$SNP)

    if (length(files) > 1) {
        for (file in files[-1]) {

            dat_new <- read_gwas_summary(file, message)
            snp_intersection <- intersect(snp_intersection, dat_new$SNP)
        }
    }

    output <- sel_dat[SNP %in% snp_intersection, c("SNP", "pval")]
    output[!duplicated(SNP), ]

}

#' Extract selected SNPs from data
#'
#' @examples
#' files <- file.path("../data", list.files("../data"))
#' tmp <- get_snp_intersection(files)
#' tmp2 <- plink_clump(tmp, "../util/plink_mac/plink", refdat = "../util/data_maf0.01_rs")
#' dat <- extract_data(files, tmp2$SNP)
#' dat <- harmonise_data_list(dat)
#'
extract_data <- function(files, SNP_list = NULL, message = TRUE) {

    dat <- list()
    for (i in 1:length(files)) {
        dat[[i]] <- read_gwas_summary(files[i], message = message)
        dat[[i]] <- dat[[i]][SNP %in% SNP_list]
        dat[[i]] <- dat[[i]][!duplicated(SNP), ]
    }
    dat
}

#' Harmonise a list of GWAS summary data tables
#'
#' @importFrom TwoSampleMR harmonise_data
#'
harmonise_data_list <- function(dat) {

    dat_new <- list()
    SNP_list <- dat[[1]]$SNP
    dat_new[[1]] <- dat[[1]]

    ## Reference dataset
    data1 <- formatData(dat[[1]], "exposure")

    stopifnot(length(dat) > 1)

    for (i in 2:length(dat)) {
        ## Format and harmonise with the reference dataset
        data2 <- formatData(dat[[i]], "outcome")
        data2 <- suppressMessages(harmonise_data(data1, data2))

        ## Extract the columns corresponding to the new dataset
        data2 <- data2[data2$mr_keep, ]
        data2$SNP <- as.character(data2$SNP)
        data2 <- data2[, c(1, grep(".outcome", colnames(data2)))]
        colnames(data2) <- gsub(".outcome", "", colnames(data2))

        dat_new[[i]] <- data.table(data2)
        SNP_list <- intersect(SNP_list, data2$SNP)
    }

    dat_new <- lapply(dat_new, function(x) x[SNP %in% SNP_list][!duplicated(SNP), ])

    ## Check the data tables in the list have the same SNP, effect_allele, and reference_allele
    check_same <- function(x) {
        prod(sapply(x, all.equal, x[[1]])) == 1
    }
    stopifnot(check_same(lapply(dat_new, function(x) x$SNP)))
    stopifnot(check_same(lapply(dat_new, function(x) x$effect_allele)))
    stopifnot(check_same(lapply(dat_new, function(x) x$reference_allele)))

    dat_new
}


## Left to be done:
## Generate Jingshu's list of (data, marker.data, cor.mat) using the functions above.
