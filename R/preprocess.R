#' Preprocess GWAS summary statistics datasets
#'
#' @description
#' This function has GWAS summary statistics data files as inputs, perform genetic instrument selection and return matrices that are ready to use for GRAPPLE
#'
#'
#' @param sel.files A vector of the GWAS summary statistics file names for the risk factors SNP selection. Each GWAS file is a ".csv" or ".txt" file containing a data frame that at least has a column "SNP" for the SNP ids and "pval" for the p-values. The length of \code{sel.files} are not required to be the same as that of \code{exp.files} and the order of the files do not matter, while we strongly suggest having one selection file for each risk factor.
#' @param exp.files A vector of length \code{k} of the GWAS summary statistics file names of the \code{k} risk factors for getting the effect sizes and standard deviations. Each GWAS file should have a column "SNP" for the SNP ids, "beta" for the effect sizes, "se" for the standard deviation, "effect_allele" for the effect allele and "other_allele" for the other allele of the SNP. 
#' @param out.files The GWAS summary statistics file name for the disease data, can be a vector of length \code{m} to allow preprocessing \code{m} diseases simultaneously. Each GWAS file should have a column "SNP" for the SNP ids, "beta" for the effect sizes, "se" for the standard deviation, "effect_allele" for the effect allele and "other_allele" for the other allele of the SNP.
#' @param plink_refdat The reference genotype files (.bed, .bim, .fam) for clumping using PLINK (loaded with --bfile). 
#' @param max.p.thres The upper threshold of the selection p-values for a SNP to be selected before clumping. It only requires that at least one of the p-values of the risk factors of the SNPs to be below the threshold. Default is \code{0.01}.
#' @param cal.cor Whether calculate the \code{(k + 1)} by \code{(k + 1)} correlation matrix between the \code{k} risk factors and the outcome. The default is TRUE
#' @param p.thres.cor The lower threshold of the p-values for a SNP to be used in calculating the correlation matrix. It only select SNPs whose p-values are above the threshold for all risk factors. Default is \code{0.5}.
#' @param get.marker.candidates Whether getting SNPs which are used for mode marker selection. Only applies to cases where the number of risk factors \code{k = 1}. Default is TRUE for \code{k = 1}.
#' @param marker.p.thres P-value threshold of p-values in the exposure files for mode markers. Default is \code{1e-5}.
#' @param marker.p.source source of p-values of mode markers, a string of either "exposure" or "selection". Default is "exposure" for obtaining more markers. 
#' @param clump_r2 The clumping r2 threshold in PLINK for genetic instrument selection. Default is set to 0.001 for selection of independent SNPs. 
#' @param clump_r2_formarkers The clumping r2 threshold in PLINK. Default is set to 0.05 for selection of candidates for the marker SNPs. 
#' 
#' @return A list of selected summary statistics, which include
#' \item{data}{A data frame of size \code{p * (3 + 2k + 2m + 1)} for the effect sizes of \code{p} number of selected independent SNPs (instruments) on \code{k} risk factors (exposures). 
#' The first three columns include the SNP rsID, the effect allele and other allele after harmonizing, 
#' the next \code{2k} columns are the estimated effect sizes and standard deviations for the \code{k} risk factors stored in \code{exp.files}, 
#' the next \code{2m} columns are the estimated effect sizes and standard deviations for the \code{m} diseases stored in \code{exp.files}
#' and the the last columns are the selection p-values obtained from \code{sel.files}}
#' \item{marker.data}{A data frame for marker candidate SNPs, which has the same columns as \code{data}}. 
#' \item{cor.mat}{The estimated \code{(k + m)} by \code{(k + m)} correlation matrix between the \code{k} risk factors and the disease (outcome) shared by SNPs. The last column is for the outcome trait.}
#'

#' @importFrom data.table fread
#' @export
getInput <- function(sel.files,
                     exp.files,
                     out.files, 
                     plink_refdat, 
                     max.p.thres = 0.01, 
					           cal.cor = T, p.thres.cor = 0.5, 
					           get.marker.candidates = T,
					           marker.p.thres = 1e-5,
					           marker.p.source = "exposure",
                     clump_r2 = 0.001, 
					           clump_r2_formarkers = 0.05) {
	if (length(exp.files) > 1) {
		if (get.marker.candidates)
			print("Marker candidates will not be obtained as number of risk factors k > 1")
		get.marker.candidates <- F
	}

	if (missing(plink_refdat))
		stop("Missing PLINK reference files.")

	
	sel.SNPs <- NULL
	pvals <- NULL
	sel.SNPs.cor <- c()
	k <- length(exp.files)
	for (file in sel.files) {
		print(paste("loading data for selection:", file, "..."))
		file.type <- strsplit(file, "[.]")[[1]][2]
		if (file.type == "rda" || file.type == ".rData")
			load(file)
		else 
			dat <- data.frame(fread(file))
		sel.snps <- dat$SNP[dat$pval < max.p.thres]

		## keep SNPs with large p-values if we need to calculate the correlation of summary statistics for overlapping cohorts
		if (cal.cor) {
			sel.snps.cor <- dat$SNP[dat$pval > p.thres.cor]
		}

		pvals.tmp <- dat$pval[dat$pval < max.p.thres]
		names(pvals.tmp) <- make.names(sel.snps, unique = T)

		if (is.null(sel.SNPs)) {
			sel.SNPs <- sel.snps
			if (cal.cor)
				sel.SNPs.cor <- sel.snps.cor
			
			pvals <- pvals.tmp
		} else {
			ori.snps <- sel.SNPs
			sel.SNPs <- union(sel.SNPs, sel.snps)
			temp <- rep(1, length(sel.SNPs))
			names(temp) <- sel.SNPs
			temp[ori.snps] <- pvals
			temp1 <- rep(1, length(sel.SNPs))
			names(temp1) <- sel.SNPs
			temp1[sel.snps] <- pvals.tmp
			## when k > 1, take the minium of k risk factors' p-value as selection p-value
			pvals <- pmin(temp, temp1)
			names(pvals) <- sel.SNPs
			rm(temp, temp1)

			if (cal.cor)
				sel.SNPs.cor <- intersect(sel.SNPs.cor, sel.snps.cor)
		}
	}

  ## Bonferroni correction on the selected p-values when there are more than one risk factor
	pvals <- pvals * length(sel.files)


	beta_exp <- NULL
	se_exp <- NULL
	ref.data.exp <- NULL
	marker.SNPs <- NULL
	marker.pvals <- NULL
	for (exp.file in exp.files) {
		print(paste("loading data from exposure:", exp.file, "..."))

		file.type <- strsplit(exp.file, "[.]")[[1]][2]
		if (file.type == "rda" || file.type == ".rData")
			load(exp.file)
		else 
			dat <- data.frame(fread(exp.file))

		sel.SNPs <- intersect(sel.SNPs, dat$SNP)
		if (cal.cor)
			sel.SNPs.cor <- intersect(sel.SNPs.cor, dat$SNP)


		if (get.marker.candidates && marker.p.source == "exposure") {
			marker.SNPs <- dat$SNP[dat$pval < marker.p.thres]
			marker.pvals <- dat$pval[dat$pval < marker.p.thres]
			names(marker.pvals) <- make.names(marker.SNPs, unique = T)
		} else {
			marker.SNPs <- c()
		}
		dat <- dat[dat$SNP %in% c(union(sel.SNPs, marker.SNPs), 
								   sel.SNPs.cor), ]


    ## harmonize one dataset by one dataset
		if (!is.null(ref.data.exp)) {
			dat <- formatData(dat, "outcome")
			ref.data.exp <- suppressMessages(harmonise_data(ref.data.exp, dat))
			ref.data.exp <- ref.data.exp[ref.data.exp$mr_keep, ]
	    	ref.data.exp$SNP <- as.character(ref.data.exp$SNP)
			dat <- ref.data.exp[, c(1, grep(".outcome", colnames(ref.data.exp)))]
			colnames(dat) <- gsub(".outcome", "", colnames(dat))
			ref.data.exp <- ref.data.exp[, c(1, grep(".exposure", colnames(ref.data.exp)))]
			SNPs.kept <- unique(ref.data.exp$SNP)
			beta_exp <- beta_exp[SNPs.kept, , drop = F]
			se_exp <- se_exp[SNPs.kept, , drop = F]
			rownames(ref.data.exp) <- make.names(ref.data.exp$SNP, unique = T)
			ref.data.exp <- ref.data.exp[SNPs.kept, , drop = F]
			rownames(dat) <- make.names(dat$SNP, unique = T)
			dat <- dat[SNPs.kept, , drop = F]
			flip <- rep(1, nrow(beta_exp))
			flip[sign(beta_exp[, 1]) != sign(ref.data.exp$beta.exposure)] <- -1
			beta_exp <- flip * beta_exp
			sel.SNPs <- intersect(sel.SNPs, SNPs.kept)
			if (cal.cor)
				sel.SNPs.cor <- intersect(sel.SNPs.cor, SNPs.kept)
		} else {
			ref.data.exp <- formatData(dat, "exposure")
			rownames(ref.data.exp) <- make.names(ref.data.exp$SNP, unique = T)
		}

		if (is.null(beta_exp)) {
			beta_exp <- data.frame(dat$beta)
			rownames(beta_exp) <- make.names(dat$SNP, unique = T)
			se_exp <- data.frame(dat$se)
			rownames(se_exp) <- make.names(dat$SNP, unique = T)
		} else {
			beta_exp <- cbind(beta_exp, dat$beta)
			se_exp <- cbind(se_exp, dat$se)
		}
	}


	## load disease (outcome) file
    beta_out <- NULL
    se_out <- NULL
    for (out.file in out.files) {
        print(paste("loading data from outcome:", out.file, "..."))

        file.type <- strsplit(out.file, "[.]")[[1]][2]
        if (file.type == "rda" || file.type == ".rData")
            load(out.file)
        else
            dat <- data.frame(fread(out.file))                                                         

        dat <- dat[dat$SNP %in% c(union(sel.SNPs, marker.SNPs), 
                                   sel.SNPs.cor), ]


        ## harmonize one dataset by one dataset
    		dat <- formatData(dat, "outcome")
    		ref.data.exp <- suppressMessages(harmonise_data(ref.data.exp, dat))
    		ref.data.exp <- ref.data.exp[ref.data.exp$mr_keep, ]
    		ref.data.exp$SNP <- as.character(ref.data.exp$SNP)
    		dat <- ref.data.exp[, c(1, grep(".outcome", colnames(ref.data.exp)))]
    		colnames(dat) <- gsub(".outcome", "", colnames(dat))
    		ref.data.exp <- ref.data.exp[, c(1, grep(".exposure", colnames(ref.data.exp)))]
    		SNPs.kept <- unique(ref.data.exp$SNP)
    		beta_exp <- beta_exp[SNPs.kept, , drop = F]
    		se_exp <- se_exp[SNPs.kept, , drop = F]
    		rownames(ref.data.exp) <- make.names(ref.data.exp$SNP, unique = T)
    		ref.data.exp <- ref.data.exp[SNPs.kept, , drop = F]
    		rownames(dat) <- make.names(dat$SNP, unique = T)                                           
    		dat <- dat[SNPs.kept, , drop = F]
    		flip <- rep(1, nrow(beta_exp))
    		flip[sign(beta_exp[, 1]) != sign(ref.data.exp$beta.exposure)] <- -1
    		#    print(sum(flip == -1))
    		beta_exp <- flip * beta_exp
    
    		sel.SNPs <- intersect(sel.SNPs, SNPs.kept)
    		
    		if (cal.cor)
    			sel.SNPs.cor <- intersect(sel.SNPs.cor, SNPs.kept)
    		if (get.marker.candidates) {
    			marker.SNPs <- intersect(marker.SNPs, SNPs.kept)
		}


        if (is.null(beta_out)) {
            beta_out <- data.frame(dat$beta)
            rownames(beta_out) <- rownames(dat)
            se_out <- data.frame(dat$se)
            rownames(se_out) <- rownames(dat)
        } else {
      			beta_out <- beta_out[SNPs.kept, , drop = F]
      			se_out <- se_out[SNPs.kept, , drop = F]
            beta_out <- cbind(beta_out, dat$beta)
            se_out <- cbind(se_out, dat$se)
        }
    }

	pvals <- pvals[sel.SNPs]
	if (get.marker.candidates) {
		if (marker.p.source == "exposure")
			marker.pvals <- marker.pvals[marker.SNPs]
		else {
			marker.SNPs <- sel.SNPs[pvals < marker.p.thres]
			marker.pvals <- pvals[pvals < marker.p.thres]
		}
	}

	dat <- dat[dat$SNP %in% c(sel.SNPs, marker.SNPs), , drop = F]


	if (cal.cor) {
		zz <- cbind(beta_exp, se_exp, beta_out, se_out)[sel.SNPs.cor,]
		z.values <- cbind(beta_exp[sel.SNPs.cor, ]/se_exp[sel.SNPs.cor, ],
						  beta_out[sel.SNPs.cor, ]/se_out[sel.SNPs.cor, ])
		colnames(z.values) <- c(paste0("exposure", 1:length(exp.files)), 
								paste0("outcome", 1:length(out.files)))

		z.values <- as.matrix(z.values)
		z.values <- z.values[rowSums(is.na(z.values)) == 0,  , drop = F]
		covv <- t(z.values) %*% z.values / nrow(z.values)
		varr <- colMeans(z.values^2, na.rm = T)
		corr <- t(covv / sqrt(varr))/sqrt(varr)
	} else
		corr <- NULL


	data.sel <- data.frame(SNP = sel.SNPs, pval = pvals)
 	
	print("Start clumping using PLINK ...")
	data.sel <- plink_clump(data.sel, "plink", 
							plink_refdat, clump_r2 = clump_r2)
	sel.SNPs <- as.character(data.sel$SNP)
	

	if (get.marker.candidates) {
		data.sel <- data.frame(SNP = marker.SNPs, pval = marker.pvals)
		data.sel <- plink_clump(data.sel, "plink", 
								plink_refdat, 
								clump_r2 = clump_r2_formarkers)
		marker.SNPs <- as.character(data.sel$SNP)
	}

	colnames(beta_exp) <- paste0("gamma_exp", 1:length(exp.files))
	colnames(se_exp) <-  paste0("se_exp", 1:length(exp.files))

	colnames(beta_out) <- paste0("gamma_out", 1:length(out.files))
	colnames(se_out) <-  paste0("se_out", 1:length(exp.files))

	
	beta_exp.marker <- beta_exp[as.character(marker.SNPs), , drop = F]
	se_exp.marker <- se_exp[as.character(marker.SNPs), , drop = F]
	beta_out.marker <- beta_out[as.character(marker.SNPs), , drop = F]
	se_out.marker <- se_out[as.character(marker.SNPs), , drop = F]
	meta_data.marker <- dat[as.character(marker.SNPs), c("SNP", "effect_allele", 
															  "other_allele")]

	beta_exp <- beta_exp[as.character(sel.SNPs), , drop = F]
	se_exp <- se_exp[as.character(sel.SNPs), , drop = F]
	beta_out <- beta_out[as.character(sel.SNPs), , drop = F]
	se_out <- se_out[as.character(sel.SNPs), , drop = F]
	meta_data <- dat[as.character(sel.SNPs), c("SNP", "effect_allele", 
													"other_allele")]

	print(paste(nrow(beta_exp), "independent genetic instruments extracted. Done!"))

	pvals <- pvals[sel.SNPs]
	names(pvals) <- sel.SNPs


	marker.pvals <- marker.pvals[marker.SNPs]
	names(marker.pvals) <- marker.SNPs

	data <- cbind(meta_data, beta_exp, se_exp, beta_out, se_out, selection_pvals = pvals)
	marker.data <- cbind(meta_data.marker, beta_exp.marker, se_exp.marker,
						 beta_out.marker, se_out.marker, selection_pvals = marker.pvals)


	gc(full = F)

	return(list(data = data, marker.data = marker.data,
				cor.mat = corr))
}





