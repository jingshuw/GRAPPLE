#' Preprocess GWAS summary statistics datasets
#'
#' @description
#' This function has GWAS summary statistics data files as inputs, perform genetic instrument selection and return matrices that are ready to use for GRAPPLE
#'
#'
#' @param sel.files A vector of length \code{k} of the GWAS summary statistics file names of the \code{k} risk factors SNP selection. Each GWAS file is a ".csv" or ".txt" file containing a data frame that at least has a column "SNP" for the SNP ids, "pval" for the p-values. 
#' @param exp.files A vector of length \code{k} of the GWAS summary statistics file names of the \code{k} risk factors for getting the effect sizes and standard deviations. Each GWAS file should have a column "SNP" for the SNP ids, "beta" for the effect sizes, "se" for the standard deviation, "effect_allele" for the effect allele of the SNP (capitalized letters) and "other_allele" for the other allele of the SNP (capitalized letters). 
#' @param out.file The GWAS summary statistics file name for the outcome data. Each GWAS file should have a column "SNP" for the SNP ids, "beta" for the effect sizes, "se" for the standard deviation, "effect_allele" for the effect allele of the SNP (capitalized letters) and "other_allele" for the other allele of the SNP (capitalized letters).
#' @param plink_exe The executable file of PLINK. PLINK should be first downloaded from \url{https://www.cog-genomics.org/plink2}. 
#' @param plink_refdat The reference genotype files (.bed, .bim, .fam) for clumping using PLINK (loaded with --bfile). 
#' @param max.p.thres The upper threshold of the selection p-values for a SNP to be selected before clumping. It only requires that at least one of the p-values of the risk factors of the SNPs to be below the threshold. Default is \code{0.01}.
#' @param cal.cor Whether calculate the \code{(k + 1)} by \code{(k + 1)} correlation matrix between the \code{k} risk factors and the outcome. The default is TRUE
#' @param p.thres.cor The lower threshold of the p-values for a SNP to be used in calculating the correlation matrix
#' @param get.marker.candidates Whether getting SNPs which are used for mode marker selection. Only applies to cases where the number of risk factors \code{k = 1}. Default is TRUE for \code{k = 1}.
#' @param marker.p.thres Selection threshold of p-values for mode markers. Default is \code{1e-5}.
#' @param clump_r2 The clumping r2 threshold in PLINK for genetic instrument selection. Default is set to 0.001 for selection of independent SNPs. 
#' @param clump_r2_formarkers The clumping r2 threshold in PLINK. Default is set to 0.05 for selection of candidates for the marker SNPs. 
#' 
#' @return A list of selected summary statistics, which includes
#' \item{b_exp}{A matrix of size \code{p * k} for the effect sizes of \code{p} number of selected independent SNPs (instruments) on \code{k} risk factors (exposures)}
#' \item{se_exp}{A matrix of size \code{p * k} for the standard deviations of each entry in \code{b_exp}}
#' \item{b_out}{A vector of size \code{p} for the effect sizes of \code{p} number of selected independent SNPs (instruments) on the disease (outcome)}
#' \item{se_out}{A vector of size \code{p} for the standard deviations of each entry in \code{b_out} 
#' \item{meta_data}{A data frame for the information of selected SNPs containing 3 columns, the SNP rsID, the effect allele and other allele after harmonizing (where the \code{b_exp} and \code{b_out} are the effects of the effect_alleles).}
#' \item{sel.pvals}{A vector of size \code{p} for p-values of selected SNPs obtained from the selection files. When there are \code{k > 1} selection files, each p-value is the Bonferronni combination of \code{k} p-values.}
#' \item{marker.data}{A list of data for marker candidate SNPs, whose elements are \code{b_exp}, \code{se_exp}, \code{b_out}, \code{se_out}, \code{meta_data} and \code{sel.pvals}.}
#' \item{cor.mat}{The estimated \code{(k + 1)} by \code{(k + 1)} correlation matrix between the \code{k} risk factors and the disease (outcome) shared by SNPs. The last column is for the outcome trait.}
#'

#' @importFrom data.table fread
#' @export
getInput <- function(sel.files,
                     exp.files,
                     out.file, 
                     plink_exe, 
					 plink_refdat, 
                     max.p.thres = 0.01, 
					 cal.cor = T, p.thres.cor = 0.5, 
					 get.marker.candidates = T,
					 sample.split.marker.candidates = F,
					 marker.p.thres = 1e-5,
                     clump_r2 = 0.001, clump_r2_formarkers = 0.05) {
	if (length(exp.files) > 1)
		get.marker.candidates <- F

	if (missing(plink_exe) || missing(plink_refdat))
		stop("Missing PLINK files.")
	
	sel.SNPs <- NULL
	pvals <- NULL
	sel.SNPs.cor <- c()
	for (file in sel.files) {
		print(paste("loading data for selection:", file, "..."))
		file.type <- strsplit(file, "[.]")[[1]][2]
		if (file.type == "rda" || file.type == ".rData")
			load(file)
		else 
			dat <- data.frame(fread(file))
		sel.snps <- dat$SNP[dat$pval < max.p.thres]


		## keep SNPs with large p-values if we need to calculate \Sigma
		if (cal.cor) {
			sel.snps.cor <- dat$SNP[dat$pval > p.thres.cor]
		}


		if (is.null(sel.SNPs)) {
			sel.SNPs <- sel.snps
			if (cal.cor)
				sel.SNPs.cor <- sel.snps.cor

			pvals <- dat$pval[dat$pval < max.p.thres]
			names(pvals) <- sel.SNPs
		} else {
			ori.snps <- sel.SNPs
			sel.SNPs <- union(sel.SNPs, sel.snps)
			temp <- rep(1, length(sel.SNPs))
			names(temp) <- sel.SNPs
			temp[ori.snps] <- pvals
			temp1 <- rep(1, length(sel.SNPs))
			names(temp1) <- sel.SNPs
			temp1[sel.snps] <- dat$pval[dat$pval < max.p.thres]
			## when k > 1, take the minium of k risk factors' p-value as selection p-value
			pvals <- pmin(temp, temp1)
			#  pvals <- temp
			names(pvals) <- sel.SNPs
			rm(temp)

			if (cal.cor)
				sel.SNPs.cor <- intersect(sel.SNPs.cor, sel.snps.cor)
		}
		rm(dat)
	}

  ## Bonferroni correction of the minium p-value across k risk factors
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

		if (get.marker.candidates && !sample.split.marker.candidates) {
			marker.snps <- dat$SNP[dat$pval < marker.p.thres]
			if (is.null(marker.SNPs)) {
				marker.SNPs <- marker.snps
				marker.pvals <- dat$pval[dat$pval < marker.p.thres]
			} else {
				ori.snps <- marker.SNPs
				marker.SNPs <- union(marker.SNPs, marker.snps)
				temp <- rep(1, length(marker.SNPs))
				names(temp) <- marker.SNPs
				temp[ori.snps] <- marker.pvals
				temp1 <- rep(1, length(marker.SNPs))
				names(temp1) <- marker.SNPs
				temp1[marker.snps] <- dat$pval[dat$pval < marker.p.thres]
				## when k > 1, take the minium of k risk factors' p-value as selection p-value
				marker.pvals <- pmin(temp, temp1)
				#  pvals <- temp
			}
			names(marker.pvals) <- marker.SNPs
		} else
			marker.SNPs <- c()
		temp <- dat[dat$SNP %in% c(union(sel.SNPs, marker.SNPs), 
								   sel.SNPs.cor), ]
		rm(dat)


    ## harmonize one dataset by one dataset
		if (!is.null(ref.data.exp)) {
			temp <- formatData(temp, "outcome")
			temp.result <- suppressMessages(harmonise_data(ref.data.exp, temp))
			temp.result <- temp.result[temp.result$mr_keep, ]
			temp <- temp.result[, c(1, grep(".outcome", colnames(temp.result)))]
			colnames(temp) <- gsub(".outcome", "", colnames(temp))
			ref.data.exp <- temp.result[, c(1, grep(".exposure", colnames(temp.result)))]
			#    colnames(data.sel) <- gsub(".exposure", "", colnames(data.sel))
			rm(temp.result)
			SNPs.kept <- unique(as.character(ref.data.exp$SNP))
			beta_exp <- beta_exp[SNPs.kept, , drop = F]
			se_exp <- se_exp[SNPs.kept, , drop = F]
			flip <- rep(1, nrow(beta_exp))
			flip[sign(beta_exp[, 1]) != sign(ref.data.exp$beta.exposure)] <- -1
			#    print(sum(flip == -1))
			beta_exp <- flip * beta_exp
			sel.SNPs <- intersect(sel.SNPs, SNPs.kept)
			if (cal.cor)
				sel.SNPs.cor <- intersect(sel.SNPs.cor, SNPs.kept)
			if (get.marker.candidates && !sample.split.marker.candidates) 
				marker.SNPs <- intersect(marker.SNPs, SNPs.kept)
		} else {
			ref.data.exp <- formatData(temp, "exposure")
		}

		if (is.null(beta_exp)) {
			beta_exp <- data.frame(temp$beta)
			rownames(beta_exp) <- make.names(temp$SNP, unique = T)
			se_exp <- data.frame(temp$se)
			rownames(se_exp) <- make.names(temp$SNP, unique = T)
		} else {
			beta_exp <- beta_exp[as.character(temp$SNP), , drop = F]
			se_exp <- se_exp[as.character(temp$SNP), , drop = F]
			beta_exp <- cbind(beta_exp, temp$beta)
			se_exp <- cbind(se_exp, temp$se)
		}
		rm(temp)
	}


	## load disease (outcome) file
	print(paste("loading data for outcome:", out.file, "..."))
	file.type <- strsplit(out.file, "[.]")[[1]][2]
	if (file.type == "rda" || file.type == ".rData")
		load(out.file)
	else 
		dat <- data.frame(fread(out.file))

	sel.SNPs <- intersect(sel.SNPs, dat$SNP)
	if (cal.cor)
		sel.SNPs.cor <- intersect(sel.SNPs.cor, dat$SNP)
	if (get.marker.candidates && !sample.split.marker.candidates) {
		marker.SNPs <- intersect(marker.SNPs, dat$SNP)
	} 
	temp <- dat[dat$SNP %in% c(union(sel.SNPs, marker.SNPs), 
							   sel.SNPs.cor), ]
	rm(dat)



  ## harmonize
	temp <- formatData(temp, "outcome")
	temp.result <- suppressWarnings(suppressMessages(harmonise_data(ref.data.exp, temp)))
	temp.result <- temp.result[temp.result$mr_keep, ]
	temp <- temp.result[, c(1, grep(".outcome", colnames(temp.result)))]
	colnames(temp) <- gsub(".outcome", "", colnames(temp))
	ref.data.exp <- temp.result[, c(1, grep(".exposure", colnames(temp.result)))]
	#    colnames(data.sel) <- gsub(".exposure", "", colnames(data.sel))
	rm(temp.result)

	SNPs.kept <- unique(as.character(ref.data.exp$SNP))
	beta_exp <- beta_exp[SNPs.kept, , drop = F]
	se_exp <- se_exp[SNPs.kept, , drop = F]

	flip <- rep(1, nrow(beta_exp))
	flip[sign(beta_exp[, 1]) != sign(ref.data.exp$beta.exposure)] <- -1
	#  print(sum(flip == -1))

	rm(ref.data.exp)

	beta_exp <- flip * beta_exp
#	sel.SNPs <- rownames(beta_exp)

	sel.SNPs <- intersect(sel.SNPs, SNPs.kept)
	if (cal.cor)
		sel.SNPs.cor <- intersect(sel.SNPs.cor, SNPs.kept)
	if (get.marker.candidates && !sample.split.marker.candidates) 
		marker.SNPs <- intersect(marker.SNPs, SNPs.kept)

	data_out <- temp
	rownames(data_out) <- make.names(data_out$SNP, unique = T)
	pvals <- pvals[sel.SNPs]
	marker.pvals <- marker.pvals[marker.SNPs]


	if (cal.cor) {
		z.values <- cbind(beta_exp[sel.SNPs.cor, ]/se_exp[sel.SNPs.cor, ],
						  data_out[sel.SNPs.cor, ]$beta/data_out[sel.SNPs.cor, ]$se)
		colnames(z.values) <- c(paste0("exposure", 1:length(exp.files)), "outcome")

		z.values <- as.matrix(z.values)
		z.values <- z.values[rowSums(is.na(z.values)) == 0,  , drop = F]
		covv <- t(z.values) %*% z.values / nrow(z.values)
		varr <- colMeans(z.values^2, na.rm = T)
		corr <- t(covv / sqrt(varr))/sqrt(varr)

	#	beta_exp <- beta_exp[-sel.SNPs.cor, , drop = F]
	#	se_exp <- se_exp[-sel.SNPs.cor, , drop = F]
	} else
		corr <- NULL

	data.sel <- data.frame(SNP = sel.SNPs, pval = pvals)
 
	data.sel <- plink_clump(data.sel, plink_exe, 
							plink_refdat, clump_r2 = clump_r2)
	sel.SNPs <- as.character(data.sel$SNP)

	
	

	
	if (get.marker.candidates) {
		if (sample.split.marker.candidates) {
			marker.pvals <- pvals[pvals < marker.p.thres]
			marker.SNPs <- names(marker.pvals)
		} else
			marker.pvals <- marker.pvals * length(exp.files)

		data.sel <- data.frame(SNP = marker.SNPs, pval = marker.pvals)
		data.sel <- plink_clump(data.sel, plink_exe, 
								plink_refdat, 
								clump_r2 = clump_r2_formarkers)
		marker.SNPs <- as.character(data.sel$SNP)
	}

	colnames(beta_exp) <- paste0("Exposure", 1:length(exp.files))
	
	beta_exp.marker <- beta_exp[as.character(marker.SNPs), , drop = F]
	se_exp.marker <- se_exp[as.character(marker.SNPs), , drop = F]
	beta_out.marker <- data_out[as.character(marker.SNPs), , drop = F]$beta
	se_out.marker <- data_out[as.character(marker.SNPs), , drop = F]$se
	meta_data.marker <- data_out[as.character(marker.SNPs), c("SNP", "effect_allele", 
															  "other_allele")]

	beta_exp <- beta_exp[as.character(sel.SNPs), , drop = F]
	se_exp <- se_exp[as.character(sel.SNPs), , drop = F]
	beta_out <- data_out[as.character(sel.SNPs), , drop = F]$beta
	se_out <- data_out[as.character(sel.SNPs), , drop = F]$se
	meta_data <- data_out[as.character(sel.SNPs), c("SNP", "effect_allele", 
													"other_allele")]

	rm(data_out)

	pvals <- pvals[sel.SNPs]
	names(pvals) <- sel.SNPs


	marker.pvals <- marker.pvals[marker.SNPs]
	names(marker.pvals) <- marker.SNPs


	gc(full = F)

  	return(list(b_exp = beta_exp, se_exp = se_exp,
              b_out = beta_out, se_out = se_out,
			  meta_data = meta_data, sel.pvals = pvals,
			  marker.data = list(b_exp = beta_exp.marker,
								 se_exp = se_exp.marker,
								 b_out = beta_out.marker,
								 se_out = se_out.marker,
								 meta_data = meta_data.marker,
								 sel.pvals = marker.pvals),
			  cor.mat = corr))
}





