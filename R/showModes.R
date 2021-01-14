#' Use the multiple modes of the robust profile likelihood function to find out multiple pathways in MR and their marker SNPs.
#'
#' THis function can be run only for \code{k = 1} when there is only one risk factor
#'
#' @inheritParams grappleRobustEst
#' @param marker.data A data frame containing the information of candidate marker genes.  
#' Default is NULL, which sets \code{marker.data} to \code{data}. 
#' Another choice is to use the \code{marker.data} element in the output of \code{getInput}. 
#' @param marker.p.thres P-value threshold for marker SNP selection. See \code{p.thres} 
#' @param mode_lmts The range of \code{beta} that the modes are searched from. Default is \code{c(-5, 5)}
#' @param k.findmodes Tuning parameters of the loss function, for loss "l2", it is NA, for loss "huber", default is 1.345 and for loss "tukey", default is 3. 
#' @param include.thres Absolute value upper threshold of the standardized test statistics of one SNP on one mode for the SNP to be included as a marker for that mode, default is 1.4
#' @param exclude.thres Absolute value lower threshold of the standardized test statistics of one SNP on other modes for the SNP to be included as a marker for that mode, default is \code{qnorm(0.975)}
#' @param map.marker Whether map each marker to the earist gene or not. Default is TRUE if multiple markers are found. It is always FALSE if there is just one mode.
#' @param ldThres the parameter passed to the \code{queryhaploReg} function. Increase to 1 when there is a "timeout" error.
#' @param npoints Number of equally spaced points chosen for grid search of modes within the range \code{mode.lmts}.
#'
#' @return A list containing the following elements:
#' \item{fun}{The profile likelihood function with argument \code{beta}}
#' \item{modes}{The position of modes. Only include modes where marker genes can be detected}
#' \item{p}{The profile likelihood plot with gene markers when there are multiple modes. The range of the x.axis depends on the distance between the maximum mode and minimum mode when there are multiple modes.}
#' \item{markers}{A data frame of marker information}
#' \item{raw.modes}{All modes of the profile likelihood function within the range of \code{mode_lmts}}
#' \item{supp_gwas}{More information about the markers.}
#'
#' @import ggplot2 
#' @importFrom ggrepel geom_text_repel
#' @importFrom haploR queryHaploreg
#' @export
findModes <- function(data,
          					  p.thres = NULL, 
          					  marker.data = NULL,      
	  		        		  marker.p.thres = NULL,	  
				          	  mode.lmts = c(-5, 5),
                      cor.mat = NULL, 
                      loss.function = c("tukey", "huber", "l2"), 
                      k.findmodes = switch(loss.function[1], 
                                 l2 = NA, huber = 1.345, 
                                 tukey = 3), 
                      include.thres = 1, exclude.thres = 2,
                      map.marker = T, ldThres = 0.9, 
					  npoints = 10000) {

  nn <- colnames(data)
  if (!is.null(p.thres)) {
	  if (length(grep("selection_pvals", nn)) == 0)
		  stop("data need to include a column named selection_pvals for the selection p-value of each SNP.")
	  data <- data[data$selection_pvals < p.thres, ]
  }
  
	b_exp <- as.matrix(data[, grep("gamma_exp", nn), drop = F])
	se_exp <-  as.matrix(data[, grep("se_exp", nn), drop = F])
	b_out <-  data[, grep("gamma_out", nn)[1]]
	se_out <-  data[, grep("se_out", nn)[1]]

  if (is.null(cor.mat))
    cor.mat <- diag(rep(1, ncol(b_exp) + 1))

  if (is.null(marker.data)) 
	  marker.data <- data

	nn <- colnames(marker.data)
	marker.b_exp <- as.matrix(marker.data[, grep("gamma_exp", nn), drop = F])
	marker.se_exp <-  as.matrix(marker.data[, grep("se_exp", nn), drop = F])
	marker.b_out <-  marker.data[, grep("gamma_out", nn)[1]]
	marker.se_out <-  marker.data[, grep("se_out", nn)[1]]

  loss.function <- match.arg(loss.function, c("tukey", "huber", "l2"))

	t_fun_marker <- function(beta, tau2 = 0) {
			upper <-marker.b_out - as.matrix(marker.b_exp) %*% beta
			temp <- t(t(cbind(as.matrix(marker.se_exp), 
							  marker.se_out)) * c(-beta, 1))
			lower <- sqrt(rowSums((temp %*% cor.mat) * temp) + tau2)
			return(upper/lower)
	}

	robust.optfun.fixtau <- robustLossFixtau(data = NULL,
											 b_exp, b_out, se_exp, se_out,
											 cor.mat, loss.function, k.findmodes)



	## Take npoints equally spaced points to do grid search to check for modes
	beta.seq <- seq(mode.lmts[1], mode.lmts[2], length.out = npoints)
	val <- try(sapply(beta.seq, function(beta) robust.optfun.fixtau(beta, 0)))
	if (class(val) == "try-error")
		stop("Possible error: finding modes is currently available only for
			 univariate MR with only one risk factor!")

	mode.pos <- findLocalModes(val)
	temp.data <- data.frame(beta = beta.seq, likelihood = val)

	llk.limit <- range(val)

	beta.mode <- beta.seq[mode.pos]
	tmp.range <- max(beta.mode) - min(beta.mode)
	tmp.lines <- data.frame(modes = beta.mode,                                              
						    mod.col = as.factor(beta.mode))
	p <- ggplot2::ggplot(aes(x = beta, y = likelihood), data = temp.data) + geom_line() +
		geom_vline(xintercept = 0) +
		geom_vline(data = tmp.lines,
				   map = aes(xintercept = modes, color = mod.col), linetype = "dashed") +
		labs(y = "Profile lieklihood") +
		theme(axis.line = element_line(linetype = "solid"),
			  axis.text.y=element_blank(),
			  axis.ticks.y=element_blank(),
			  legend.position="none",
			  panel.background=element_blank(),
			  panel.border=element_blank(),
			  panel.grid.major=element_blank(),
			  panel.grid.minor=element_blank(),
			  plot.background=element_blank(),
			  plot.title=element_text(size = 11))



	res.mat <- sapply(beta.mode, function(beta) tt <- t_fun_marker(beta,0))

	if (length(beta.mode) == 1)
		res.mat <- as.matrix(res.mat)

	colnames(res.mat) <- beta.mode

	rownames(res.mat) <- rownames(marker.b_exp)

	ss <- which(rowSums(abs(res.mat) < exclude.thres) == 1 & 
				rowSums(abs(res.mat) < include.thres) > 0)

	res.mat <- res.mat[ss, , drop = F]
	markers <- as.data.frame(abs(res.mat) < include.thres)

	keep.mode <- colSums(markers) > 0
	markers <- markers[, keep.mode, drop = F]
	res.mat <- res.mat[, keep.mode, drop = F]


	if (ncol(markers) <= 1)
		map.marker <- F
	if (map.marker) {

	

		p <- p +  xlim(max(mode.lmts[1], min(beta.mode) -  2 * tmp.range), 
					   min(mode.lmts[2], max(beta.mode) + 2 * tmp.range)) 


		snp_ids <- rownames(markers)

  	## map to HaploReg
	
    results <- queryHaploreg(query = snp_ids, ldThres = ldThres)
	  results <- results[, c("rsID", "chr", "pos_hg38", "GENCODE_name", 
										 "gwas", "dbSNP_functional_annotation",
										 "is_query_snp", "r2")]
  	tt <- strsplit(results$gwas, split = ";")
  	tt <- sapply(tt, function(traits) {
  					 trait <- strsplit(traits, split = ",")
  					 if (trait[[1]][1] == ".")
  						 return(traits)
  					 trait_name <- sapply(trait, function(item) item[2])
  					 pvalues <- sapply(trait, function(item) as.numeric(item[3]))
  					 return(paste0(unique(trait_name[sort(pvalues, index.return = T)$ix]), collapse = ","))				 
  							 })
  	results$gwas_short <- tt
  	results.supp <- results[results$is_query_snp == 0 & results$gwas != ".", ]
  	results <- results[results$is_query_snp == 1, ]
  	ss <- subset(results, select = c(gwas, dbSNP_functional_annotation))
  	results <- subset(results, select = -c(gwas, dbSNP_functional_annotation, is_query_snp, r2))
  	results.supp <- subset(results.supp, select = -is_query_snp)
  	results.supp <- results.supp[order(results.supp$GENCODE_name), ]
  
  	markers <- cbind(markers, res.mat)
  	colnames(markers) <- c(paste0("Mode", 1:sum(keep.mode), "_marker"), 
  						   paste0("Mode", 1:sum(keep.mode), "_stats"))
      markers <- cbind(results, data.frame(markers)[results$rsID, , drop = F], ss)
  	res.mat <- data.frame(res.mat)[markers$rsID, , drop = F]

  } else {
	  markers	<- cbind(markers, res.mat)
	  if (nrow(markers) != 0)
		  colnames(markers) <- c(paste0("Mode", 1:sum(keep.mode), "_marker"),
								 paste0("Mode", 1:sum(keep.mode), "_stats"))
	  results.supp <- NULL
  }


  res.mat.ranking <- data.frame(res.mat)
  res.mat.ranking[abs(res.mat) > include.thres] <- Inf
  res.mat.ranking <- abs(res.mat.ranking)

  try(markers <- markers[do.call("order", 
								 res.mat.ranking[, 1:sum(keep.mode), drop = F]),, drop = F])

  try(res.mat <- res.mat[do.call("order", 
								 res.mat.ranking[, 1:sum(keep.mode), drop = F]),, drop = F])

 
  if (map.marker) {
	  names(marker.b_out) <- rownames(marker.b_exp)
	  marker.ratio <- marker.b_out[markers$rsID] / marker.b_exp[markers$rsID, 1]
	  marker.mode <- as.matrix(markers[, 6:(5 + sum(keep.mode)), drop = F]) %*% 
		  as.matrix(beta.mode[keep.mode])
  	  markers <- cbind(markers[, 1:4], marker.est = marker.ratio,
						 mode = marker.mode, markers[, -(1:4)])

	  tmp.est <- data.frame(x = marker.ratio, col = as.factor(marker.mode),
							y = rep(0.9 * llk.limit[1] + 0.1 * llk.limit[2], length(marker.ratio)),
							genes = markers$GENCODE_name)
	  idx <- rep(T, nrow(tmp.est))
	  for (i in 2:nrow(tmp.est)) {
		  ii <- which(tmp.est$genes[1:(i - 1)] == tmp.est$genes[i])
	  	if (length(ii) >= 1 && tmp.est$col[i] == tmp.est$col[max(ii)])
			idx[i] <- F
	  }
	  tmp.est1 <- tmp.est[idx, ]

	  p <- p + geom_point(data = tmp.est, mapping = aes(x = x, y= y, col = col, 
														fill = col), shape = "|", size = 2) + 
		  geom_text_repel(data = tmp.est1, mapping = aes(x= x, y= y, label = genes, 
														 col = col), inherit.aes = F)
		   
  }



  return(list(fun = robust.optfun.fixtau,
              modes = beta.mode[keep.mode],
              raw.modes = beta.mode, 
              p = p,
              markers = markers,
			  supp_gwas = results.supp))
}


#' findLocal maximum modes of a series of points
#'
#' @keywords internal
#' 
findLocalModes <- function(v, npts = 21) {
	v.temp <- c(rep(0, npts), v, rep(0, npts))

	idx <- c(1:npts, 1:npts + npts + 1)
	#  print(idx)
	judge <- sapply(idx, function(i) v - v.temp[i:(length(v) + i -1)] > 0)
	modes <- rowSums(judge) == ncol(judge)
	return(which(modes))
}

#' Calculate the robustified profile likelihood function
#'
#' @inheritParams findModes
#' @return The robustified profile likelihood function
#' @export
#' 
robustLossFixtau <- function(data = NULL,                                                              
							 b_exp = NULL, b_out = NULL,
							 se_exp = NULL, se_out = NULL, 
							 cor.mat = NULL,
							 loss.function = c("tukey", "huber", "l2"),
							 k.findmodes = switch(loss.function[1],
												  l2 = NA, huber = 1.345,
												  tukey = 3)) {

	if (!is.null(data)) {
		nn <- colnames(data)
		b_exp <- as.matrix(data[, grep("gamma_exp", nn), drop = F])
		se_exp <-  as.matrix(data[, grep("se_exp", nn), drop = F])
		b_out <-  data[, grep("gamma_out", nn)[1]]
		se_out <-  data[, grep("se_out", nn)[1]]
    } else {
        if (is.null(b_exp) || is.null(b_out) || is.null(se_exp) || is.null(se_out))
            stop("Require providing either data or all values of b_exp, b_out,
                 se_exp and se_out")
    }

    b_exp <- as.matrix(b_exp)
    se_exp <- as.matrix(se_exp)
    b_out <- as.vector(b_out)
    se_out <- as.vector(se_out)


    loss.function <- match.arg(loss.function, c("tukey", "huber", "l2"))
    if (is.null(cor.mat))
      cor.mat <- diag(rep(1, ncol(b_exp) + 1))
    rho <- switch(loss.function,
                  l2 = function(r, ...) rho.l2(r, ...),
                  huber = function(r, ...) rho.huber(r, k.findmodes, ...),
                  tukey = function(r, ...) rho.tukey(r, k.findmodes, ...))

   t_fun <- function(beta, tau2 = 0) {
        upper <- b_out - b_exp %*% beta
        temp <- t(t(cbind(se_exp, se_out)) * c(-beta, 1))
        lower <- sqrt(rowSums((temp %*% cor.mat) * temp) + tau2)
        return(upper/lower)
    }

   robust.optfun.fixtau <- function(beta, tau2 = 0) {
	   -sum(rho(t_fun(beta, tau2)))
    }
   return(robust.optfun.fixtau)
}




