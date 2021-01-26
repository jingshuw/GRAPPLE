#' Compute the conditional Q-statistic for assessing instrument strength
#'
#' @param dat.list Object returned from \code{getInput}
#' @inheritParams grappleRobustEst
#'
#' @return The conditional Q-statistics, degrees of freedom, and the corresponding p-values
#'
#' @references Eleanor Sanderson, George Davey Smith, Frank Windmeijer, Jack Bowden, An examination of multivariable Mendelian randomization in the single-sample and two-sample summary data settings, International Journal of Epidemiology, Volume 48, Issue 3, June 2019, Pages 713–727, https://doi.org/10.1093/ije/dyy262.
#'
computeQ <- function(dat.list, p.thres = NULL) {

    data <- dat.list$data
    if (!is.null(p.thres))
        data <- data[data$selection_pvals < p.thres, ]

    nsnps <- nrow(data)
    cor.mat <- dat.list$cor.mat
    nn <- colnames(data)
    b_exp <- as.matrix(data[, grep("gamma_exp", nn), drop = F])
    se_exp <-  as.matrix(data[, grep("se_exp", nn), drop = F])
    k <- ncol(b_exp)

    Q.stats <- sapply(1:k, function(i) {
        print(paste("Computing the conditional Q-stat for exposure ", i))
        temp.dat <- cbind(b_exp[, -i], b_exp[, i], se_exp[, -i], se_exp[, i])
        colnames(temp.dat) <- c(paste0("gamma_exp", 1:(k-1)), "gamma_out1",
                                paste0("se_exp", 1:(k-1)), "se_out1")
        ss <- 1:k
        temp.cor <- cor.mat[c(ss[-i], i), c(ss[-i], i)]

        delta <- grappleRobustEst(temp.dat, cor.mat = temp.cor,
                                  diagnosis = F)$beta.hat

        upper <- (b_exp[, i] - b_exp[, -i, drop = F] %*% delta)^2
        temp <- t(t(cbind(se_exp[, -i], se_exp[, i])) * c(-delta, 1))
    	lower <- rowSums((temp %*% temp.cor) * temp)
        Q.stat <- sum(upper/lower)
        return(Q.stat)
    })
    names(Q.stats) <- paste0("exposure", 1:k)

    df <- nrow(data) - 1
    Q.stats.pval <- pchisq(Q.stats, df, lower.tail = F)
    Q.stats.log.pval <- pchisq(Q.stats, df, lower.tail = F, log.p = T)
    return(list(Q.stats = Q.stats,
                df = df,
                Q.stats.pval = Q.stats.pval,
                Q.stats.log.pval = Q.stats.log.pval))
}
