#' Robust multivariate MR estimation
#'
#' The main function of GRAPPLE to estimate causal effects of risk factors \code{beta} under a random effect model of the pleiotropic effects.
#'
#' @param data A data frame containing the information of the selected genetic instruments.
#' One can simply take the \code{data} element from the output of function \code{getInput}, or provide their own data frame.
#' The required columns include \code{SNP} for the SNP rsID, the columns \code{gamma_exp1} to \code{gamma_expk} for the
#' estimated effect sizes of the SNPs on risk factors \code{1} to \code{k}, the columns \code{se_exp1} to \code{se_expk} for their
#' standard deviations, the columns \code{gamma_out1} to \code{gamma_outm} for the
#' estimated effect sizes of the SNPs on diseases \code{1} to \code{m},the columns \code{se_out1} to \code{se_outm} for their
#' standard deviations.
#' @param p.thres The p-value threshold for SNP selection. The SNPs whose \code{selection_pvals} are less than \code{p.thres} are selected.
#' The default value is NULL, which is to include all SNPs in \code{data}. If it is not NULL, then \code{data} should have a column \code{selection_pvals}
#' that stores the selection p-value of each SNP.
#' @param cor.mat Either NULL or a \code{k + 1} by \code{k + 1} symmetric matrix. The shared correlation matrix for \code{(b_exp[j], b_out[j])} across SNP j. Used only when \code{input.list} is NULL and the default value is NULL, for the identity matrix.
#' @param tau2 The dispersion parameter. The default value is NULL, which is to be estimated by the function
#' @param loss.function Loss function used, one of "tukey", "huber" or "l2". Default is "tukey", which is robust to outlier SNPs with large pleiotropic effects
#' @param k Tuning parameters of the loss function, for loss "l2", it is NA, for loss "huber", default is 1.345 and for loss "tukey", default is 4.685
#' @param niter Number of maximum iterations allowed for optimization. Default is 20
#' @param tol Tolerance for convergence, default is the square root of the smallest positive floating number depending on the machine R is running on
#' @param opt.method the optimization used, which is one of choices the R function \code{optim} accepts. Default value is "L-BFGS-B".
#' @param diagnosis Run diagnosis analysis based on the residuals or not, default is FALSE
#' @param plot.it Whether show the QQ-plot or not if diagnosis if performed. Default is TRUE.
#'
#' @return A list with elements
#' \item{beta.hat}{Point estimates of \code{beta}}
#' \item{tau2.hat}{Point estimates of the pleiotropic effect variance \code{tau2} if the argument \code{tau2} is set to NULL}
#' \item{beta.variance}{Estimated covariance matrix of \code{beta.hat}}
#' \item{tau2.se}{Estimated standard deviation of \code{tau2.hat}}
#' \item{beta.p.vaue}{A vector of p-values where the kth element is the p-value for whether \code{beta_k = 0}}
#' \item{std.resid}{Returned if \code{diagnosis} is TRUE. A vector of standardized residuals of each SNP}
#'
#' @import nortest
#' @import ggplot2
#' @export
grappleRobustEst <- function(data,
                             p.thres = NULL,
                             cor.mat = NULL,
                             tau2 = NULL,
                             loss.function = c("tukey", "huber", "l2"),
                             k = switch(loss.function[1],
                                        l2 = NA, huber = 1.345,
                                        tukey = 4.685),
                             niter = 20,
                             tol = .Machine$double.eps^0.5,
                             opt.method = "L-BFGS-B",
                             diagnosis = T, plot.it = T) {

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


    loss.function <- match.arg(loss.function, c("tukey", "huber", "l2"))
    if (is.null(cor.mat))
      cor.mat <- diag(rep(1, ncol(b_exp) + 1))
    rho <- switch(loss.function,
                  l2 = function(r, ...) rho.l2(r, ...),
                  huber = function(r, ...) rho.huber(r, k, ...),
                  tukey = function(r, ...) rho.tukey(r, k, ...))

    delta <- integrate(function(x)  rho(x) * dnorm(x), -Inf, Inf)$value

    c1 <- integrate(function(x) rho(x, deriv = 1)^2 * dnorm(x), -Inf, Inf)$value
    c2 <- integrate(function(x) rho(x)^2 * dnorm(x), -Inf, Inf)$value - delta^2
    c4 <- integrate(function(x) rho(x, deriv = 1) * x * dnorm(x), -Inf, Inf)$value
    c3 <- c4


    ## First, calculate t for each SNP
    t_fun <- function(beta, tau2) {
      upper <- b_out - b_exp %*% beta
      temp <- t(t(cbind(se_exp, se_out)) * c(-beta, 1))
      lower <- sqrt(rowSums((temp %*% cor.mat) * temp) + tau2)

      return(upper/lower)
    }


    robust.optfun.fixtau <- function(beta, tau2) {
     -sum(rho(t_fun(beta, tau2)))
    }

    robust.E <- function(beta, tau2) {
      return(sum(rho(t_fun(beta, tau2))) - (length(b_out) - 1) * delta)
     }



    ## Initialize
    if (is.null(tau2))
      tau2.hat <- 0
    else
      tau2.hat <- tau2
    bound.beta <- apply(abs(b_out / b_exp), 2, function(v)quantile(v[is.finite(v)],
																   probs = 0.95, na.rm = T)) * 2

    bound.tau2 <- median(b_out^2) * 2
    if (ncol(b_exp) == 1) {
        beta.seq <- seq(-bound.beta, bound.beta, length.out = 5000)
      beta.hat <- beta.seq[which.max(sapply(beta.seq, robust.optfun.fixtau, tau2 = 0))]
    } else {
      beta.hat <- as.vector(lm(b_out ~ b_exp + 0)$coef)
      temp.fun <- function(bb, idx) {
        beta <- rep(0, length(beta.hat))
        beta[idx] <- bb
        return(robust.optfun.fixtau(beta, tau2=0))
      }

            for (i in 1:length(beta.hat)) {
        beta.seq <- seq(-bound.beta[i], bound.beta[i], length.out = 5000)
        beta.hat[i] <- beta.seq[which.max(sapply(beta.seq, temp.fun, idx = i))]
      }
    }

    for (iter in 1:niter) {
#		print(iter)
      beta.hat.old <- beta.hat
      tau2.hat.old <- tau2.hat
      if (is.null(tau2)) {
        temp <- robust.E(beta.hat, 0)

        if (temp < 0)
          tau2.hat <- 0
        else {
          tau2.hat <- tryCatch(uniroot(function(tau2) sum(robust.E(beta.hat, tau2)),
                                       bound.tau2 * c(0, 1),
									   extendInt = "downX",
                                       tol = .Machine$double.eps^0.25)$root,
                               error = function(e) {warning("Did not find a solution for tau2."); 0})
		}
      }


        if (opt.method == "L-BFGS-B") {
          beta.hat <- optim(beta.hat, function(beta) robust.optfun.fixtau(beta, tau2.hat),
                            method = opt.method, lower = -bound.beta, upper = bound.beta,
                            control = list(fnscale = -1))$par
        } else
          beta.hat <- optim(beta.hat, function(beta) robust.optfun.fixtau(beta, tau2.hat),
                            method = opt.method, #lower = -bound.beta, upper = bound.beta,
                            control = list(fnscale = -1))$par
        if (length(beta.hat) == 1) {
          beta.hat.temp <- beta.seq[which.max(sapply(beta.seq, robust.optfun.fixtau, tau2 = tau2.hat))]
          if (robust.optfun.fixtau(beta.hat.temp, tau2.hat) > robust.optfun.fixtau(beta.hat, tau2.hat))
            beta.hat <- beta.hat.temp
        }

        int.extend <- 0
        while (abs(beta.hat) > 0.95 * bound.beta && int.extend <= niter && opt.method == "L-BGFS-B") {
            int.extend <- int.extend + 1
            bound.beta <- bound.beta * 2
	    beta.hat <- optim(beta.hat, function(beta) robust.optfun.fixtau(beta, tau2.hat),
                            method = opt.method, lower = -bound.beta, upper = bound.beta,
                            control = list(fnscale = -1))$par
          #  beta.hat <- optim(function(beta) robust.optfun.fixtau(beta, tau2.hat),
          #                       bound.beta * c(-1, 1), maximum = TRUE)$maximum
        }
        if (int.extend == niter) {
            stop("Failed to find beta.")
        }
        if (max(abs(beta.hat.old - beta.hat) / abs(beta.hat + 1e-10)) +
            abs(tau2.hat.old - tau2.hat) / abs(tau2.hat + 1e-10) <= tol) {
            break
        }
    }

    if (max(abs(beta.hat.old - beta.hat) / abs(beta.hat + 1e-10)) +
        abs(tau2.hat.old - tau2.hat) / abs(tau2.hat + 1e-10) > tol) {
      warning("Did not converge when solving the estimating equations.
              Consider to increase niter or decrease tol.")
    }


    ## final round of coordinate-wise grid-search optimization to make the non-convex optimization
    ## more likely reach a global optimal point

    for(s in 1:5) {
      beta.hat.old <- beta.hat
      tau2.hat.old <- tau2.hat
      if (is.null(tau2)) {
        temp <- robust.E(beta.hat, 0)
        if (temp < 0)
          tau2.hat <- 0
        else
          tau2.hat <- tryCatch(uniroot(function(tau2) sum(robust.E(beta.hat, tau2)),
                                       bound.tau2 * c(0, 1), extendInt = "downX",
                                       tol = bound.tau2 * .Machine$double.eps^0.25)$root,
                               error = function(e) {warning("Did not find a solution for tau2."); 0})
      }
      temp.fun <- function(bb, idx, beta.hat.temp = beta.hat) {
        beta.hat.temp[idx] <- bb
        return(robust.optfun.fixtau(beta.hat.temp, tau2=tau2.hat))
      }
      for (h in 1:5) {
        beta.temp <- beta.hat
        for (i in 1:length(beta.hat)) {
          bb <- beta.hat[i] * 50
          beta.seq <- seq(-abs(bb), abs(bb), length.out = 5000)
          tt <- beta.hat
          tt[i] <- beta.seq[which.max(sapply(beta.seq, temp.fun, idx = i, beta.hat.temp = beta.hat))]
          if (robust.optfun.fixtau(tt, tau2.hat) > robust.optfun.fixtau(beta.hat, tau2.hat))
            beta.hat <- tt
        }
        if (max(abs(beta.hat.old - beta.hat) / abs(beta.hat + 1e-10)) <= tol)
          break
      }
      if (max(abs(beta.hat.old - beta.hat) / abs(beta.hat + 1e-10)) +
          abs(tau2.hat.old - tau2.hat) / abs(tau2.hat + 1e-10) <= tol) {
        break
      }
    }

    if (s == 5)
      warning("Did not converge in the final step. May not find the global optimal point.")


    ### calculating asymptotic variance
    ## calculate the derivatives of t_fun at beta.hat, tau2.hat
    r <- ncol(b_exp)
    p <- nrow(b_exp)
    res <- as.vector(b_out - b_exp %*% beta.hat) ## size p
    temp <- t(t(cbind(se_exp, se_out)) * c(-beta.hat, 1))  ## size p * (r + 1)
    var.res <- rowSums((temp %*% cor.mat) * temp) + tau2.hat ## size p
    t.values <- res/sqrt(var.res) ## size p
    ## calculate each row as Sigma_X_j beta - Sigma_{X_jY_j}
    coefs <- se_exp * t(cor.mat[1:r, 1:r] %*% (beta.hat * t(se_exp))) -
             t(cor.mat[1:r, r + 1] * t(se_exp * se_out)) ## size p * r
    ## matrix of u_j
    u.values <- - (var.res * b_exp + res * coefs) / var.res^{3/2}
    S <- t(u.values) %*% u.values

    ## calculate B
    B <- matrix(0, r + 1, r + 1)
    B[1:r, 1:r] <- c1 * S
    B[r + 1, r + 1] <- p * c2

    ## calculate A
    A <- matrix(0, r + 1, r + 1)
    A[r + 1, r+1] <- - c3/2 * sum(1/var.res)
    A[1:r, r+1] <- c4 * colSums(coefs / var.res^2)


    rho.dev.var <- rho(t.values, deriv = 1) / var.res^{3/2} ## size p
    temp <- array(0, c(p, r, r))
    rho.se <- rho.dev.var * se_exp
    for (j in 1:p)
      temp[j, , ] <- res[j] * t(rho.se[j, ] * cor.mat[1:r, 1:r]) * se_exp[j, ]
    temp <- colSums(temp, dim = 1)
    temp1 <- t(rho.dev.var * b_exp) %*% coefs
    A[1:r, 1:r] <- c4 * S  - (2 * temp1 - t(temp1) + temp)

    A.inv <- solve(A)
    asymp.var <- A.inv %*% B %*% t(A.inv)

   if (diagnosis) {
        std.resid <- res / sqrt(var.res)
   		qq.result <- qqDiagnosis(std.resid, plot.it = plot.it)
   }

   if (is.null(tau2))
     tau2.se <- sqrt(asymp.var[nrow(asymp.var), nrow(asymp.var)])
   else
     tau2.se <- NA

    out <- list(beta.hat = beta.hat,
                tau2.hat = tau2.hat,
                beta.var = asymp.var[1:ncol(b_exp), 1:ncol(b_exp)],
                tau2.se = tau2.se, # / sqrt(efficiency),
                beta.p.value = pmin(1, 2 * pnorm(abs(beta.hat) / sqrt(diag(asymp.var)[1:r]),
                                                 lower.tail = F)))# / sqrt(efficiency),
   if (diagnosis) {
        out$std.resid <- std.resid
   		out$p <- qq.result$p
		out$outliers <- qq.result$outliers
   }


    out

}

#' Huber loss function and its derivatives
#'
#' @keywords internal
#'
rho.huber <- function(r, k = 1.345, deriv = 0) {
    if (deriv == 0) {
        return(ifelse(abs(r) <= k,
                      r^2/2,
                      k * (abs(r) - k/2)))
    } else if (deriv == 1) {
        return(ifelse(abs(r) <= k,
                      r,
                      k * sign(r)))
    } else if (deriv == 2) {
        return(ifelse(abs(r) <= k,
                      1,
                      0))
    } else {
        stop("deriv must be 0, 1, or 2.")
    }
}

#' Tukey's beweight loss function and its derivatives
#'
#' @param r Function value
#' @param k Tuning parameter, default value is 4.685
#' @param deriv The derivative function to calculate. 0 is the Tukey's loss function, 1 is for the first derivative and 2 for the second derivative function
#'
#' @return The value of the corresponding function at \code{r}
#' @export
#'
rho.tukey <- function(r, k = 4.685, deriv = 0) {
    if (deriv == 0) {
        pmin(1 - (1 - (r/k)^2)^3, 1)
    } else if (deriv == 1) {
        r * (1 - (r / k)^2)^2 * (abs(r) <= k)
    } else if (deriv == 2) {
        t <- (r/k)^2
        ifelse(t < 1, (1 - t) * (1 - 5 * t), 0)
    } else {
		stop("Higher order derivatives not implemented")
	}
}

#' L2 loss function and its derivatives
#'
#' @keywords internal
#'
rho.l2 <- function(r, deriv = 0) {
    if (deriv == 0) {
        return(r^2/2)
    } else if (deriv == 1) {
        return(r)
    } else if (deriv == 2) {
        return(1)
    } else {
        stop("deriv must be 0, 1, or 2.")
    }
}
