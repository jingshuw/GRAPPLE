#' QQ-plot diagnosis and outlier detection from standardized residuals
#'
#' @inheritParams grappleRobustEst
#' @param std.residuals A vector of standardized residuals, can be the output element \code{std.resid} from function \code{grappleRobustEst}
#' @param outlier.quantile The quantile threshold for outliers. 
#'
#' @return A list with two elements
#' \item{p}{The QQ plot}
#' \item{outliers}{The data frame for the detected outliers}
#'
#' @import ggplot2
#' @import nortest
#' @import dplyr
#' @export
qqDiagnosis <- function(std.residuals, outlier.quantile = 0.1/length(std.residuals),
					   plot.it = T) {
	thres <- qnorm(1 - outlier.quantile/2)
	qtls <- qqnorm(std.residuals, plot.it = F)
	dd <- data.frame(theoretical = qtls$x, 
					 sample = qtls$y)
	rownames(dd) <- names(qtls$y)
	p <- ggplot(dd, aes(x = theoretical, y = sample)) + geom_point(alpha = 0.3) + 
		geom_abline(intercept = 0, slope = 1, color = "blue") + 
		theme_classic() + 
		ggtitle("QQ plot of standardized residuals")

	print("Test of normality for the standardized residuals:")
	print(paste("Anderson-Darling test: p-value =", 
				signif(ad.test(std.residuals)$p, 3)))
	print(paste("Shapiro-Wilk test: p-value =",
                signif(shapiro.test(std.residuals)$p, 3)))

	outliers <- which(abs(std.residuals) > qnorm(1 - outlier.quantile/2))
	outliers <- dd %>% filter(abs(sample) > thres)

	print(paste(nrow(outliers), "outliers detected!"))

	p <- p + geom_hline(yintercept = c(thres, -thres), color = "red", linetype = "dashed") + 
		geom_point(data = outliers, aes(x = theoretical, y = sample), color = "red")

	if (plot.it)
		print(p)

	return(list(p = p, outliers = outliers))


}
