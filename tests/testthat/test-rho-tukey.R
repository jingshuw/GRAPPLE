test_that("rho.tukey derivatives are the derivatives of rho.tukey", {
  k <- 4.685
  r <- c(-6, -4, -2.5, -1, -0.3, 0, 0.4, 1.5, 3, 4.5, 7)
  h <- 1e-5
  d1 <- (rho.tukey(r + h, k) - rho.tukey(r - h, k)) / (2 * h)
  d2 <- (rho.tukey(r + h, k, deriv = 1) - rho.tukey(r - h, k, deriv = 1)) / (2 * h)
  expect_equal(rho.tukey(r, k, deriv = 1), d1, tolerance = 1e-6)
  expect_equal(rho.tukey(r, k, deriv = 2), d2, tolerance = 1e-6)
})

test_that("Tukey tau2.se is on the same scale as the L2 and Huber ones", {
  ## No outliers, so the three losses should give similar tau2 and tau2.se.
  set.seed(1)
  m <- 40
  gam <- runif(m, 0.05, 0.2)
  sx <- runif(m, 0.01, 0.03)
  sy <- runif(m, 0.01, 0.03)
  bx <- rnorm(m, gam, sx)
  by <- rnorm(m, 0.5 * gam + rnorm(m, 0, 0.02), sy)
  dat <- data.frame(SNP = paste0("v", 1:m), gamma_exp1 = bx, se_exp1 = sx,
                    gamma_out1 = by, se_out1 = sy)
  fit <- lapply(c(tukey = "tukey", huber = "huber", l2 = "l2"), function(loss)
    suppressWarnings(grappleRobustEst(dat, loss.function = loss,
                                      diagnosis = FALSE, plot.it = FALSE)))
  se <- sapply(fit, function(f) f$tau2.se)
  expect_true(all(sapply(fit, function(f) f$tau2.hat) > 0))
  expect_gt(se[["tukey"]] / se[["l2"]], 0.7)
  expect_lt(se[["tukey"]] / se[["l2"]], 1.4)
  expect_gt(se[["tukey"]] / se[["huber"]], 0.7)
  expect_lt(se[["tukey"]] / se[["huber"]], 1.4)
})
