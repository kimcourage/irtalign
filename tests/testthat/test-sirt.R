# sirt::invariance.alignment() with meth = 3 uses the same simplicity function
# and the product constraint on the variances. Means and variances relative to
# the first group must agree.

test_that("alignment agrees with sirt::invariance.alignment()", {
  skip_if_not_installed("sirt")
  skip_on_cran()
  e <- new.env()
  utils::data("data.ess2005", package = "sirt", envir = e)
  lambda <- e$data.ess2005$lambda
  nu <- e$data.ess2005$nu
  G <- nrow(lambda)
  rownames(lambda) <- rownames(nu) <- paste0("C", seq_len(G))
  ours <- align_pars(lambda, nu, n = rep(1000, G), ref = 1, eps = 1e-3, nstarts = 30, seed = 1,
                     se = FALSE)
  s <- suppressWarnings(sirt::invariance.alignment(
    lambda, nu, wgt = matrix(1, G, ncol(nu)), align.pow = c(0.5, 0.5),
    align.scale = c(1, 1), eps = 1e-3, meth = 3, fixed = TRUE))
  rel <- function(al, ps) cbind((al - al[1]) / sqrt(ps[1]), ps / ps[1])
  r_ours <- rel(ours$optimization$alpha, exp(ours$optimization$lpsi))
  r_sirt <- rel(s$pars$alpha0, s$pars$psi0^2)
  expect_equal(unname(r_ours), unname(r_sirt), tolerance = 1e-3)
})
