# Standard errors rest on d theta / d m1 = -H22^{-1} H21. Check it against
# re-optimizing the simplicity function after perturbing configural estimates.

# Solve the first-order conditions again by Newton steps from the previous
# solution (the Hessian itself is checked against numerical derivatives in
# test-loss.R).
reoptimize <- function(inp, opt, a0 = inp$a0, c0 = inp$c0, eps, wc) {
  inp2 <- inp
  inp2$a0 <- a0
  inp2$c0 <- c0
  wa <- rep(1, ncol(a0))
  loss <- irtalign:::.loss_fun(inp2, opt$W, wa, wc, eps, 0.5)
  ob <- irtalign:::.objective(loss, opt$map)
  th <- opt$theta
  for (it in 1:50) {
    g <- ob$gradient(th)
    if (max(abs(g)) < 1e-13) break
    f <- irtalign:::.full_from_theta(th, opt$map)
    H <- irtalign:::align_hess_cpp(inp2$a0, inp2$c0, as.integer(inp2$citem - 1L), opt$W,
                                   wa, wc, f$alpha, f$lpsi, eps, 0.5, inp2$ia, inp2$ic,
                                   inp2$M1)$H22
    Hr <- crossprod(opt$map$J, H %*% opt$map$J)
    th <- th - solve(Hr, g)
  }
  th
}

test_that("implicit derivatives match re-optimization", {
  ti <- toy_input()
  inp <- irtalign:::.new_input(ti$a0, ti$c_list, ti$n, rep("custom", 4), ti$vcov)
  wc <- rep(1, ncol(inp$c0))
  eps <- 0.01
  for (align in c("fixed", "free")) {
    opt <- irtalign:::.optimize_alignment(inp, align, 2L, eps, 0.5, 10, rep(1, 4), wc, seed = 1)
    vc <- irtalign:::.align_vcov(inp, opt, rep(1, 4), wc, eps, 0.5)
    expect_true(vc$ok)
    h <- 1e-5
    for (gi in list(c(1, 1), c(3, 4), c(5, 3))) {
      g <- gi[1]
      i <- gi[2]
      ap <- am <- inp$a0
      ap[g, i] <- ap[g, i] + h
      am[g, i] <- am[g, i] - h
      d <- (reoptimize(inp, opt, a0 = ap, eps = eps, wc = wc) -
              reoptimize(inp, opt, a0 = am, eps = eps, wc = wc)) / (2 * h)
      expect_equal(vc$A[, inp$ia[g, i] + 1], d, tolerance = 1e-4)
    }
    for (gp in list(c(2, 1), c(4, 5))) {
      g <- gp[1]
      p <- gp[2]
      cp <- cm <- inp$c0
      cp[g, p] <- cp[g, p] + h
      cm[g, p] <- cm[g, p] - h
      d <- (reoptimize(inp, opt, c0 = cp, eps = eps, wc = wc) -
              reoptimize(inp, opt, c0 = cm, eps = eps, wc = wc)) / (2 * h)
      expect_equal(vc$A[, inp$ic[g, p] + 1], d, tolerance = 1e-4)
    }
  }
})

test_that("delta-method covariances match brute-force propagation", {
  ti <- toy_input()
  eps <- 0.01
  fit <- align_pars(ti$a0, ti$c_list, n = ti$n, vcov = ti$vcov, eps = eps, ref = 2,
                    nstarts = 10, seed = 1)
  inp <- fit$input
  opt <- irtalign:::.optimize_alignment(inp, "fixed", 2L, eps, 0.5, 10, rep(1, 4),
                                        rep(1, ncol(inp$c0)), seed = 1)
  wc <- rep(1, ncol(inp$c0))
  # reported quantities as a function of the stacked configural vector m1
  m1 <- numeric(inp$M1)
  for (g in seq_along(inp$groups)) {
    for (i in seq_along(inp$items)) if (inp$ia[g, i] >= 0) m1[inp$ia[g, i] + 1] <- inp$a0[g, i]
    for (p in seq_len(ncol(inp$c0))) if (inp$ic[g, p] >= 0) m1[inp$ic[g, p] + 1] <- inp$c0[g, p]
  }
  quantities <- function(m) {
    a0 <- inp$a0
    c0 <- inp$c0
    for (g in seq_along(inp$groups)) {
      for (i in seq_along(inp$items)) if (inp$ia[g, i] >= 0) a0[g, i] <- m[inp$ia[g, i] + 1]
      for (p in seq_len(ncol(c0))) if (inp$ic[g, p] >= 0) c0[g, p] <- m[inp$ic[g, p] + 1]
    }
    th <- reoptimize(inp, opt, a0 = a0, c0 = c0, eps = eps, wc = wc)
    f <- irtalign:::.full_from_theta(th, opt$map)
    o2 <- list(alpha = f$alpha, lpsi = f$lpsi)
    inp2 <- inp
    inp2$a0 <- a0
    inp2$c0 <- c0
    rp <- irtalign:::.reported(inp2, o2, "refgroup", 2L)
    c(rp$alpha, rp$psi, rp$a[, "i3"], rp$c[, "i4:c2"])
  }
  Jn <- sapply(seq_along(m1), function(k) {
    e <- replace(numeric(length(m1)), k, 1e-5)
    (quantities(m1 + e) - quantities(m1 - e)) / 2e-5
  })
  V1 <- matrix(0, inp$M1, inp$M1)
  for (g in seq_along(inp$groups)) V1[inp$m1_index[[g]], inp$m1_index[[g]]] <- inp$vcov[[g]]
  se_brute <- unname(sqrt(diag(Jn %*% V1 %*% t(Jn))))
  pr <- fit$params
  se_pkg <- unname(c(fit$means$se_mean, fit$means$se_variance,
              pr$se[pr$item == "i3" & pr$param == "a"],
              pr$se[pr$item == "i4" & pr$param == "c2"]))
  expect_equal(se_pkg, se_brute, tolerance = 1e-4)
})
