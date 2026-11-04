test_that("connected components are found", {
  adj <- matrix(FALSE, 6, 6)
  adj[cbind(c(1, 2, 4), c(2, 3, 5))] <- TRUE
  adj <- adj | t(adj)
  diag(adj) <- TRUE
  comp <- irtalign:::.components(adj)
  expect_equal(comp, c(1L, 1L, 1L, 2L, 2L, 3L))
})

test_that("the invariant set excludes clearly different groups", {
  x <- c(0, 0.01, -0.02, 1, 0.015, -0.9)
  S <- diag(0.02^2, 6)
  n <- c(100, 200, 300, 400, 500, 600)
  d <- irtalign:::.detect_one(x, S, n, 0.01, 0.001)
  expect_equal(d$invariant, c(TRUE, TRUE, TRUE, FALSE, TRUE, FALSE))
  w <- n[d$invariant] / sum(n[d$invariant])
  expect_equal(d$mean, sum(w * x[d$invariant]))
})

test_that("ties between starting sets go to the larger total sample size", {
  x <- c(0, 0.01, 1, 1.01)
  S <- diag(0.01^2, 4)
  d <- irtalign:::.detect_one(x, S, c(100, 100, 500, 500), 0.01, 0.001)
  expect_equal(d$invariant, c(FALSE, FALSE, TRUE, TRUE))
})

test_that("planted noninvariance is detected in simulated data", {
  fit <- get_sim_fit()
  tab <- invariance_table(fit)
  flagged <- tab[!tab$invariant, c("item", "param", "group")]
  expect_true(any(flagged$item == "I1" & flagged$param == "c1" & flagged$group == "g3"))
  expect_true(any(flagged$item == "I5" & flagged$param == "a" & flagged$group == "g4"))
  expect_lte(nrow(flagged), 4)
})

test_that("R-squared is one for exactly invariant parameters", {
  G <- 4
  alpha <- c(0, 0.5, -0.3, 0.2)
  psi <- c(1, 1.4, 0.7, 1.1)
  a <- c(1.2, 0.8, 1.5)
  cc <- c(0.3, -0.4, 1)
  a0 <- outer(sqrt(psi), a)
  c0 <- matrix(cc, G, 3, byrow = TRUE) + outer(alpha, a)
  dimnames(a0) <- dimnames(c0) <- list(paste0("g", 1:G), paste0("i", 1:3))
  fit <- align_pars(a0, c0, n = rep(500, G), ref = 1, se = FALSE, nstarts = 5, seed = 1)
  expect_equal(coef(fit)$mean, alpha, tolerance = 1e-6)
  expect_equal(coef(fit)$variance, psi, tolerance = 1e-6)
  expect_equal(fit$item_summary$r2, rep(1, 6), tolerance = 1e-8)
})

test_that("factor mean comparison lists groups with smaller means", {
  fit <- get_sim_fit()
  mc <- mean_comparison(fit)
  expect_equal(mc$group[1], fit$means$group[which.max(fit$means$mean)])
  expect_true(all(diff(mc$mean) <= 0))
  for (r in seq_len(nrow(mc))) {
    sm <- mc$smaller[[r]]
    expect_true(all(fit$means$mean[match(sm, fit$means$group)] < mc$mean[r]))
  }
})
