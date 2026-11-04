toy_fit <- function(...) {
  ti <- toy_input()
  align_pars(ti$a0, ti$c_list, n = ti$n, vcov = ti$vcov, seed = 1, nstarts = 10, ...)
}

test_that("FIXED alignment fixes the reference mean and variance", {
  fit <- toy_fit(align = "fixed", ref = "g3")
  m <- coef(fit)
  expect_equal(m$mean[3], 0)
  expect_equal(m$variance[3], 1)
  expect_equal(m$se_mean[3], 0)
  expect_equal(m$se_variance[3], 0)
  expect_equal(fit$settings$ref, "g3")
})

test_that("the two metrics are rescalings of the same solution", {
  f1 <- toy_fit(align = "free", metric = "refgroup", ref = 2)
  f2 <- toy_fit(align = "free", metric = "product", ref = 2)
  expect_equal(f1$optimization$value, f2$optimization$value)
  expect_equal(prod(coef(f2)$variance), 1, tolerance = 1e-10)
  k <- coef(f2)$variance[2]
  expect_equal(coef(f1)$variance, coef(f2)$variance / k, tolerance = 1e-8)
  expect_equal(coef(f1)$mean, coef(f2)$mean / sqrt(k), tolerance = 1e-8)
  # aligned intercepts do not depend on the metric, slopes scale with sqrt(k)
  p1 <- f1$params
  p2 <- f2$params
  ic <- p1$param != "a"
  expect_equal(p1$aligned[ic], p2$aligned[ic], tolerance = 1e-8)
  expect_equal(p1$aligned[!ic], p2$aligned[!ic] * sqrt(k), tolerance = 1e-8)
  # intercept tests and R2 are metric free
  expect_equal(p1$invariant[ic], p2$invariant[ic])
  expect_equal(p1$se[ic], p2$se[ic], tolerance = 1e-8)
  expect_equal(f1$item_summary$r2, f2$item_summary$r2, tolerance = 1e-8)
})

test_that("FREE alignment does not depend on the order of the groups", {
  ti <- toy_input()
  f1 <- align_pars(ti$a0, ti$c_list, n = ti$n, align = "free", metric = "product",
                   eps = 0.01, se = FALSE, nstarts = 40, seed = 1)
  o <- c(4, 1, 5, 3, 2)
  f2 <- align_pars(ti$a0[o, ], lapply(ti$c_list, function(m) m[o, , drop = FALSE]),
                   n = ti$n[o], align = "free", metric = "product", eps = 0.01,
                   se = FALSE, nstarts = 40, seed = 1)
  expect_equal(f1$optimization$value, f2$optimization$value, tolerance = 1e-8)
  expect_equal(coef(f1)$mean[o], coef(f2)$mean, tolerance = 1e-5)
  expect_equal(coef(f1)$variance[o], coef(f2)$variance, tolerance = 1e-5)
})

test_that("ref = 'auto' uses the group with the free mean closest to zero", {
  ti <- toy_input()
  free <- align_pars(ti$a0, ti$c_list, n = ti$n, align = "free", metric = "product",
                     se = FALSE, nstarts = 10, seed = 1)
  auto <- align_pars(ti$a0, ti$c_list, n = ti$n, align = "fixed", ref = "auto",
                     se = FALSE, nstarts = 10, seed = 1)
  expect_equal(auto$settings$ref_index, which.min(abs(free$optimization$alpha)))
})

test_that("starting values are reproducible and keep the caller's random state", {
  ti <- toy_input()
  set.seed(99)
  before <- stats::runif(1)
  set.seed(99)
  f1 <- align_pars(ti$a0, ti$c_list, n = ti$n, ref = 1, se = FALSE, nstarts = 10, seed = 1)
  after <- stats::runif(1)
  expect_equal(before, after)
  f2 <- align_pars(ti$a0, ti$c_list, n = ti$n, ref = 1, se = FALSE, nstarts = 10, seed = 1)
  expect_identical(f1$optimization$starts, f2$optimization$starts)
  expect_equal(nrow(f1$optimization$starts), 10)
  expect_gte(f1$optimization$n_best, 1)
})

test_that("item weights change only the intercept weights", {
  f1 <- toy_fit(item_weight = "parameter", se = FALSE)
  f2 <- toy_fit(item_weight = "item", se = FALSE)
  s1 <- f1$item_summary
  s2 <- f2$item_summary
  expect_equal(s1$contribution[s1$param == "a"] > 0, s2$contribution[s2$param == "a"] > 0)
  expect_false(isTRUE(all.equal(f1$optimization$value, f2$optimization$value)))
})

test_that("FREE alignment warns with two groups", {
  ti <- toy_input()
  expect_warning(align_pars(ti$a0[c(1, 3), ],
                            lapply(ti$c_list, function(m) m[c(1, 3), , drop = FALSE]),
                            n = ti$n[c(1, 3)], align = "free", se = FALSE, nstarts = 2),
                 "not identified")
})

test_that("ref = 'auto' avoids sacrificing the reference group", {
  # small slope noninvariance everywhere makes the location of the means
  # matter; a fixed reference far from the free solution distorts it
  set.seed(92)
  G <- 30
  I <- 10
  alpha <- stats::rnorm(G, 0, 0.5)
  alpha[1] <- 1.2
  psi <- exp(stats::rnorm(G, 0, 0.3))
  a <- stats::runif(I, 0.8, 2)
  cc <- stats::rnorm(I)
  a0 <- outer(sqrt(psi), a) * matrix(exp(stats::rnorm(G * I, 0, 0.08)), G, I)
  c0 <- matrix(cc, G, I, byrow = TRUE) + outer(alpha, a) +
    matrix(stats::rnorm(G * I, 0, 0.03), G, I)
  dimnames(a0) <- dimnames(c0) <- list(paste0("g", seq_len(G)), paste0("i", seq_len(I)))
  truth <- (alpha - alpha[1]) / sqrt(psi[1])
  auto <- align_pars(a0, c0, n = rep(1000, G), metric = "product", se = FALSE,
                     nstarts = 10, seed = 1)
  est <- coef(auto)$mean
  rel <- (est - est[1])
  expect_gt(stats::cor(rel, truth), 0.99)
  expect_lt(max(abs(rel / sqrt(coef(auto)$variance[1]) - truth)), 0.15)
})
