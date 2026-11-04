# Large-sample recovery of known means and variances with mixed 2PL, graded
# and generalized partial credit items, one noninvariant intercept and one
# noninvariant slope.

test_that("known means and variances are recovered", {
  skip_on_cran()
  alpha <- c(0, 0.5, -0.4, 0.3)
  psi <- c(1, 1.5, 0.7, 1.2)
  s <- sim_mixed(5000, alpha, psi, seed = 2026)
  fit <- align_irt(s$data, group = "grp", itemtype = s$itemtype, nstarts = 10, seed = 1)
  m <- coef(fit)
  expect_lt(max(abs(m$mean - alpha)), 0.1)
  expect_lt(max(abs(m$variance - psi)), 0.2)
  z <- (m$mean - alpha)[-1] / m$se_mean[-1]
  expect_true(all(abs(z) < 3))
  tab <- invariance_table(fit)
  bad <- tab[!tab$invariant, ]
  expect_true(any(bad$item == "I1" & bad$param == "c1" & bad$group == "g3"))
  expect_true(any(bad$item == "I5" & bad$param == "a" & bad$group == "g4"))
  # aligned parameters of invariant items recover the generating values
  gp <- tab[tab$item == "I6" & tab$param == "a", ]
  expect_lt(max(abs(gp$aligned - 0.9)), 0.15)
  gc <- tab[tab$item == "I4", ]
  expect_lt(max(abs(gc$aligned[gc$param == "c2"] - 0)), 0.15)
})
