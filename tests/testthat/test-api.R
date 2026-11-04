sci <- function() {
  dat <- mirt::Science
  dat$grp <- rep(c("A", "B", "C"), length.out = nrow(dat))
  dat
}

test_that("align_fits() reproduces align_irt()", {
  s <- sim_mixed(800, alpha = c(0, 0.3, -0.2), psi = c(1, 1.2, 0.8), seed = 5,
                 noninvariance = FALSE)
  f1 <- align_irt(s$data, group = "grp", itemtype = s$itemtype, nstarts = 5, seed = 1)
  fits <- lapply(split(s$data[, -1], s$data$grp), function(d)
    mirt::mirt(d, 1, itemtype = s$itemtype, SE = TRUE, verbose = FALSE))
  f2 <- align_fits(fits, nstarts = 5, seed = 1)
  expect_equal(coef(f1), coef(f2), tolerance = 1e-8)
  expect_equal(f1$params, f2$params, tolerance = 1e-8)
})

test_that("print and summary show the main blocks", {
  fit <- get_sim_fit()
  expect_output(print(fit), "Alignment: FIXED")
  out <- capture.output(print(summary(fit)))
  expect_true(any(grepl("Approximate measurement \\(non\\)invariance", out)))
  expect_true(any(grepl("Factor mean comparison", out)))
  expect_true(any(grepl("\\(g4\\)", out)))
  expect_true(any(grepl("Item parameters", out)))
  expect_s3_class(summary(fit), "summary.irtalign")
})

test_that("graded items may have empty categories in a group", {
  s <- sim_mixed(800, alpha = c(0, 0.3, -0.2), psi = c(1, 1.2, 0.8), seed = 5,
                 noninvariance = FALSE)
  dat <- s$data
  # categories 0..3 of I4; group g2 has no 1s. Its first threshold is then the
  # cumulative logit of Y >= 2, i.e. boundary c2.
  dat$I4[dat$grp == "g2" & dat$I4 == 1] <- 2
  fit <- align_irt(dat, group = "grp", itemtype = s$itemtype, nstarts = 5, seed = 1)
  pr <- invariance_table(fit)
  i4 <- pr[pr$item == "I4", ]
  expect_false(any(i4$group == "g2" & i4$param == "c1"))
  expect_true(all(c("g1", "g3") %in% i4$group[i4$param == "c1"]))
  expect_true(all(c("c2", "c3") %in% i4$param[i4$group == "g2"]))
})

test_that("gpcm items with an empty middle category are rejected", {
  s <- sim_mixed(800, alpha = c(0, 0.3, -0.2), psi = c(1, 1.2, 0.8), seed = 5,
                 noninvariance = FALSE)
  dat <- s$data
  dat$I6[dat$grp == "g2" & dat$I6 == 1] <- 2
  expect_error(suppressWarnings(align_irt(dat, group = "grp", itemtype = s$itemtype,
                                          nstarts = 2)),
               "empty middle category")
})

test_that("items without variation in a group are left out of that group", {
  dat <- sci()
  dat$Work[dat$grp == "C"] <- 3
  expect_message(fit <- align_irt(dat, group = "grp", nstarts = 5, seed = 1),
                 "C:Work")
  pr <- invariance_table(fit)
  expect_false(any(pr$item == "Work" & pr$group == "C"))
  expect_true(any(pr$item == "Work" & pr$group == "A"))
})

test_that("unsupported models and inputs give clear errors", {
  dat <- sci()
  expect_error(align_irt(dat, group = "grp", itemtype = "Rasch"), "Unsupported itemtype")
  expect_error(align_irt(dat, group = "nope"), "column name")
  fits <- list(A = mirt::mirt(mirt::Science[1:200, ], 1, verbose = FALSE))
  expect_error(align_fits(fits), "at least two")
  two <- list(A = mirt::mirt(mirt::Science[1:200, ], 2, verbose = FALSE),
              B = mirt::mirt(mirt::Science[201:392, ], 2, verbose = FALSE))
  expect_error(align_fits(two), "unidimensional")
  a <- matrix(1, 3, 2)
  expect_error(align_pars(a, matrix(0, 2, 2), n = rep(100, 3)), "same dimensions")
  expect_error(align_pars(a, matrix(0, 3, 2), n = c(100, 100)), "sample size")
})

test_that("align_fits() matches categories by value when a range is shifted", {
  s <- sim_mixed(1500, alpha = c(0, 0.3, -0.2), psi = c(1, 1.2, 0.8), seed = 5,
                 noninvariance = FALSE)
  dat <- s$data
  # gpcm item I6 (categories 0..3): groups g1 and g3 use 0..2, group g2 uses 1..3
  dat$I6[dat$grp != "g2" & dat$I6 == 3] <- 2
  dat$I6[dat$grp == "g2" & dat$I6 == 0] <- 1
  f1 <- align_irt(dat, group = "grp", itemtype = s$itemtype, nstarts = 5, seed = 1)
  fits <- lapply(split(dat[, -1], dat$grp), function(d)
    mirt::mirt(d, 1, itemtype = s$itemtype, SE = TRUE, verbose = FALSE))
  f2 <- align_fits(fits, nstarts = 5, seed = 1)
  expect_equal(coef(f1), coef(f2), tolerance = 1e-8)
  p2 <- invariance_table(f2)
  i6 <- p2[p2$item == "I6", ]
  expect_setequal(i6$param[i6$group == "g2"], c("a", "c2", "c3"))
  expect_setequal(i6$param[i6$group == "g1"], c("a", "c1", "c2"))
})

test_that("a boundary that no group observes is dropped", {
  s <- sim_mixed(800, alpha = c(0, 0.3, -0.2), psi = c(1, 1.2, 0.8), seed = 5,
                 noninvariance = FALSE)
  dat <- s$data
  # graded item I4 (0..3): g1 never uses 1, g2 and g3 never use 0, so no group
  # has a threshold for Y >= 1
  dat$I4[dat$grp == "g1" & dat$I4 == 1] <- 2
  dat$I4[dat$grp != "g1" & dat$I4 == 0] <- 1
  fit <- align_irt(dat, group = "grp", itemtype = s$itemtype, nstarts = 5, seed = 1)
  i4 <- invariance_table(fit)
  i4 <- i4[i4$item == "I4", ]
  expect_false("c1" %in% i4$param)
  expect_true(all(c("c2", "c3") %in% i4$param))
  ti <- toy_input()
  cl <- ti$c_list
  cl$i3[, 2] <- NA
  vc <- lapply(ti$vcov, function(V) {
    keep <- rownames(V) != "i3:c2"
    V[keep, keep]
  })
  f <- align_pars(ti$a0, cl, n = ti$n, vcov = vc, nstarts = 3, seed = 1)
  expect_false(any(f$params$item == "i3" & f$params$param == "c2"))
  expect_true(any(f$params$item == "i3" & f$params$param == "c3"))
})

test_that("align_pars() matches named n, vcov and intercept rows to the groups", {
  ti <- toy_input()
  names(ti$n) <- names(ti$vcov) <- rownames(ti$a0)
  f1 <- align_pars(ti$a0, ti$c_list, n = ti$n, vcov = ti$vcov, ref = 1, nstarts = 3, seed = 1)
  o <- 5:1
  cl <- lapply(ti$c_list, function(m) m[o, , drop = FALSE])
  f2 <- align_pars(ti$a0, cl, n = ti$n[o], vcov = ti$vcov[o], ref = 1, nstarts = 3, seed = 1)
  expect_equal(coef(f1), coef(f2), tolerance = 1e-10)
  names(ti$n)[1] <- "zz"
  expect_error(align_pars(ti$a0, ti$c_list, n = ti$n, nstarts = 1), "do not match")
})

test_that("align_pars() works without covariance matrices", {
  ti <- toy_input()
  fit <- align_pars(ti$a0, ti$c_list, n = ti$n, nstarts = 3, seed = 1)
  expect_true(all(is.na(fit$params$se)))
  expect_true(all(is.na(fit$params$invariant)))
  expect_null(fit$mean_comparison)
  expect_output(print(summary(fit)), "Item parameters")
})
