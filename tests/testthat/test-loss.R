test_that("compiled simplicity function matches the R reference", {
  ti <- toy_input()
  inp <- irtalign:::.new_input(ti$a0, ti$c_list, ti$n, rep("custom", 4), ti$vcov)
  W <- irtalign:::.pair_weights(inp$n)
  set.seed(1)
  al <- stats::rnorm(5, sd = 0.5)
  lp <- stats::rnorm(5, sd = 0.3)
  for (pw in c(0.5, 0.25)) {
    for (eps in c(1e-4, 1e-2)) {
      cpp <- irtalign:::align_loss_cpp(inp$a0, inp$c0, inp$citem - 1L, W, rep(1, 4),
                                       rep(1, ncol(inp$c0)), al, lp, eps, pw, FALSE)$value
      ref <- loss_ref(inp$a0, inp$c0, inp$citem, inp$n, al, lp, eps, pw)
      expect_equal(cpp, ref, tolerance = 1e-12)
    }
  }
})

test_that("analytic gradient matches numerical derivatives", {
  ti <- toy_input()
  inp <- irtalign:::.new_input(ti$a0, ti$c_list, ti$n, rep("custom", 4), ti$vcov)
  W <- irtalign:::.pair_weights(inp$n)
  wc <- 1 / inp$K[inp$citem]
  set.seed(2)
  x <- c(stats::rnorm(5, sd = 0.5), stats::rnorm(5, sd = 0.3))
  f <- function(x) irtalign:::align_loss_cpp(inp$a0, inp$c0, inp$citem - 1L, W, rep(1, 4),
                                             wc, x[1:5], x[6:10], 0.01, 0.5, FALSE)$value
  g <- irtalign:::align_loss_cpp(inp$a0, inp$c0, inp$citem - 1L, W, rep(1, 4), wc,
                                 x[1:5], x[6:10], 0.01, 0.5, TRUE)$gradient
  expect_equal(g, num_grad(f, x), tolerance = 1e-6)
})

test_that("second derivatives match numerical derivatives of the gradient", {
  ti <- toy_input()
  inp <- irtalign:::.new_input(ti$a0, ti$c_list, ti$n, rep("custom", 4), ti$vcov)
  W <- irtalign:::.pair_weights(inp$n)
  ci <- inp$citem - 1L
  set.seed(4)
  al <- stats::rnorm(5, sd = 0.5)
  lp <- stats::rnorm(5, sd = 0.3)
  for (pw in c(0.5, 0.25)) {
    grad_at <- function(a0, c0, al, lp)
      irtalign:::align_loss_cpp(a0, c0, ci, W, rep(1, 4), rep(1, ncol(c0)), al, lp,
                                0.01, pw, TRUE)$gradient
    H <- irtalign:::align_hess_cpp(inp$a0, inp$c0, ci, W, rep(1, 4), rep(1, ncol(inp$c0)),
                                   al, lp, 0.01, pw, inp$ia, inp$ic, inp$M1)
    # H22 by differencing the gradient in (alpha, lpsi)
    x <- c(al, lp)
    Hn <- sapply(seq_along(x), function(k) {
      e <- replace(numeric(10), k, 1e-6)
      (grad_at(inp$a0, inp$c0, (x + e)[1:5], (x + e)[6:10]) -
         grad_at(inp$a0, inp$c0, (x - e)[1:5], (x - e)[6:10])) / 2e-6
    })
    expect_equal(H$H22, Hn, tolerance = 1e-5)
    # H21 for every configural estimate
    for (g in 1:5) for (i in 1:4) {
      if (inp$ia[g, i] < 0) next
      ap <- am <- inp$a0
      ap[g, i] <- ap[g, i] + 1e-6
      am[g, i] <- am[g, i] - 1e-6
      d <- (grad_at(ap, inp$c0, al, lp) - grad_at(am, inp$c0, al, lp)) / 2e-6
      expect_equal(H$H21[, inp$ia[g, i] + 1], d, tolerance = 1e-5)
    }
    for (g in 1:5) for (p in seq_len(ncol(inp$c0))) {
      if (inp$ic[g, p] < 0) next
      cp <- cm <- inp$c0
      cp[g, p] <- cp[g, p] + 1e-6
      cm[g, p] <- cm[g, p] - 1e-6
      d <- (grad_at(inp$a0, cp, al, lp) - grad_at(inp$a0, cm, al, lp)) / 2e-6
      expect_equal(H$H21[, inp$ic[g, p] + 1], d, tolerance = 1e-5)
    }
  }
})

test_that("the optimizer's Hessian equals the full second-derivative routine", {
  ti <- toy_input()
  inp <- irtalign:::.new_input(ti$a0, ti$c_list, ti$n, rep("custom", 4), ti$vcov)
  W <- irtalign:::.pair_weights(inp$n)
  ci <- inp$citem - 1L
  wc <- 1 / inp$K[inp$citem]
  set.seed(8)
  al <- stats::rnorm(5, sd = 0.5)
  lp <- stats::rnorm(5, sd = 0.3)
  for (pw in c(0.5, 0.25)) {
    full <- irtalign:::align_hess_cpp(inp$a0, inp$c0, ci, W, rep(1, 4), wc, al, lp, 1e-3, pw,
                                      inp$ia, inp$ic, inp$M1)$H22
    lean <- irtalign:::align_hess22_cpp(inp$a0, inp$c0, ci, W, rep(1, 4), wc, al, lp, 1e-3, pw)
    expect_equal(lean, full, tolerance = 1e-12)
  }
})

test_that("group contributions add up to the simplicity function", {
  fit <- get_sim_fit()
  expect_equal(sum(fit$group_contribution$contribution), fit$optimization$value,
               tolerance = 1e-10)
  expect_equal(sum(fit$item_summary$contribution), fit$optimization$value,
               tolerance = 1e-10)
})
