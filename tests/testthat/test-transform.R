# The aligned model must reproduce the configural likelihood in every group:
# loadings, intercepts, mean and variance in the reported metric are a
# reparameterization of the configural model.

aligned_loglik <- function(fit, g) {
  inp <- fit$input
  mod <- fit$fits[[g]]
  dat <- mirt::extract.mirt(mod, "data")
  itype <- mirt::extract.mirt(mod, "itemtype")
  sv <- mirt::mod2values(mod)
  pr <- fit$params[fit$params$group == g, ]
  for (j in seq_along(itype)) {
    it <- colnames(dat)[j]
    a <- pr$aligned[pr$item == it & pr$param == "a"]
    cc <- pr$aligned[pr$item == it & pr$param != "a"]
    sv$value[sv$item == it & sv$name == "a1"] <- a
    if (itype[j] %in% c("2PL", "3PL")) {
      sv$value[sv$item == it & sv$name == "d"] <- cc
    } else {
      d <- if (itype[j] == "graded") cc else cumsum(cc)
      for (k in seq_along(d)) sv$value[sv$item == it & sv$name == paste0("d", k)] <- d[k]
    }
  }
  m <- fit$means[fit$means$group == g, ]
  sv$value[sv$name == "MEAN_1"] <- m$mean
  sv$value[sv$name == "COV_11"] <- m$variance
  sv$est <- FALSE
  sv0 <- mirt::mod2values(mod)
  sv0$est <- FALSE
  # both evaluated on the same wide quadrature grid
  ev <- function(pars) suppressWarnings(suppressMessages(mirt::extract.mirt(
    mirt::mirt(dat, 1, itemtype = itype, pars = pars, TOL = NaN, verbose = FALSE,
               quadpts = 201, technical = list(theta_lim = c(-12, 12))), "logLik")))
  c(aligned = ev(sv), configural = ev(sv0))
}

test_that("aligned parameters reproduce the configural likelihood (2PL, graded, gpcm)", {
  fit <- get_sim_fit()
  for (g in fit$groups) {
    ll <- aligned_loglik(fit, g)
    expect_equal(unname(ll["aligned"]), unname(ll["configural"]), tolerance = 1e-6)
  }
})

test_that("3PL slopes and intercepts are aligned and the lower asymptote is kept", {
  skip_on_cran()
  set.seed(21)
  G <- 3
  a <- c(1.4, 1.1, 1.8, 0.9, 1.3)
  d <- c(0.5, -0.2, 1, -0.8, 0.1)
  dat <- do.call(rbind, lapply(seq_len(G), function(g) {
    th <- stats::rnorm(3000, c(0, 0.4, -0.3)[g], sqrt(c(1, 1.2, 0.8)[g]))
    p <- 0.15 + 0.85 * stats::plogis(outer(th, a) + matrix(d, 3000, 5, byrow = TRUE))
    data.frame(grp = g, matrix(as.integer(stats::runif(3000 * 5) < p), 3000, 5))
  }))
  fit <- suppressWarnings(align_irt(dat, group = "grp", itemtype = c("3PL", rep("2PL", 4)),
                                    nstarts = 5, seed = 1))
  expect_equal(dim(fit$guess), c(3L, 5L))
  expect_true(all(!is.na(fit$guess[, 1])) && all(is.na(fit$guess[, -1])))
  for (g in fit$groups) {
    ll <- aligned_loglik(fit, g)
    expect_equal(unname(ll["aligned"]), unname(ll["configural"]), tolerance = 1e-6)
  }
})
