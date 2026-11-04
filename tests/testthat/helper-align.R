# Shared helpers for the tests

# Reference implementation of the simplicity function in plain R
loss_ref <- function(a0, c0, citem, n, alpha, lpsi, eps, power = 0.5,
                     wa = rep(1, ncol(a0)), wc = rep(1, ncol(c0))) {
  w <- sqrt(n) / mean(sqrt(n))
  W <- outer(w, w)
  s <- exp(-lpsi / 2)
  a1 <- a0 * s
  c1 <- c0 - (alpha * a1)[, citem, drop = FALSE]
  f <- function(x) (x^2 + eps)^(power / 2)
  tot <- 0
  for (q in seq_len(ncol(a1))) {
    d <- outer(a1[, q], a1[, q], "-")
    ut <- upper.tri(d) & !is.na(d)
    tot <- tot + wa[q] * sum(W[ut] * f(d[ut]))
  }
  for (q in seq_len(ncol(c1))) {
    d <- outer(c1[, q], c1[, q], "-")
    ut <- upper.tri(d) & !is.na(d)
    tot <- tot + wc[q] * sum(W[ut] * f(d[ut]))
  }
  tot
}

num_grad <- function(f, x, h = 1e-6) {
  vapply(seq_along(x), function(k) {
    e <- replace(numeric(length(x)), k, h)
    (f(x + e) - f(x - e)) / (2 * h)
  }, 1)
}

# Small random standard-form input with a missing item, a missing intercept
# and a 3-intercept item
toy_input <- function(G = 5, seed = 3) {
  set.seed(seed)
  items <- c("i1", "i2", "i3", "i4")
  a0 <- matrix(stats::runif(G * 4, 0.6, 1.8), G, 4,
               dimnames = list(paste0("g", seq_len(G)), items))
  c_list <- list(
    i1 = matrix(stats::rnorm(G), G, 1),
    i2 = matrix(stats::rnorm(G), G, 1),
    i3 = matrix(sort(stats::rnorm(3), decreasing = TRUE), G, 3, byrow = TRUE) +
      matrix(stats::rnorm(G * 3, sd = 0.3), G, 3),
    i4 = matrix(stats::rnorm(G * 2), G, 2)
  )
  c_list <- lapply(c_list, function(m) { rownames(m) <- rownames(a0); m })
  a0[2, "i2"] <- NA
  c_list$i3[4, 2] <- NA
  vc <- lapply(seq_len(G), function(g) {
    labs <- character()
    for (it in items) {
      if (is.na(a0[g, it])) next
      ok <- which(!is.na(c_list[[it]][g, ]))
      labs <- c(labs, paste0(it, ":a"), paste0(it, ":c", ok))
    }
    m <- length(labs)
    X <- matrix(stats::rnorm(m * (m + 3)), m + 3, m)
    V <- crossprod(X) / 400
    dimnames(V) <- list(labs, labs)
    V
  })
  list(a0 = a0, c_list = c_list, n = c(800, 1200, 500, 1000, 700)[seq_len(G)], vcov = vc)
}

# Mixed-format data for G groups from known means and variances
sim_mixed <- function(n, alpha, psi, seed, noninvariance = TRUE) {
  G <- length(alpha)
  mk <- function(type, a, c) list(type = type, a = rep(a, G),
                                  c = matrix(rep(c, each = G), nrow = G))
  pars <- list(I1 = mk("2PL", 1.2, 0.2), I2 = mk("2PL", 0.9, -0.5),
               I3 = mk("2PL", 1.5, 0.8), I4 = mk("graded", 1.1, c(1.5, 0, -1.2)),
               I5 = mk("graded", 1.4, c(2, 0.5, -0.8)),
               I6 = mk("gpcm", 0.9, c(0.5, 0.2, -0.6)), I7 = mk("gpcm", 1.2, c(1, -0.3)))
  if (noninvariance) {
    pars$I1$c[3, 1] <- pars$I1$c[3, 1] + 0.8
    pars$I5$a[min(4, G)] <- 0.6
  }
  set.seed(seed)
  resp <- irtalign:::.simulate_responses(pars, alpha, psi, rep(n, length.out = G))
  dat <- do.call(rbind, lapply(seq_len(G), function(g)
    data.frame(grp = paste0("g", g), resp[[g]])))
  list(data = dat, pars = pars,
       itemtype = c("2PL", "2PL", "2PL", "graded", "graded", "gpcm", "gpcm"))
}

fit_cache <- new.env()
get_sim_fit <- function() {
  if (is.null(fit_cache$fit)) {
    s <- sim_mixed(1500, alpha = c(0, 0.3, -0.4, 0.8, -0.2),
                   psi = c(1, 1.3, 0.8, 1.1, 0.7), seed = 11)
    fit_cache$fit <- suppressMessages(
      align_irt(s$data, group = "grp", itemtype = s$itemtype, nstarts = 10, seed = 1))
  }
  fit_cache$fit
}
