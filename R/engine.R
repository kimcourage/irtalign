# Alignment optimization --------------------------------------------------
#
# Full alignment parameters: alpha_1..alpha_G (means) and l_1..l_G (log
# variances) in the optimization metric, where sum(l) = 0. The reduced vector
# theta holds the free means (all groups for FREE, all but the reference group
# for FIXED) and l_2..l_G, with l_1 = -sum(l_2..l_G).

.theta_map <- function(G, align, ref) {
  afree <- if (align == "free") seq_len(G) else setdiff(seq_len(G), ref)
  k <- length(afree)
  q <- k + G - 1L
  J <- matrix(0, 2 * G, q)
  J[cbind(afree, seq_len(k))] <- 1
  for (g in 2:G) {
    J[G + g, k + g - 1L] <- 1
    J[G + 1L, k + g - 1L] <- -1
  }
  list(J = J, afree = afree, q = q, G = G)
}

.theta_from_full <- function(alpha, lpsi, map) {
  lpsi <- lpsi - mean(lpsi)
  c(alpha[map$afree], lpsi[-1])
}

.full_from_theta <- function(theta, map) {
  full <- drop(map$J %*% theta)
  G <- map$G
  list(alpha = full[seq_len(G)], lpsi = full[G + seq_len(G)])
}

# Least-squares linking start: log slopes give the variances, then means are
# found by alternating least squares on the intercepts.
.ls_start <- function(inp) {
  la <- log(pmax(abs(inp$a0), 1e-3))
  la <- sweep(la, 2, colMeans(la, na.rm = TRUE))
  lpsi <- 2 * rowMeans(la, na.rm = TRUE)
  lpsi[!is.finite(lpsi)] <- 0
  lpsi <- lpsi - mean(lpsi)
  a1 <- inp$a0 * exp(-lpsi / 2)
  a1c <- a1[, inp$citem, drop = FALSE]
  alpha <- rep(0, nrow(inp$a0))
  for (iter in seq_len(50)) {
    nu <- colMeans(inp$c0 - alpha * a1c, na.rm = TRUE)
    res <- sweep(inp$c0, 2, nu)
    num <- rowSums(a1c * res, na.rm = TRUE)
    den <- rowSums(a1c^2 * !is.na(res), na.rm = TRUE)
    alpha_new <- ifelse(den > 0, num / den, 0)
    alpha_new <- alpha_new - mean(alpha_new)
    if (max(abs(alpha_new - alpha)) < 1e-8) {
      alpha <- alpha_new
      break
    }
    alpha <- alpha_new
  }
  list(alpha = alpha, lpsi = lpsi)
}

.loss_fun <- function(inp, W, wa, wc, eps, power) {
  a0 <- inp$a0
  c0 <- inp$c0
  ci <- as.integer(inp$citem - 1L)
  function(alpha, lpsi, gradient = FALSE) {
    align_loss_cpp(a0, c0, ci, W, wa, wc, alpha, lpsi, eps, power, gradient)
  }
}

# Objective and gradient in theta, sharing one evaluation per point.
.objective <- function(loss, map) {
  cache <- new.env(parent = emptyenv())
  eval_at <- function(th) {
    if (!identical(th, cache$th)) {
      f <- .full_from_theta(th, map)
      r <- loss(f$alpha, f$lpsi, TRUE)
      cache$th <- th
      cache$value <- r$value
      cache$grad <- drop(crossprod(map$J, r$gradient))
    }
  }
  list(value = function(th) { eval_at(th); cache$value },
       gradient = function(th) { eval_at(th); cache$grad })
}

.optimize_alignment <- function(inp, align, ref, eps, power, nstarts, wa, wc,
                                seed = NULL, control = list(), extra_start = NULL) {
  G <- length(inp$groups)
  W <- .pair_weights(inp$n)
  map <- .theta_map(G, align, ref)
  loss <- .loss_fun(inp, W, wa, wc, eps, power)
  ob <- .objective(loss, map)
  obj <- ob$value
  grd <- ob$gradient
  ci <- as.integer(inp$citem - 1L)
  hes <- function(th) {
    f <- .full_from_theta(th, map)
    H <- align_hess22_cpp(inp$a0, inp$c0, ci, W, wa, wc, f$alpha, f$lpsi, eps, power)
    crossprod(map$J, H %*% map$J)
  }

  ls <- .ls_start(inp)
  starts <- .with_seed(seed, {
    lst <- list(
      list(alpha = ls$alpha, lpsi = ls$lpsi),
      list(alpha = rep(0, G), lpsi = rep(0, G))
    )
    extra <- max(nstarts - 2L, 0L)
    for (s in seq_len(extra)) {
      lst[[length(lst) + 1L]] <- list(alpha = stats::rnorm(G), lpsi = stats::rnorm(G, sd = 0.5))
    }
    lst[seq_len(max(nstarts, 1L))]
  })
  if (!is.null(extra_start)) starts <- c(starts, list(extra_start))

  ctrl <- utils::modifyList(list(iter.max = 2000, eval.max = 4000, rel.tol = 1e-12,
                                 x.tol = 1e-10), control)
  runs <- lapply(seq_along(starts), function(s) {
    st <- starts[[s]]
    alpha <- st$alpha
    if (align == "fixed") alpha <- alpha - alpha[ref]
    th0 <- .theta_from_full(alpha, st$lpsi, map)
    o <- tryCatch(stats::nlminb(th0, obj, grd, hes, control = ctrl),
                  error = function(e) NULL)
    if (is.null(o) || !is.finite(o$objective))
      return(list(theta = th0, value = NA_real_, convergence = 1L, iterations = NA_integer_))
    list(theta = o$par, value = o$objective, convergence = o$convergence,
         iterations = o$iterations)
  })
  vals <- vapply(runs, `[[`, 1, "value")
  if (all(is.na(vals))) stop("Alignment optimization failed from every start.", call. = FALSE)
  best <- which.min(vals)
  fmin <- vals[best]
  tol <- 1e-6 * max(1, abs(fmin))
  n_best <- sum(abs(vals - fmin) <= tol, na.rm = TRUE)
  th <- runs[[best]]$theta
  full <- .full_from_theta(th, map)
  list(theta = th, alpha = full$alpha, lpsi = full$lpsi, value = fmin,
       n_best = n_best, map = map, W = W,
       starts = data.frame(start = seq_along(runs), value = vals,
                           convergence = vapply(runs, `[[`, 1L, "convergence"),
                           iterations = vapply(runs, function(r) as.integer(r$iterations), 1L)),
       converged = runs[[best]]$convergence == 0)
}
