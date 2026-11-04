# Reported metric and standard errors ----------------------------------------
#
# Reported metric: metric = "refgroup" sets the reference group variance to 1
# (kappa = l_ref); metric = "product" keeps the product of variances at 1
# (kappa = 0). With psi' = exp(l - kappa) and alpha' = alpha * exp(-kappa / 2):
#   slope      a' = a0 / sqrt(psi')
#   intercept  c' = c0 - alpha' * a0 / sqrt(psi')

.kappa <- function(lpsi, metric, ref) if (metric == "refgroup") lpsi[ref] else 0

.reported <- function(inp, opt, metric, ref) {
  kap <- .kappa(opt$lpsi, metric, ref)
  psi <- exp(opt$lpsi - kap)
  alpha <- opt$alpha * exp(-kap / 2)
  a <- inp$a0 / sqrt(psi)
  slope_g <- inp$a0 * exp(-opt$lpsi / 2)
  cc <- inp$c0 - (opt$alpha * slope_g)[, inp$citem, drop = FALSE]
  list(alpha = alpha, psi = psi, a = a, c = cc)
}

# Covariance pieces from the implicit function theorem:
#   d theta / d m1 = A = -H22^{-1} H21,  Var(theta) = A V1 A',  B = V1 A'.
.align_vcov <- function(inp, opt, wa, wc, eps, power) {
  map <- opt$map
  H <- align_hess_cpp(inp$a0, inp$c0, as.integer(inp$citem - 1L), opt$W, wa, wc,
                      opt$alpha, opt$lpsi, eps, power, inp$ia, inp$ic, inp$M1)
  H22 <- crossprod(map$J, H$H22 %*% map$J)
  H21 <- crossprod(map$J, H$H21)
  ev <- eigen(H22, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= 1e-12 * max(abs(ev)))
    return(list(ok = FALSE, eigen = ev))
  A <- -solve(H22, H21)
  B <- matrix(0, inp$M1, map$q)
  for (g in seq_along(inp$vcov)) {
    idx <- inp$m1_index[[g]]
    B[idx, ] <- inp$vcov[[g]] %*% t(A[, idx, drop = FALSE])
  }
  C22 <- A %*% B
  C22 <- (C22 + t(C22)) / 2
  list(ok = TRUE, A = A, B = B, C22 = C22, H22 = H22, eigen = ev,
       cond = max(ev) / min(ev))
}

# Covariance of quantities x. J1 is the Jacobian with respect to the
# configural estimates in sparse form (one row per group: m1 positions `idx`,
# values `val`, group `grp`); K is the Jacobian with respect to theta.
.qcov <- function(J1, K, vc, inp) {
  S <- K %*% vc$C22 %*% t(K)
  if (!is.null(J1)) {
    n <- nrow(J1$idx)
    Tm <- matrix(0, n, ncol(vc$B))
    D <- numeric(n)
    for (k in seq_len(ncol(J1$idx))) {
      ok <- !is.na(J1$idx[, k])
      Tm[ok, ] <- Tm[ok, ] + J1$val[ok, k] * vc$B[J1$idx[ok, k], , drop = FALSE]
    }
    for (r in seq_len(n)) {
      ok <- !is.na(J1$idx[r, ])
      g <- J1$grp[r]
      loc <- match(J1$idx[r, ok], inp$m1_index[[g]])
      v <- J1$val[r, ok]
      D[r] <- sum(v * (inp$vcov[[g]][loc, loc, drop = FALSE] %*% v))
    }
    S <- S + diag(D, n) + Tm %*% t(K) + K %*% t(Tm)
  }
  (S + t(S)) / 2
}

# Jacobians in the reported metric. Each returns list(J1, K, rows) for the
# groups in `rows`.
.jac_means <- function(inp, opt, metric, ref) {
  G <- length(inp$groups)
  kap <- .kappa(opt$lpsi, metric, ref)
  J2 <- matrix(0, G, 2 * G)
  al <- opt$alpha * exp(-kap / 2)
  J2[cbind(seq_len(G), seq_len(G))] <- exp(-kap / 2)
  if (metric == "refgroup") J2[, G + ref] <- J2[, G + ref] - al / 2
  list(J1 = NULL, K = J2 %*% opt$map$J)
}

.jac_vars <- function(inp, opt, metric, ref) {
  G <- length(inp$groups)
  kap <- .kappa(opt$lpsi, metric, ref)
  psi <- exp(opt$lpsi - kap)
  J2 <- matrix(0, G, 2 * G)
  J2[cbind(seq_len(G), G + seq_len(G))] <- psi
  if (metric == "refgroup") J2[, G + ref] <- J2[, G + ref] - psi
  list(J1 = NULL, K = J2 %*% opt$map$J)
}

.jac_slope <- function(inp, opt, metric, ref, i) {
  G <- length(inp$groups)
  rows <- which(!is.na(inp$a0[, i]))
  kap <- .kappa(opt$lpsi, metric, ref)
  n <- length(rows)
  idx <- matrix(inp$ia[rows, i] + 1L, n, 1)
  sc <- exp(-(opt$lpsi[rows] - kap) / 2)
  ap <- inp$a0[rows, i] * sc
  J2 <- matrix(0, n, 2 * G)
  J2[cbind(seq_len(n), G + rows)] <- -ap / 2
  if (metric == "refgroup") J2[, G + ref] <- J2[, G + ref] + ap / 2
  list(J1 = list(idx = idx, val = matrix(sc, n, 1), grp = rows),
       K = J2 %*% opt$map$J, rows = rows)
}

.jac_intercept <- function(inp, opt, metric, ref, p) {
  G <- length(inp$groups)
  i <- inp$citem[p]
  rows <- which(!is.na(inp$c0[, p]))
  n <- length(rows)
  s <- exp(-opt$lpsi[rows] / 2)
  a0 <- inp$a0[rows, i]
  al <- opt$alpha[rows]
  idx <- cbind(inp$ic[rows, p] + 1L, inp$ia[rows, i] + 1L)
  val <- cbind(1, -al * s)
  J2 <- matrix(0, n, 2 * G)
  J2[cbind(seq_len(n), rows)] <- -a0 * s
  J2[cbind(seq_len(n), G + rows)] <- al * a0 * s / 2
  list(J1 = list(idx = idx, val = val, grp = rows), K = J2 %*% opt$map$J, rows = rows)
}
