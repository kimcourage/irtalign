# Noninvariance detection and summaries ------------------------------------

.pairwise_p <- function(x, S) {
  dv <- outer(diag(S), diag(S), "+") - 2 * S
  dx <- outer(x, x, "-")
  z <- abs(dx) / sqrt(pmax(dv, 0))
  z[dv <= 0 & dx == 0] <- 0
  z[dv <= 0 & dx != 0] <- Inf
  p <- 2 * stats::pnorm(-z)
  diag(p) <- 1
  p
}

.components <- function(adj) {
  m <- nrow(adj)
  comp <- rep(NA_integer_, m)
  k <- 0L
  for (s in seq_len(m)) {
    if (!is.na(comp[s])) next
    k <- k + 1L
    queue <- s
    comp[s] <- k
    while (length(queue)) {
      v <- queue[1]
      queue <- queue[-1]
      nb <- which(adj[v, ] & is.na(comp))
      comp[nb] <- k
      queue <- c(queue, nb)
    }
  }
  comp
}

# Invariant set for one item parameter (Asparouhov & Muthen, 2014): groups are
# connected when their pairwise difference has p > p_connect; the largest
# connected set starts the search. Groups are then added or removed by
# comparing each group with the sample-size weighted mean of the current set
# (p > p_member) until the set no longer changes.
.detect_one <- function(x, S, n, p_connect, p_member, maxit = 100L) {
  m <- length(x)
  if (m == 0)
    return(list(invariant = logical(0), mean = NA_real_, iterations = 0L, cycle = FALSE))
  if (m < 2)
    return(list(invariant = TRUE, mean = x, iterations = 0L, cycle = FALSE))
  pv <- .pairwise_p(x, S)
  comp <- .components(pv > p_connect)
  sizes <- tabulate(comp)
  cand <- which(sizes == max(sizes))
  if (length(cand) > 1) {
    tot <- vapply(cand, function(k) sum(n[comp == k]), 1)
    cand <- cand[which.max(tot)]
  }
  set <- which(comp == cand)
  seen <- list()
  cycle <- FALSE
  iter <- 0L
  repeat {
    iter <- iter + 1L
    w <- numeric(m)
    w[set] <- n[set] / sum(n[set])
    mbar <- sum(w * x)
    C <- diag(m) - matrix(w, m, m, byrow = TRUE)
    v <- rowSums((C %*% S) * C)
    d <- x - mbar
    z <- abs(d) / sqrt(pmax(v, 0))
    z[v <= 0 & d == 0] <- 0
    z[v <= 0 & d != 0] <- Inf
    p <- 2 * stats::pnorm(-z)
    new <- which(p > p_member)
    if (setequal(new, set) || iter >= maxit) break
    key <- paste(new, collapse = ",")
    if (key %in% seen || !length(new)) {
      cycle <- length(new) > 0
      if (length(new)) set <- new
      break
    }
    seen <- c(seen, paste(set, collapse = ","))
    set <- new
  }
  inv <- seq_len(m) %in% set
  w <- n[set] / sum(n[set])
  list(invariant = inv, mean = if (length(set)) sum(w * x[set]) else NA_real_,
       iterations = iter, cycle = cycle)
}

# R-squared of invariance (Asparouhov & Muthen, 2014, eqs. 13-14), truncated
# at zero. Means are unweighted over the groups that have the parameter.
.r2_slope <- function(a0, a_rep, psi) {
  if (length(a0) < 2 || stats::var(a0) <= 0) return(NA_real_)
  max(0, 1 - stats::var(a0 - sqrt(psi) * mean(a_rep)) / stats::var(a0))
}

.r2_intercept <- function(c0, c_rep, alpha, a_rep_mean) {
  if (length(c0) < 2 || stats::var(c0) <= 0) return(NA_real_)
  max(0, 1 - stats::var(c0 - mean(c_rep) - alpha * a_rep_mean) / stats::var(c0))
}

.mean_comparison <- function(groups, alpha, S, level) {
  ord <- order(alpha, decreasing = TRUE)
  pv <- .pairwise_p(alpha, S)
  smaller <- lapply(ord, function(g) {
    k <- ord[alpha[ord] < alpha[g] & pv[g, ord] < level]
    groups[k]
  })
  out <- data.frame(rank = seq_along(ord), group = groups[ord], mean = alpha[ord],
                    stringsAsFactors = FALSE)
  out$smaller <- smaller
  out
}
