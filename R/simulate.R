# Response simulation from aligned (standard-form) parameters --------------
#
# pars: list over items of list(type, a = G-vector, c = G x K matrix) on the
# latent scale eta ~ N(alpha_g, psi_g). Types: "2PL", "graded" (cumulative
# intercepts, decreasing) and "gpcm" (step intercepts). Responses are coded
# 0, 1, ..., K.

.simulate_responses <- function(pars, alpha, psi, n) {
  G <- length(alpha)
  out <- vector("list", G)
  for (g in seq_len(G)) {
    eta <- stats::rnorm(n[g], alpha[g], sqrt(psi[g]))
    resp <- matrix(NA_integer_, n[g], length(pars))
    for (j in seq_along(pars)) {
      it <- pars[[j]]
      a <- it$a[g]
      cc <- matrix(it$c, nrow = G)[g, ]
      if (is.na(a)) next
      if (it$type == "2PL") {
        resp[, j] <- as.integer(stats::runif(n[g]) < stats::plogis(a * eta + cc))
      } else if (it$type == "graded") {
        cum <- sapply(cc, function(ck) stats::plogis(a * eta + ck))
        cum <- matrix(cum, nrow = n[g])
        u <- stats::runif(n[g])
        resp[, j] <- rowSums(u < cum)
      } else if (it$type == "gpcm") {
        K <- length(cc)
        num <- sapply(0:K, function(k) k * a * eta + c(0, cumsum(cc))[k + 1])
        num <- matrix(num, nrow = n[g])
        num <- exp(num - apply(num, 1, max))
        pr <- num / rowSums(num)
        u <- stats::runif(n[g])
        resp[, j] <- rowSums(u > t(apply(pr, 1, cumsum))[, seq_len(K), drop = FALSE])
      } else {
        stop("Unknown item type: ", it$type, call. = FALSE)
      }
    }
    colnames(resp) <- names(pars) %||% paste0("I", seq_along(pars))
    out[[g]] <- resp
  }
  out
}
