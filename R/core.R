# Assemble an alignment analysis from a standard-form input ---------------

.align_args <- function(align, ref, metric, simplicity, eps, nstarts, item_weight,
                        se, p_connect, p_member, level, seed, control) {
  list(align = match.arg(align, c("fixed", "free")),
       ref = ref,
       metric = match.arg(metric, c("refgroup", "product")),
       simplicity = match.arg(simplicity, c("sqrt", "fourthrt")),
       eps = eps, nstarts = as.integer(nstarts),
       item_weight = match.arg(item_weight, c("parameter", "item")),
       se = isTRUE(se), p_connect = p_connect, p_member = p_member,
       level = level, seed = seed, control = control)
}

.align_core <- function(inp, args, fits = NULL, call = NULL) {
  G <- length(inp$groups)
  I <- length(inp$items)
  P <- ncol(inp$c0)
  power <- .simplicity_power(args$simplicity)
  wa <- rep(1, I)
  wc <- if (args$item_weight == "item") 1 / inp$K[inp$citem] else rep(1, P)
  if (!is.numeric(args$eps) || length(args$eps) != 1 || args$eps <= 0)
    stop("`eps` must be a positive number.", call. = FALSE)
  if (args$nstarts < 1) stop("`nstarts` must be at least 1.", call. = FALSE)

  # ref = "auto": the group whose mean is closest to zero in a free alignment.
  # Its free solution, shifted, is added to the starting values.
  extra <- NULL
  ref_auto <- identical(args$ref, "auto")
  if (ref_auto && G > 2) {
    free <- .optimize_alignment(inp, "free", 1L, args$eps, power, args$nstarts,
                                wa, wc, args$seed, args$control)
    ref <- which.min(abs(free$alpha))
    extra <- list(alpha = free$alpha - free$alpha[ref], lpsi = free$lpsi)
  } else if (ref_auto) {
    ref <- 1L
  } else {
    ref <- .resolve_group(args$ref, inp$groups)
  }
  if (args$align == "free" && G == 2)
    warning("FREE alignment is not identified with two groups; use align = \"fixed\".",
            call. = FALSE)

  opt <- .optimize_alignment(inp, args$align, ref, args$eps, power, args$nstarts,
                             wa, wc, args$seed, args$control,
                             extra_start = if (args$align == "fixed") extra)
  rp <- .reported(inp, opt, args$metric, ref)

  vc <- NULL
  if (args$se && !is.null(inp$vcov)) {
    vc <- .align_vcov(inp, opt, wa, wc, args$eps, power)
    if (!vc$ok) {
      warning("The simplicity function is not locally convex at the solution; ",
              "standard errors and invariance tests are not available.", call. = FALSE)
      vc <- NULL
    }
  }
  if (!is.null(vc) && args$align == "free" && vc$cond > 1e8)
    warning("FREE alignment may be poorly identified; consider align = \"fixed\" ",
            "with ref = \"auto\".", call. = FALSE)

  # group means and variances
  se_mean <- se_var <- rep(NA_real_, G)
  Smean <- NULL
  if (!is.null(vc)) {
    jm <- .jac_means(inp, opt, args$metric, ref)
    Smean <- .qcov(NULL, jm$K, vc, inp)
    jv <- .jac_vars(inp, opt, args$metric, ref)
    Svar <- .qcov(NULL, jv$K, vc, inp)
    se_mean <- sqrt(pmax(diag(Smean), 0))
    se_var <- sqrt(pmax(diag(Svar), 0))
    dimnames(Smean) <- list(inp$groups, inp$groups)
  }
  means <- data.frame(group = inp$groups, n = inp$n, mean = rp$alpha, se_mean = se_mean,
                      variance = rp$psi, se_variance = se_var, stringsAsFactors = FALSE)

  # aligned item parameters, invariance and summaries
  slope1 <- inp$a0 * exp(-opt$lpsi / 2)
  int1 <- inp$c0 - (opt$alpha * slope1)[, inp$citem, drop = FALSE]
  contrib <- align_contrib_cpp(cbind(slope1, int1), opt$W, c(wa, wc), args$eps, power)
  long <- list()
  summ <- list()
  for (q in seq_len(I + P)) {
    if (q <= I) {
      i <- q
      rows <- which(!is.na(inp$a0[, i]))
      param <- "a"
      conf <- inp$a0[rows, i]
      est <- rp$a[rows, i]
      r2 <- .r2_slope(conf, est, rp$psi[rows])
      jac <- if (!is.null(vc)) .jac_slope(inp, opt, args$metric, ref, i)
    } else {
      p <- q - I
      i <- inp$citem[p]
      rows <- which(!is.na(inp$c0[, p]))
      param <- paste0("c", inp$ck[p])
      conf <- inp$c0[rows, p]
      est <- rp$c[rows, p]
      r2 <- .r2_intercept(conf, est, rp$alpha[rows], mean(rp$a[!is.na(rp$a[, i]), i]))
      jac <- if (!is.null(vc)) .jac_intercept(inp, opt, args$metric, ref, p)
    }
    se <- rep(NA_real_, length(rows))
    inv <- rep(NA, length(rows))
    inv_mean <- NA_real_
    if (!is.null(vc)) {
      S <- .qcov(jac$J1, jac$K, vc, inp)
      se <- sqrt(pmax(diag(S), 0))
      det <- .detect_one(est, S, inp$n[rows], args$p_connect, args$p_member)
      inv <- det$invariant
      inv_mean <- det$mean
    }
    long[[q]] <- data.frame(item = inp$items[i], param = param, group = inp$groups[rows],
                            configural = conf, aligned = est, se = se, invariant = inv,
                            stringsAsFactors = FALSE)
    summ[[q]] <- data.frame(item = inp$items[i], param = param, groups = length(rows),
                            noninvariant = if (is.null(vc)) NA_integer_ else sum(!inv),
                            invariant_mean = inv_mean, r2 = r2,
                            contribution = sum(contrib[, q]),
                            variance = if (length(est) > 1) stats::var(est) else NA_real_,
                            stringsAsFactors = FALSE)
  }
  params <- do.call(rbind, long)
  item_summary <- do.call(rbind, summ)
  rownames(params) <- rownames(item_summary) <- NULL
  ord <- order(match(params$item, inp$items), params$param != "a", params$param)
  params <- params[ord, ]
  item_summary <- item_summary[order(match(item_summary$item, inp$items),
                                     item_summary$param != "a", item_summary$param), ]
  rownames(params) <- rownames(item_summary) <- NULL

  mc <- if (!is.null(Smean)) .mean_comparison(inp$groups, rp$alpha, Smean, args$level)

  structure(list(
    call = call,
    settings = list(align = args$align, ref = inp$groups[ref], ref_index = ref,
                    ref_auto = ref_auto,
                    metric = args$metric, simplicity = args$simplicity, eps = args$eps,
                    nstarts = args$nstarts, item_weight = args$item_weight,
                    p_connect = args$p_connect, p_member = args$p_member,
                    level = args$level),
    groups = inp$groups, items = inp$items, itemtype = inp$itemtype,
    means = means, params = params, item_summary = item_summary,
    group_contribution = data.frame(group = inp$groups, contribution = rowSums(contrib),
                                    stringsAsFactors = FALSE),
    mean_comparison = mc, vcov_means = Smean,
    optimization = list(value = opt$value, n_best = opt$n_best, starts = opt$starts,
                        converged = opt$converged, theta = opt$theta,
                        alpha = opt$alpha, lpsi = opt$lpsi,
                        condition = if (!is.null(vc)) vc$cond else NA_real_),
    guess = inp$guess, input = inp, fits = fits),
    class = "irtalign")
}
