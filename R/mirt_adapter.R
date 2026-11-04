# Conversion of single-group mirt fits to the standard form -----------------

.supported_itemtypes <- c("2PL", "3PL", "graded", "gpcm")

# Standard form of one fitted model: for each item one slope `a`, intercepts
# `c` (local boundary order), lower asymptote `g` and the covariance matrix of
# all standard-form estimates (NULL when mirt has no usable vcov).
.std_from_mirt <- function(fit, group) {
  if (!inherits(fit, "SingleGroupClass"))
    stop("Group ", group, ": expected a single-group model from mirt::mirt().",
         call. = FALSE)
  if (mirt::extract.mirt(fit, "nfact") != 1)
    stop("Group ", group, ": only unidimensional models are supported.",
         call. = FALSE)
  items <- mirt::extract.mirt(fit, "itemnames")
  itype <- mirt::extract.mirt(fit, "itemtype")
  bad <- setdiff(unique(itype), .supported_itemtypes)
  if (length(bad))
    stop("Group ", group, ": unsupported itemtype ", paste(bad, collapse = ", "),
         ". Supported: ", paste(.supported_itemtypes, collapse = ", "), ".",
         call. = FALSE)
  mv <- mirt::mod2values(fit)
  gp <- mv[mv$class == "GroupPars", , drop = FALSE]
  mu <- gp$value[gp$name == "MEAN_1"]
  s2 <- gp$value[gp$name == "COV_11"]
  if (any(gp$est) || abs(mu) > 1e-8 || abs(s2 - 1) > 1e-8)
    stop("Group ", group, ": the latent mean and variance must be fixed at 0 and 1.",
         call. = FALSE)

  mv$label <- paste0(mv$name, ".", mv$parnum)
  V <- tryCatch(mirt::extract.mirt(fit, "vcov"), error = function(e) NULL)
  if (!is.matrix(V) || !length(V)) V <- NULL

  out <- vector("list", length(items))
  names(out) <- items
  std_labs <- character()
  map_rows <- list()
  for (j in seq_along(items)) {
    it <- items[j]
    rows <- mv[mv$item == it, , drop = FALSE]
    val <- stats::setNames(rows$value, rows$name)
    lab <- stats::setNames(rows$label, rows$name)
    a <- unname(val["a1"])
    g <- NA_real_
    if (itype[j] %in% c("2PL", "3PL")) {
      cvec <- unname(val["d"])
      maps <- list(stats::setNames(1, lab["d"]))
      if (itype[j] == "3PL") g <- unname(val["g"])
    } else if (itype[j] == "graded") {
      dn <- grep("^d[0-9]+$", rows$name, value = TRUE)
      dn <- dn[order(as.integer(sub("^d", "", dn)))]
      cvec <- unname(val[dn])
      maps <- lapply(dn, function(nm) stats::setNames(1, lab[nm]))
    } else {
      akn <- grep("^ak[0-9]+$", rows$name, value = TRUE)
      akn <- akn[order(as.integer(sub("^ak", "", akn)))]
      if (!isTRUE(all.equal(unname(val[akn]), seq_along(akn) - 1)))
        stop("Group ", group, ", item ", it,
             ": gpcm items must use the default scoring 0, 1, 2, ...", call. = FALSE)
      dn <- grep("^d[0-9]+$", rows$name, value = TRUE)
      dn <- dn[order(as.integer(sub("^d", "", dn)))]
      d <- unname(val[dn])
      if (abs(d[1]) > 1e-10)
        stop("Group ", group, ", item ", it, ": gpcm intercept d0 must be 0.",
             call. = FALSE)
      cvec <- diff(d)
      maps <- lapply(seq_len(length(dn) - 1), function(k) {
        m <- stats::setNames(1, lab[dn[k + 1]])
        if (k > 1) m <- c(m, stats::setNames(-1, lab[dn[k]]))
        m
      })
    }
    out[[j]] <- list(item = it, itemtype = itype[j], a = a, c = cvec, g = g)
    std_labs <- c(std_labs, paste0(it, ":a"), paste0(it, ":c", seq_along(cvec)))
    map_rows <- c(map_rows, list(stats::setNames(1, lab["a1"])), maps)
  }

  Vstd <- NULL
  if (!is.null(V)) {
    # only the entries of the aligned parameters must be finite (a 3PL lower
    # asymptote at its bound may have none)
    used <- intersect(unique(unlist(lapply(map_rows, names))), colnames(V))
    V <- V[used, used, drop = FALSE]
    if (all(is.finite(V))) {
      Tm <- matrix(0, length(std_labs), ncol(V), dimnames = list(std_labs, colnames(V)))
      for (r in seq_along(map_rows)) {
        m <- map_rows[[r]]
        m <- m[names(m) %in% colnames(V)]
        if (length(m)) Tm[r, names(m)] <- m
      }
      Vstd <- Tm %*% V %*% t(Tm)
    }
  }
  list(items = out, vcov = Vstd, n = nrow(mirt::extract.mirt(fit, "data")),
       converged = isTRUE(mirt::extract.mirt(fit, "converged")))
}

# Combine group-wise standard forms into one input. `cats` optionally gives,
# for each group, a named list of observed category values per item.
.input_from_std <- function(std, groups, n, cats = NULL) {
  G <- length(std)
  items <- unique(unlist(lapply(std, function(s) names(s$items))))
  I <- length(items)
  itemtype <- stats::setNames(rep(NA_character_, I), items)
  K_local <- matrix(NA_integer_, G, I, dimnames = list(groups, items))
  for (g in seq_len(G)) {
    for (it in names(std[[g]]$items)) {
      s <- std[[g]]$items[[it]]
      if (!is.na(itemtype[it]) && itemtype[it] != s$itemtype)
        stop("Item ", it, " has different itemtypes across groups.", call. = FALSE)
      itemtype[it] <- s$itemtype
      K_local[g, it] <- length(s$c)
    }
  }

  # global category values per item
  gcats <- stats::setNames(vector("list", I), items)
  for (it in items) {
    if (!is.null(cats)) {
      gcats[[it]] <- sort(unique(unlist(lapply(cats, function(cg) cg[[it]]))))
    } else {
      k <- unique(stats::na.omit(K_local[, it]))
      if (length(k) != 1)
        stop("Item ", it, " has a different number of categories across groups. ",
             "Use align_irt() with the raw data so that categories can be matched.",
             call. = FALSE)
      gcats[[it]] <- seq_len(k + 1) - 1
    }
  }
  K <- vapply(gcats, function(x) length(x) - 1L, 1L)

  a0 <- matrix(NA_real_, G, I, dimnames = list(groups, items))
  guess <- a0
  c_list <- lapply(K, function(k) matrix(NA_real_, G, k, dimnames = list(groups, NULL)))
  vc <- vector("list", G)
  for (g in seq_len(G)) {
    relab <- character()
    for (it in names(std[[g]]$items)) {
      s <- std[[g]]$items[[it]]
      obs <- if (is.null(cats)) gcats[[it]] else sort(cats[[g]][[it]])
      if (length(obs) != length(s$c) + 1)
        stop("Group ", groups[g], ", item ", it,
             ": observed categories do not match the fitted model.", call. = FALSE)
      pos <- match(obs[-1], gcats[[it]][-1])
      if (s$itemtype == "gpcm") {
        adjacent <- match(obs[-length(obs)], gcats[[it]]) + 1L == match(obs[-1], gcats[[it]])
        if (!all(adjacent))
          stop("Group ", groups[g], ", item ", it, ": a gpcm item with an empty ",
               "middle category cannot be aligned. Collapse the category or use ",
               "the graded model.", call. = FALSE)
      }
      a0[g, it] <- s$a
      guess[g, it] <- s$g
      c_list[[it]][g, pos] <- s$c
      relab <- c(relab, stats::setNames(paste0(it, ":c", pos),
                                        paste0(it, ":c", seq_along(s$c))))
    }
    if (!is.null(std[[g]]$vcov)) {
      Vg <- std[[g]]$vcov
      nm <- rownames(Vg)
      hit <- nm %in% names(relab)
      nm[hit] <- relab[nm[hit]]
      dimnames(Vg) <- list(nm, nm)
      vc[[g]] <- Vg
    }
  }
  has_v <- all(vapply(vc, Negate(is.null), TRUE))
  if (all(is.na(guess))) guess <- NULL
  .new_input(a0, c_list, n = n, itemtype = unname(itemtype),
             vcov = if (has_v) vc else NULL, guess = guess, cats = gcats)
}
