# Standard-form input shared by all entry points ---------------------------
#
# a0      G x I configural slopes (NA where an item is absent in a group)
# c0      G x P configural intercepts, items stacked: item i has K_i columns
# citem   item index (1-based) of each intercept column
# ck      boundary index of each intercept column within its item
# vcov    list of G covariance matrices of the group's standard-form
#         estimates, ordered item by item (slope, then intercepts); NULL if
#         unavailable
# ia, ic  0-based positions of each estimate in the stacked parameter vector
#         m1 (-1 when absent); m1_index lists the 1-based positions per group

.new_input <- function(a0, c_list, n, itemtype, vcov = NULL, guess = NULL,
                       cats = NULL) {
  groups <- rownames(a0)
  items <- colnames(a0)
  G <- nrow(a0)
  I <- ncol(a0)
  if (G < 2) stop("At least two groups are required.", call. = FALSE)
  if (length(n) != G || any(!is.finite(n)) || any(n <= 0))
    stop("`n` must give a positive sample size for each group.", call. = FALSE)
  if (!identical(names(c_list), items))
    stop("Intercepts must be given for the same items as the slopes.", call. = FALSE)

  # boundaries keep their position k (label "item:ck"); a boundary that no
  # group has is dropped
  c_list <- lapply(seq_along(c_list), function(i) {
    m <- c_list[[i]]
    m[is.na(a0[, i]), ] <- NA_real_
    colnames(m) <- seq_len(ncol(m))
    m[, colSums(!is.na(m)) > 0, drop = FALSE]
  })
  K <- vapply(c_list, ncol, 1L)
  if (any(K == 0))
    stop("Item ", items[K == 0][1], " has no intercepts.", call. = FALSE)
  ck <- as.integer(unlist(lapply(c_list, colnames)))
  c0 <- do.call(cbind, c_list)
  rownames(c0) <- groups
  colnames(c0) <- paste0(rep(items, K), ":c", ck)
  citem <- rep(seq_len(I), K)

  for (i in seq_len(I)) {
    absent <- is.na(a0[, i])
    cols <- which(citem == i)
    c0[absent, cols] <- NA_real_
    if (any(!absent & rowSums(!is.na(c0[, cols, drop = FALSE])) == 0))
      stop("Item ", items[i], " has a slope but no intercept in some group.",
           call. = FALSE)
  }
  present_groups <- colSums(!is.na(a0))
  if (any(present_groups < 2))
    stop("Each item must be present in at least two groups: ",
         paste(items[present_groups < 2], collapse = ", "), call. = FALSE)

  inp <- list(groups = groups, n = as.numeric(n), items = items,
              itemtype = itemtype, K = K, a0 = a0, c0 = c0, citem = citem,
              ck = ck, guess = guess, cats = cats, vcov = NULL)
  if (!is.null(vcov)) inp <- .attach_vcov(inp, vcov)
  inp
}

.attach_vcov <- function(inp, vcov) {
  G <- length(inp$groups)
  if (!is.list(vcov) || length(vcov) != G)
    stop("`vcov` must be a list with one covariance matrix per group.", call. = FALSE)
  ia <- matrix(-1L, G, length(inp$items))
  ic <- matrix(-1L, G, ncol(inp$c0))
  m1_index <- vector("list", G)
  V <- vector("list", G)
  off <- 0L
  for (g in seq_len(G)) {
    labs <- character()
    for (i in seq_along(inp$items)) {
      if (is.na(inp$a0[g, i])) next
      cols <- which(inp$citem == i & !is.na(inp$c0[g, ]))
      labs <- c(labs, paste0(inp$items[i], ":a"), colnames(inp$c0)[cols])
    }
    Vg <- vcov[[g]]
    if (is.null(Vg) || !is.matrix(Vg) || is.null(rownames(Vg)))
      stop("Covariance matrix for group ", inp$groups[g],
           " must be a matrix with parameter labels.", call. = FALSE)
    miss <- setdiff(labs, rownames(Vg))
    if (length(miss))
      stop("Covariance matrix for group ", inp$groups[g], " lacks: ",
           paste(miss, collapse = ", "), call. = FALSE)
    Vg <- Vg[labs, labs, drop = FALSE]
    if (any(!is.finite(Vg)))
      stop("Covariance matrix for group ", inp$groups[g],
           " has missing or infinite entries.", call. = FALSE)
    V[[g]] <- (Vg + t(Vg)) / 2
    pos <- off + seq_along(labs) - 1L
    names(pos) <- labs
    for (i in seq_along(inp$items)) {
      lab <- paste0(inp$items[i], ":a")
      if (lab %in% labs) ia[g, i] <- pos[[lab]]
    }
    for (p in seq_len(ncol(inp$c0))) {
      lab <- colnames(inp$c0)[p]
      if (lab %in% labs) ic[g, p] <- pos[[lab]]
    }
    m1_index[[g]] <- pos + 1L
    off <- off + length(labs)
  }
  inp$vcov <- V
  inp$ia <- ia
  inp$ic <- ic
  inp$m1_index <- m1_index
  inp$M1 <- off
  inp
}
