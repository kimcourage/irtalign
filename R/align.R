#' Alignment optimization for multiple-group item response models
#'
#' `align_irt()` fits a unidimensional item response model separately in each
#' group with [mirt::mirt()] (latent mean 0 and variance 1 in every group) and
#' aligns the group-specific estimates. `align_fits()` aligns models that have
#' already been fitted, and `align_pars()` aligns parameter estimates supplied
#' directly.
#'
#' Every item is written in a standard form with one slope `a` and intercepts
#' `c1, ..., cK` on the logit scale: `a * theta + c` for 2PL and 3PL items,
#' `a * theta + ck` for the cumulative logits of graded items, and the step
#' intercepts `dk - d(k-1)` for generalized partial credit items. For a group
#' with latent mean \eqn{\alpha_g} and variance \eqn{\psi_g} the aligned
#' parameters are
#' \deqn{a_{g} = a_{g0} / \sqrt{\psi_g}, \qquad c_{g} = c_{g0} - \alpha_g a_{g0} / \sqrt{\psi_g},}
#' where \eqn{a_{g0}} and \eqn{c_{g0}} are the configural estimates. The means
#' and variances minimize the simplicity function
#' \deqn{F = \sum_{\textrm{parameters}} \sum_{g < h} w_{gh} f(x_g - x_h), \qquad f(x) = (x^2 + \epsilon)^{1/4},}
#' over all slopes and intercepts, with \eqn{w_{gh} = w_g w_h} and
#' \eqn{w_g = \sqrt{N_g} / \textrm{mean}(\sqrt{N})} (Asparouhov & Muthen, 2014).
#' The optimization uses the product of the group variances equal to one;
#' `align = "fixed"` also fixes the mean of the reference group at zero, while
#' `align = "free"` estimates all means. Results are reported with the
#' reference group variance equal to one (`metric = "refgroup"`) or with the
#' product of the variances equal to one (`metric = "product"`).
#'
#' Standard errors of the means, variances and aligned parameters follow from
#' the implicit function theorem applied to the first-order conditions of the
#' simplicity function and the delta method. For each slope and intercept the
#' invariant set of groups is found as in Asparouhov and Muthen (2014): groups
#' whose pairwise difference has a p-value above `p_connect` are connected,
#' the largest connected set is the starting set, and each group is then kept
#' in the set when its difference from the sample-size weighted mean of the set
#' has a p-value above `p_member`, until the set does not change. Groups
#' outside the set are noninvariant for that parameter.
#'
#' The lower asymptote of 3PL items is not aligned and is reported in
#' `$guess`. In `align_irt()` an item that shows a single response category in
#' a group is left out of that group. Graded items may have empty categories
#' in some groups; the intercepts of the missing boundaries are omitted for
#' those groups.
#'
#' @param data A data frame or matrix of item responses.
#' @param items Column names of the items. By default all columns except
#'   `group`.
#' @param group Name of the grouping column in `data`, or a vector of group
#'   labels with one element per row.
#' @param itemtype Item types passed to [mirt::mirt()]: `"2PL"`, `"3PL"`,
#'   `"graded"` or `"gpcm"`, either one value for all items or one per item. By
#'   default items with two categories are `"2PL"` and others `"graded"`.
#' @param align `"fixed"` (mean of the reference group fixed at zero) or
#'   `"free"`.
#' @param ref Reference group as a label or index, or `"auto"` (default) for
#'   the group whose mean is closest to zero in a free alignment; the free
#'   solution is also used as an additional starting value. With two groups
#'   `"auto"` selects the first group.
#' @param metric `"refgroup"` or `"product"`; see Details.
#' @param simplicity `"sqrt"` for \eqn{f(x) = (x^2 + \epsilon)^{1/4}} or
#'   `"fourthrt"` for \eqn{f(x) = (x^2 + \epsilon)^{1/8}}.
#' @param eps The constant \eqn{\epsilon} of the simplicity function.
#' @param nstarts Number of starting values. The first two are a least-squares
#'   linking solution and zero means with unit variances; the others are
#'   random. Each start is optimized by Newton steps with the analytic Hessian
#'   ([stats::nlminb()]).
#' @param item_weight `"parameter"` gives every slope and intercept weight one;
#'   `"item"` divides the intercept weights of an item by its number of
#'   intercepts.
#' @param se Compute standard errors and invariance tests.
#' @param p_connect,p_member P-value cutoffs of the invariance search.
#' @param level Significance level of the factor mean comparisons.
#' @param seed Optional seed for the random starting values. The global random
#'   number state is restored afterwards.
#' @param control List of control settings passed to [stats::nlminb()].
#' @param ... Further arguments passed to [mirt::mirt()], for example
#'   `quadpts` or `SE.type`.
#'
#' @return An object of class `irtalign` with elements
#'   * `means`: group sample sizes, aligned means and variances with standard
#'     errors;
#'   * `params`: configural and aligned slopes and intercepts with standard
#'     errors and invariance status by group;
#'   * `item_summary`: for each slope and intercept the number of
#'     noninvariant groups, the weighted mean over the invariant set,
#'     invariance R-squared, contribution to the simplicity function and the
#'     variance of the aligned values;
#'   * `mean_comparison`: groups ranked by mean with the groups that have a
#'     significantly smaller mean;
#'   * `optimization`: the minimum of the simplicity function, the number of
#'     starts that reached it and the values from all starts;
#'   * `settings`, `group_contribution`, `vcov_means`, `guess`, `input` and
#'     `fits`.
#'
#' @references
#' Asparouhov, T., & Muthen, B. (2014). Multiple-group factor analysis
#' alignment. *Structural Equation Modeling, 21*(4), 495-508.
#' \doi{10.1080/10705511.2014.919210}
#'
#' Muthen, B., & Asparouhov, T. (2014). IRT studies of many groups: The
#' alignment method. *Frontiers in Psychology, 5*, 978.
#' \doi{10.3389/fpsyg.2014.00978}
#'
#' @seealso [summary.irtalign()], [invariance_table()], [mean_comparison()]
#'
#' @examples
#' \donttest{
#' dat <- mirt::Science
#' dat$group <- rep(c("A", "B", "C"), length.out = nrow(dat))
#' fit <- align_irt(dat, group = "group", nstarts = 10, seed = 1)
#' fit
#' summary(fit)
#' }
#' @export
align_irt <- function(data, items = NULL, group, itemtype = NULL,
                      align = c("fixed", "free"), ref = "auto",
                      metric = c("refgroup", "product"),
                      simplicity = c("sqrt", "fourthrt"), eps = 1e-4, nstarts = 30,
                      item_weight = c("parameter", "item"), se = TRUE,
                      p_connect = 0.01, p_member = 0.001, level = 0.05,
                      seed = NULL, control = list(), ...) {
  args <- .align_args(align, ref, metric, simplicity, eps, nstarts, item_weight,
                      se, p_connect, p_member, level, seed, control)
  data <- as.data.frame(data)
  if (is.character(group) && length(group) == 1 && group %in% names(data)) {
    gv <- data[[group]]
    if (is.null(items)) items <- setdiff(names(data), group)
  } else {
    if (length(group) != nrow(data))
      stop("`group` must be a column name or have one value per row.", call. = FALSE)
    gv <- group
    if (is.null(items)) items <- names(data)
  }
  if (anyNA(gv)) stop("`group` has missing values.", call. = FALSE)
  miss <- setdiff(items, names(data))
  if (length(miss)) stop("Unknown items: ", paste(miss, collapse = ", "), call. = FALSE)
  resp <- data[, items, drop = FALSE]
  if (!all(vapply(resp, is.numeric, TRUE)))
    stop("Item responses must be numeric.", call. = FALSE)

  groups <- if (is.factor(gv)) levels(droplevels(gv)) else sort(unique(as.character(gv)))
  gv <- as.character(gv)
  ncat <- vapply(resp, function(x) length(unique(stats::na.omit(x))), 1L)
  if (is.null(itemtype)) {
    itemtype <- ifelse(ncat <= 2, "2PL", "graded")
  } else if (length(itemtype) == 1) {
    itemtype <- rep(itemtype, length(items))
  } else if (!is.null(names(itemtype))) {
    itemtype <- unname(itemtype[items])
  }
  if (length(itemtype) != length(items) || anyNA(itemtype))
    stop("`itemtype` must have one value or one value per item.", call. = FALSE)
  bad <- setdiff(unique(itemtype), .supported_itemtypes)
  if (length(bad))
    stop("Unsupported itemtype: ", paste(bad, collapse = ", "), ". Supported: ",
         paste(.supported_itemtypes, collapse = ", "), ".", call. = FALSE)
  names(itemtype) <- items

  fits <- vector("list", length(groups))
  names(fits) <- groups
  cats <- fits
  n <- numeric(length(groups))
  dropped <- character()
  for (k in seq_along(groups)) {
    g <- groups[k]
    dg <- resp[gv == g, , drop = FALSE]
    dg <- dg[rowSums(!is.na(dg)) > 0, , drop = FALSE]
    keep <- vapply(dg, function(x) length(unique(stats::na.omit(x))) >= 2, TRUE)
    if (!any(keep)) stop("Group ", g, " has no item with two or more categories.",
                         call. = FALSE)
    if (!all(keep)) dropped <- c(dropped, paste0(g, ":", items[!keep]))
    dg <- dg[, keep, drop = FALSE]
    dg <- dg[rowSums(!is.na(dg)) > 0, , drop = FALSE]
    n[k] <- nrow(dg)
    cats[[k]] <- lapply(dg, function(x) sort(unique(stats::na.omit(x))))
    fits[[k]] <- tryCatch(
      mirt::mirt(dg, 1, itemtype = unname(itemtype[keep]), SE = args$se,
                 verbose = FALSE, ...),
      error = function(e) stop("Group ", g, ": ", conditionMessage(e), call. = FALSE))
  }
  if (length(dropped))
    message("Items without response variation were left out of these groups: ",
            paste(dropped, collapse = ", "))
  .align_from_fits(fits, n, cats, args, match.call())
}

#' @rdname align_irt
#' @param fits A named list of single-group, unidimensional models fitted with
#'   [mirt::mirt()] with the latent mean and variance fixed at 0 and 1. Names
#'   are used as group labels.
#' @param n Sample size of each group, matched by name when named. By default
#'   the number of rows used in each fit.
#' @export
align_fits <- function(fits, n = NULL, align = c("fixed", "free"), ref = "auto",
                       metric = c("refgroup", "product"),
                       simplicity = c("sqrt", "fourthrt"), eps = 1e-4, nstarts = 30,
                       item_weight = c("parameter", "item"), se = TRUE,
                       p_connect = 0.01, p_member = 0.001, level = 0.05,
                       seed = NULL, control = list()) {
  args <- .align_args(align, ref, metric, simplicity, eps, nstarts, item_weight,
                      se, p_connect, p_member, level, seed, control)
  if (!is.list(fits) || length(fits) < 2)
    stop("`fits` must be a list of at least two fitted models.", call. = FALSE)
  if (is.null(names(fits))) names(fits) <- paste0("G", seq_along(fits))
  n <- .by_group(n, names(fits), "n")
  # category values of each item as stored with the fitted data
  cats <- lapply(names(fits), function(g) {
    d <- tryCatch(mirt::extract.mirt(fits[[g]], "data"), error = function(e) NULL)
    if (is.null(d)) return(NULL)
    lapply(as.data.frame(d), function(x) sort(unique(stats::na.omit(x))))
  })
  if (any(vapply(cats, is.null, TRUE))) cats <- NULL
  .align_from_fits(fits, n, cats, args, match.call())
}

# Put a per-group vector or list in the order of `groups`: by name when named,
# otherwise by position.
.by_group <- function(x, groups, arg) {
  if (is.null(x) || is.null(names(x))) return(x)
  if (!setequal(names(x), groups) || anyDuplicated(names(x)))
    stop("Names of `", arg, "` do not match the groups.", call. = FALSE)
  x[groups]
}

.align_from_fits <- function(fits, n, cats, args, call) {
  groups <- names(fits)
  std <- lapply(seq_along(fits), function(k) .std_from_mirt(fits[[k]], groups[k]))
  if (is.null(n)) n <- vapply(std, `[[`, 1, "n")
  nc <- !vapply(std, `[[`, TRUE, "converged")
  if (any(nc))
    warning("The model did not converge in group(s): ", paste(groups[nc], collapse = ", "),
            call. = FALSE)
  if (args$se && any(vapply(std, function(s) is.null(s$vcov), TRUE))) {
    warning("Standard errors are not available for every group; invariance tests ",
            "are skipped.", call. = FALSE)
  }
  inp <- .input_from_std(std, groups, n, cats)
  .align_core(inp, args, fits = fits, call = call)
}

#' @rdname align_irt
#' @param a Matrix of configural slopes with groups in rows and items in
#'   columns (`NA` for an item that is absent in a group).
#' @param c Configural intercepts: a matrix shaped like `a` for items with one
#'   intercept, or a list with one matrix of intercepts (groups by boundaries)
#'   per item, named by item.
#' @param vcov Optional list with one covariance matrix per group of that
#'   group's estimates, labelled `"item:a"` and `"item:c1"`, `"item:c2"`, ...
#'   Named elements of `n` and `vcov`, and row names of `c`, are matched to the
#'   row names of `a`; otherwise groups are taken in order.
#' @export
align_pars <- function(a, c, n, vcov = NULL, align = c("fixed", "free"), ref = "auto",
                       metric = c("refgroup", "product"),
                       simplicity = c("sqrt", "fourthrt"), eps = 1e-4, nstarts = 30,
                       item_weight = c("parameter", "item"), se = TRUE,
                       p_connect = 0.01, p_member = 0.001, level = 0.05,
                       seed = NULL, control = list()) {
  args <- .align_args(align, ref, metric, simplicity, eps, nstarts, item_weight,
                      se, p_connect, p_member, level, seed, control)
  a <- as.matrix(a)
  if (is.null(rownames(a))) rownames(a) <- paste0("G", seq_len(nrow(a)))
  if (is.null(colnames(a))) colnames(a) <- paste0("I", seq_len(ncol(a)))
  if (is.matrix(c) || is.data.frame(c)) {
    c <- as.matrix(c)
    if (!all(dim(c) == dim(a)))
      stop("`c` must have the same dimensions as `a`.", call. = FALSE)
    c <- lapply(seq_len(ncol(a)), function(i) c[, i, drop = FALSE])
    names(c) <- colnames(a)
  }
  if (!is.list(c) || !setequal(names(c), colnames(a)))
    stop("`c` must be a matrix or a list named by item.", call. = FALSE)
  c <- lapply(c[colnames(a)], function(m) {
    m <- as.matrix(m)
    if (nrow(m) != nrow(a)) stop("Each intercept matrix needs one row per group.",
                                 call. = FALSE)
    if (!is.null(rownames(m))) {
      if (!setequal(rownames(m), rownames(a)))
        stop("Row names of the intercepts do not match the groups.", call. = FALSE)
      m <- m[rownames(a), , drop = FALSE]
    }
    dimnames(m) <- list(rownames(a), NULL)
    m
  })
  n <- .by_group(n, rownames(a), "n")
  vcov <- .by_group(vcov, rownames(a), "vcov")
  inp <- .new_input(a, c, n = n, itemtype = rep("custom", ncol(a)), vcov = vcov)
  .align_core(inp, args, call = match.call())
}
