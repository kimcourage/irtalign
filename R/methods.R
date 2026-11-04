#' Summaries of an alignment analysis
#'
#' `summary()` collects the invariance results, group means and variances,
#' the factor mean comparison and the item parameter summary. Printing it
#' shows the same blocks as the alignment output of Mplus: groups that are
#' noninvariant for a parameter are shown in parentheses.
#'
#' @param object,x An object of class `irtalign` (or `summary.irtalign` for
#'   printing).
#' @param digits Number of decimals shown.
#' @param width Line width for lists of groups.
#' @param ... Not used.
#'
#' @return `summary()` returns an object of class `summary.irtalign`;
#'   `coef()` returns the table of group means and variances.
#' @export
summary.irtalign <- function(object, ...) {
  s <- object$item_summary
  tested <- !is.na(s$noninvariant)
  pct <- if (any(tested)) {
    100 * sum(s$noninvariant[tested]) / sum(s$groups[tested])
  } else NA_real_
  structure(list(settings = object$settings, optimization = object$optimization,
                 groups = object$groups, n = object$means$n,
                 means = object$means, params = object$params, items = s,
                 mean_comparison = object$mean_comparison,
                 noninvariant_percent = pct),
            class = "summary.irtalign")
}

.fmt <- function(x, digits) formatC(x, format = "f", digits = digits)

.settings_lines <- function(st, opt, G, I, N) {
  how <- if (isTRUE(st$ref_auto)) " (automatic)" else ""
  ref <- if (st$align == "fixed") paste0(", reference group ", st$ref, how) else
    paste0(", metric reference group ", st$ref, how)
  metric <- if (st$metric == "refgroup") "reference group variance 1" else
    "product of variances 1"
  c(sprintf("Alignment: %s%s", toupper(st$align), ref),
    sprintf("Groups: %d   Items: %d   Observations: %s", G, I, format(sum(N), big.mark = ",")),
    sprintf("Simplicity function: %s, eps = %s; metric: %s",
            toupper(st$simplicity), format(st$eps), metric),
    sprintf("Simplicity function value: %s (reached by %d of %d starts)",
            .fmt(opt$value, 4), opt$n_best, nrow(opt$starts)))
}

#' @rdname summary.irtalign
#' @export
print.summary.irtalign <- function(x, digits = 3, width = 72, ...) {
  cat(.settings_lines(x$settings, x$optimization, length(x$groups), length(unique(x$items$item)),
                      x$n), sep = "\n")

  if (!all(is.na(x$params$invariant))) {
    cat("\nApproximate measurement (non)invariance for groups\n")
    for (blk in c("a", "c")) {
      sel <- if (blk == "a") x$items$param == "a" else x$items$param != "a"
      if (!any(sel)) next
      cat(if (blk == "a") " Slopes\n" else " Intercepts\n")
      for (r in which(sel)) {
        it <- x$items$item[r]
        pa <- x$items$param[r]
        pr <- x$params[x$params$item == it & x$params$param == pa, ]
        lab <- ifelse(pr$invariant, pr$group, paste0("(", pr$group, ")"))
        name <- if (blk == "a") it else paste(it, pa)
        lines <- .wrap_labels(lab, width - 18, indent = 18)
        cat(sprintf("   %-15s %s\n", name, lines[1]))
        if (length(lines) > 1) cat(paste0(lines[-1], "\n"), sep = "")
      }
    }
    cat(sprintf("\nNoninvariant parameter-group combinations: %s%%\n",
                .fmt(x$noninvariant_percent, 1)))
  }

  cat("\nGroup means and variances\n")
  m <- x$means
  tab <- data.frame(Group = m$group, N = m$n, Mean = .fmt(m$mean, digits),
                    SE = .fmt(m$se_mean, digits), Variance = .fmt(m$variance, digits),
                    `SE ` = .fmt(m$se_variance, digits), check.names = FALSE)
  tab$SE[is.na(m$se_mean)] <- ""
  tab$`SE `[is.na(m$se_variance)] <- ""
  print(tab, row.names = FALSE, right = TRUE)

  if (!is.null(x$mean_comparison)) {
    cat(sprintf("\nFactor mean comparison at the %s%% level in descending order\n",
                format(100 * x$settings$level)))
    mc <- x$mean_comparison
    gw <- max(nchar(c("Group", mc$group)))
    cat(sprintf(" %4s  %-*s  %8s   %s\n", "Rank", gw, "Group", "Mean",
                "Groups with significantly smaller mean"))
    for (r in seq_len(nrow(mc))) {
      lines <- .wrap_labels(mc$smaller[[r]], width - gw - 20, indent = gw + 20)
      cat(sprintf(" %4d  %-*s  %8s   %s\n", mc$rank[r], gw, mc$group[r],
                  .fmt(mc$mean[r], digits), lines[1]))
      if (length(lines) > 1) cat(paste0(lines[-1], "\n"), sep = "")
    }
  }

  cat("\nItem parameters\n")
  s <- x$items
  tab <- data.frame(Item = s$item, Parameter = s$param,
                    Noninvariant = ifelse(is.na(s$noninvariant), "",
                                          paste0(s$noninvariant, "/", s$groups)),
                    `Invariant mean` = ifelse(is.na(s$invariant_mean), "",
                                              .fmt(s$invariant_mean, digits)),
                    R2 = ifelse(is.na(s$r2), "", .fmt(s$r2, digits)),
                    Contribution = .fmt(s$contribution, digits),
                    Variance = ifelse(is.na(s$variance), "", .fmt(s$variance, digits)),
                    check.names = FALSE)
  print(tab, row.names = FALSE, right = TRUE)
  invisible(x)
}

#' @rdname summary.irtalign
#' @export
print.irtalign <- function(x, digits = 3, ...) {
  cat(.settings_lines(x$settings, x$optimization, length(x$groups), length(x$items),
                      x$means$n), sep = "\n")
  cat("\n")
  m <- x$means
  tab <- data.frame(Group = m$group, N = m$n, Mean = .fmt(m$mean, digits),
                    SE = ifelse(is.na(m$se_mean), "", .fmt(m$se_mean, digits)),
                    Variance = .fmt(m$variance, digits),
                    `SE ` = ifelse(is.na(m$se_variance), "", .fmt(m$se_variance, digits)),
                    check.names = FALSE)
  print(tab, row.names = FALSE, right = TRUE)
  invisible(x)
}

#' @rdname summary.irtalign
#' @export
coef.irtalign <- function(object, ...) {
  object$means
}

#' Invariance results by item parameter and group
#'
#' @param object An object of class `irtalign`.
#' @return A data frame with one row per item parameter and group: configural
#'   and aligned estimates in the reported metric, standard errors and whether
#'   the group belongs to the invariant set of the parameter.
#' @export
invariance_table <- function(object) {
  stopifnot(inherits(object, "irtalign"))
  object$params
}

#' Factor mean comparison
#'
#' @param object An object of class `irtalign`.
#' @return A data frame with groups ranked by aligned mean and a list column
#'   `smaller` with the groups whose mean is significantly smaller.
#' @export
mean_comparison <- function(object) {
  stopifnot(inherits(object, "irtalign"))
  object$mean_comparison
}
