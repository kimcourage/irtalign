# Internal helpers -----------------------------------------------------------

`%||%` <- function(x, y) if (is.null(x)) y else x

# Run expr with a fixed seed and restore the caller's RNG state afterwards.
.with_seed <- function(seed, expr) {
  if (is.null(seed)) return(expr)
  has_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (has_seed) old <- get(".Random.seed", envir = globalenv())
  on.exit({
    if (has_seed) {
      assign(".Random.seed", old, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  })
  set.seed(seed)
  expr
}

# Pair weights w_gh = w_g * w_h with w_g = sqrt(N_g) / mean(sqrt(N)).
.pair_weights <- function(n) {
  w <- sqrt(n) / mean(sqrt(n))
  outer(w, w)
}

.simplicity_power <- function(simplicity) {
  switch(simplicity, sqrt = 0.5, fourthrt = 0.25)
}

# Wrap a vector of labels into lines no wider than `width`.
.wrap_labels <- function(x, width = 70, indent = 0) {
  if (!length(x)) return("")
  lines <- character()
  cur <- ""
  for (lab in x) {
    cand <- if (nzchar(cur)) paste(cur, lab) else lab
    if (nchar(cand) > width && nzchar(cur)) {
      lines <- c(lines, cur)
      cur <- lab
    } else {
      cur <- cand
    }
  }
  lines <- c(lines, cur)
  pad <- strrep(" ", indent)
  paste0(c("", rep(pad, length(lines) - 1)), lines)
}

.resolve_group <- function(ref, groups, arg = "ref") {
  if (is.numeric(ref)) {
    if (length(ref) != 1 || ref < 1 || ref > length(groups) || ref != round(ref))
      stop("`", arg, "` must be a group index between 1 and ", length(groups), ".",
           call. = FALSE)
    return(as.integer(ref))
  }
  k <- match(as.character(ref), groups)
  if (is.na(k)) stop("`", arg, "` does not match a group: ", ref, call. = FALSE)
  k
}

.param_label <- function(item, param) {
  paste(item, param)
}
