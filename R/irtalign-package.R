#' irtalign: alignment optimization for multiple-group item response models
#'
#' Estimates latent means and variances across many groups without exact
#' measurement invariance (Asparouhov & Muthen, 2014). Dichotomous and
#' polytomous items can be mixed.
#'
#' Main functions:
#' * [align_irt()]: fits each group with [mirt::mirt()] and aligns the
#'   estimates;
#' * [align_fits()], [align_pars()]: align fitted models or parameter
#'   estimates;
#' * [summary.irtalign()], [invariance_table()], [mean_comparison()]: results.
#'
#' @keywords internal
#' @useDynLib irtalign, .registration = TRUE
#' @importFrom Rcpp sourceCpp
"_PACKAGE"
