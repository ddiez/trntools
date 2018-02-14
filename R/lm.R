#' Fits a linear model to all genes using a matrix of regulators' activities as predictors.
#'
#' @param x matrix of regulators' activities.
#' @param y matrix of genes' expression levels.
#'
#' @export
trn_lm <- function(x, y) {
  n <- ncol(x)
  mod <- lm(y ~ x)

  sum <- summary(mod)

  coefs <- vapply(sum, function(m) {
    coef(m)[, 1]
  }, rep(NA_real_, n + 1))

  pvals <- vapply(sum, function(m) {
    coef(m)[, 4]
  }, rep(NA_real_, n + 1))


  rownames(coefs) <- c("bias", colnames(x))
  rownames(pvals) <- c("bias", colnames(x))

  colnames(coefs) <- colnames(y)
  colnames(pvals) <- colnames(y)

  z <- list(coef = coefs, pval = pvals)
  class(z) <- "trn"
  z
}

#' Fits a linear model to all genes using a matrix of regulators' activities as predictors
#' and a grouping variable.
#'
#' @param x matrix of regulators' activities.
#' @param y matrix of genes' expression levels.
#' @param group grouping variable.
#'
#' @export
trn_lm_by_group <- function(x, y, group) {
  if (length(group) != nrow(x)) stop("length of group and nrow of x don't agree.")
  if (length(group) != nrow(y)) stop("length of group and nrow of y don't agree.")

  groups <- unique(group)
  lapply(groups, function(g) {
    sel.group <- group == g
    xg <- x[sel.group, , drop = FALSE]
    yg <- y[sel.group, , drop = FALSE]
    trn_lm(xg, yg)
  })
}
