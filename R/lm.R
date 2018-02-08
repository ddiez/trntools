#' Fits a linear model to all genes using a matrix of regulators' activities as predictors.
#'
#' @param x matrix of regulators' activities.
#' @param y matrix of genes' expression levels.
#'
#' @export
trm_lm <- function(x, y) {
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

  list(coef = coefs, pval = pvals)
}
