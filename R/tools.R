#' Convert incidence matrix to tbl_graph.
#'
#' @param x incidence matrix
#'
#' @export
trn2tbl_graph <- function(x) {
  x <- igraph::graph_from_incidence_matrix(x, directed = TRUE, mode = "out", weighted = TRUE)

  x <- x %>% as_tbl_graph %>%
    activate(nodes) %>%
    mutate(regulator = !type, type = !type)

  if (x %>% activate(edges) %>% as_tibble() %>% nrow()) {
    x <- x %>%
      activate(edges) %>%
      mutate(activity = abs(weight), role = factor(
        c("repressor", "activator")[(weight > 0) + 1],
        levels = c("repressor", "activator")
      )) %>%
      select(-weight)
  }

  x
}


#' Computes RMSD between two TRN models.
#'
#' Takes a real and predicted model and computes the root mean squared distance.
#'
#' @param x a TRN object.
#' @param y a TRN object.
#'
#' @export
trn_rmsd <- function(x, y) {
  mean(sqrt((x - y)^2))
}


#' Computes percentage of coherent directions between two TRN models.
#'
#' @param x a TRN object (truth).
#' @param y a TRN object (predicted).
#' @param direction whether to use all edges, empty edges, positive edges or negative edges
#'
#' @export
trn_compare <- function(x, y, direction = "all") {
  match.arg(direction, c("all", "empty", "positive", "negative"))

  switch(direction,
         all = check_direction(x, y),
         empty = check_empty(x, y),
         positive = check_positive(x, y),
         negative = check_negative(x, y))
}

check_direction <- function(x, y) {
  n <- length(x)
  ndiff <- sum((sign(x) - sign(y)) != 0, na.rm = TRUE)
  (n - ndiff) / n
}

check_empty <- function(x, y) {
  sel <- x == 0
  xs <- sum(sel, na.rm = TRUE)
  ys <- sum(y[sel] == 0, na.rm = TRUE)
  ys / xs
}

check_positive <- function(x, y) {
  sel <- x > 0
  xs <- sum(sel, na.rm = TRUE)
  ys <- sum(y[sel] > 0, na.rm = TRUE)
  ys / xs
}

check_negative <- function(x, y) {
  sel <- x < 0
  xs <- sum(sel, na.rm = TRUE)
  ys <- sum(y[sel] < 0, na.rm = TRUE)
  ys / xs
}


#' Filter a TRN at a given FDR.
#'
#' @param x a TRN model.
#' @param fdr FDR level.
#' @param method method for p.adjust().
#'
#' @export
trn_filter <- function(x, fdr = 0.01, method = "bonferroni") {
  xc <- x[["coef"]][-1, , drop = FALSE]
  xp <- x[["pval"]][-1, , drop = FALSE]
  xp <- p.adjust(xp, method = method)
  xc[ xp > fdr] <- 0
  xc
}
