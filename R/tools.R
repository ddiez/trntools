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
#' @param direction whether to use all edges, only positive one or negative ones.
#'
#' @export
trn_compare <- function(x, y, direction = "all") {
  match.arg(direction, c("all", "positive", "negative"))

  switch(direction,
         all = check_direction(x, y),
         positive = check_positive(x, y),
         negative = check_negative(x, y))
}

check_direction <- function(x, y) {
  n <- length(x)
  ndiff <- sum((sign(x) - sign(y)) != 0)
  (n - ndiff) / n
}

check_positive <- function(x, y) {
  sel <- x > 0
  xs <- sum(sel)
  ys <- sum(y[sel] > 0)
  ys / xs
}

check_negative <- function(x, y) {
  sel <- x < 0
  xs <- sum(sel)
  ys <- sum(y[sel] < 0)
  ys / xs
}
