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