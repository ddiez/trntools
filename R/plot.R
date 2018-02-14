trn2tbl_graph <- function(x) {
  x <- graph_from_incidence_matrix(x, directed = TRUE, mode = "out", weighted = TRUE)

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

#' Plot transcriptional regulatory networks.
#'
#' @param x a network object with attributes.
#'
#' @export
#' @rdname plot_trn
plot_trn <- function(x, ...) {
  UseMethod("plot_trn")
}

#' @rdname plot_trn
#' @export
plot_trn.data.frame <- function(x, node.size = 10, label.size = 3, regulator.color = "palevioletred1", target.color = "turquoise3") {
  g <- ggraph(x) + theme_graph(base_family = "Arial", base_size = 12)

  label <- function(x) {
    l <- rep("repressor", length(x))
    l[x == "red"] <- "activator"
    l
  }


  g <- g + geom_node_point(
    aes(shape = regulator, fill = regulator, color = regulator),
    size = node.size) +
    geom_node_text(aes(label = name), size = label.size) +
    scale_fill_manual(name = "regulator", values = c(target.color, regulator.color)) +
    scale_color_manual(name = "regulator", values = c(target.color, regulator.color)) +
    scale_shape_manual(name = "regulator", values = c(21, 22)) +
    guides(fill = guide_legend(order = 1), color = guide_legend(order = 1), shape = guide_legend(order = 1))

  ne <- attr(x, "graph") %>% activate(edges) %>% as_tibble() %>% nrow()

  if (ne > 0) {
    g  <- g +
      geom_edge_fan(
        aes(width = abs(activity), colour = c("steelblue1", "palevioletred1")[role]),
        start_cap = circle(10, 'points'),
        end_cap = circle(10, 'points'),
        arrow = arrow(type = "closed", length = unit(10, 'points'))) +
      scale_edge_color_identity("role", labels = label, guide = "legend") +
      scale_edge_width("strength", range = c(1, 2), guide = "legend")
  }
  g
}

#' @rdname plot_trn
#' @export
plot_trn.tbl_graph <- function(x, layout = "nicely", ...) {
  plot_trn(create_layout(x, layout = layout), ...)
}
