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
plot_trn.data.frame <- function(x, node.size = 10, label.size = 3, regulator.color = "palevioletred1", target.color = "turquoise3", color.by = NULL, color.low = "blue", color.mid = "white", color.high = "red") {
  g <- ggraph(x) + theme_graph()

  label <- function(x) {
    l <- rep("repressor", length(x))
    l[x == "palevioletred1"] <- "activator"
    l
  }

  ne <- attr(x, "graph") %>% activate(edges) %>% as_tibble() %>% nrow()

  if (ne > 0) {
    g  <- g +
      geom_edge_fan(
        aes(width = abs(activity), colour = c("steelblue1", "palevioletred1")[role]),
        start_cap = circle(node.size + 2, 'points'),
        end_cap = circle(node.size + 2, 'points'),
        arrow = arrow(type = "closed", length = unit(5, 'points'))) +
      scale_edge_color_identity("role", labels = label, guide = "legend") +
      scale_edge_width("strength", range = c(.2, 1), guide = "legend")
  }

  if (is.null(color.by)) {
    g <- g +
      geom_node_point(
        aes(shape = regulator, fill = regulator, color = regulator),
        size = node.size) +
      geom_node_text(aes(label = name), size = label.size) +
      scale_fill_manual(name = "regulator", values = c(target.color, regulator.color)) +
      scale_color_manual(name = "regulator", values = c(target.color, regulator.color)) +
      scale_shape_manual(name = "regulator", values = c(21, 22)) +
      guides(fill = guide_legend(order = 1), color = guide_legend(order = 1), shape = guide_legend(order = 1))
  } else {
    g <- g +
      geom_node_point(
        aes(shape = regulator, fill = .data[[color.by]], color = .data[[color.by]]),
        size = node.size) +
      geom_node_text(aes(label = name), size = label.size) +
      scale_color_gradient2(name = "expression", low = color.low, mid = color.mid, high = color.high) +
      scale_fill_gradient2(name = "expression", low = color.low, mid = color.mid, high = color.high) +
      scale_shape_manual(name = "regulator", values = c(21, 22)) +
      guides(shape = guide_legend(order = 1))
  }
  g
}

#' @rdname plot_trn
#' @export
plot_trn.tbl_graph <- function(x, layout = "nicely", ...) {
  plot_trn(create_layout(x, layout = layout), ...)
}


#' Plot a heatmap of a trn.
#'
#' @param x a TRN object.
#' @param guide logical; whether to plot the guide.
#' @param limits set limits to scale_fill_gradient2().
#'
#' @export
plot_heatmap <- function(x, guide = TRUE, limits = NULL) {
  UseMethod("plot_heatmap")
}

#' @rdname plot_heatmap
#' @export
plot_heatmap.matrix <- function(x, guide = TRUE, limits = NULL) {
  d <- x %>% as.data.frame() %>% rownames_to_column("regulator") %>%
    gather("gene", "activity", -regulator) %>%
    mutate(regulator = factor(regulator, rownames(x))) %>%
    mutate(gene = factor(gene, levels = colnames(x)))

  plot_heatmap(d, guide = guide, limits = limits)
}

#' @rdname plot_heatmap
#' @export
plot_heatmap.data.frame <- function(x, guide = TRUE, limits = NULL) {
  p <- ggplot(x, aes(x = gene, y = regulator, fill = activity)) + geom_tile() +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_gradient2(low = "steelblue1", mid = "white", high = "palevioletred1", midpoint = 0, guide = guide_legend(reverse = TRUE), limits = limits) +
    labs(x = NULL, y = NULL, title = NULL)

  if (!guide)
    p <- p + guides(fill = FALSE)

  if (ncol(x) > 20)
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

  if (nrow(x) > 20)
    p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

  p
}



#' plot_trn_heatmap
#'
#' @param x
#'
#' @export
#'
plot_trn_heatmap <- function(x) {
  mycol <- circlize::colorRamp2(c(-1 * max(abs(range(x))), 0, max(abs(range(x)))), c("blue", "white", "red"))
  ComplexHeatmap::Heatmap(x, show_column_names = FALSE, name = "coef", col = mycol, cluster_rows = FALSE, column_title = "Target genes", row_title = "Transcription factors")
}

#' plot_trn_graph
#'
#' @param x
#'
#' @export
#'
plot_trn_graph <- function(x, cutoff = NULL, layout = "nicely", color.by = NULL, seed = 42, ...) {
  set.seed(seed)

  if (!is.null(cutoff)) {
    x[abs(x) < cutoff] <- 0
    x <- x[rowSums(x) != 0, colSums(x) != 0, drop = FALSE]
  }

  x <- trn2tbl_graph(x)
  l <- create_layout(x, layout)

  plot_trn(l, color.by = color.by, ...)
}
