library(tidygraph)
library(dplyr)
library(tidyr)
library(ggraph)
library(Matrix)

#' +.tbl_graph
#'
#' Combine two TRN models.
#'
#' @param e1 a TRN object.
#' @param e2 a TRN object.
#'
#' @return a TRN object.
#'
#' @export
`+.tbl_graph` <- function(e1, e2) {
  e1 %>% bind_edges(e2 %>% activate("edges") %>% as_tibble())
}


trn_empty <- function(N = 10, nr = 3, group = NULL) {
  ng <- N - nr
  gn <- paste0("G", seq_len(ng))
  rn <- paste0("R", seq_len(nr))

  create_empty(N, directed = TRUE) %>% activate(nodes) %>%
    mutate(name = c(rn, gn), regulator = rep(c(TRUE, FALSE), times = c(nr, ng))) %>%
    activate(edges) %>%
    mutate(from = numeric(), to = numeric(), activity = numeric(), role = character(), group = character())
}

evolve_random <- function(x, prob = c(0.1, 1, 0.1), group = NULL) {
  nodes <- x %>%
    to_nodes() %>%
    rowid_to_column("id") %>%
    select("id", "regulator")

  from_ids <- nodes %>% filter(regulator) %>% pull("id")
  to_ids <- nodes %>% filter(! regulator) %>% pull("id")

  net <- expand.grid(from = from_ids, to = to_ids)
  #net <- net %>% sample_n(n, replace = FALSE)
  #net <- net %>% mutate(activity = sample(c(-1, 1), n, replace = TRUE))

  net <- net %>% left_join(edges, by = c("from", "to")) %>% replace_na(list(activity = 0, role = "neutral", group = "group"))
  net <- net %>% mutate(activity = sample(c(-1, 0, 1), n(), prob = prob, replace = TRUE), role = to_role(activity)) %>% filter(activity != 0)

  x %>% add_edge(from = net[["from"]], to = net[["to"]], activity = net[["activity"]], group = group)
}

to_nodes <- function(x) {
  x %>% activate(nodes) %>% as_tibble()
}

to_edges <- function(x) {
  x %>% activate(edges) %>% as_tibble()
}

get_regulators <- function(x) {
  x <- to_nodes(x) %>% filter(regulator)
  x[["name"]]
}

get_genes <- function(x, regulators = FALSE) {
  x <- to_nodes(x)

  if (! regulators)
    x <- x %>% filter(! regulator)

  x[["name"]]
}

get_node_index <- function(x, name) {
  nodes <- to_nodes(x)
  nodes %>% filter(name == !!name) %>% pull(id)
}

to_role <- function(x) {
  c("repressor", "neutral", "activator")[sign(x) + 2]
}

add_edge <- function(x, from, to, activity, group = NULL) {
  nodes <- to_nodes(x)
  edges <- to_edges(x)

  if (is.null(group)) group <- "group"

  if (is.character(from))
    from <- get_node_index(x, from)
  if (is.character(to))
    to <- get_node_index(x, to)

  new <- tibble(from = from, to = to, activity = activity, role = to_role(activity), group = group)

  edges <- edges %>% bind_rows(new)

  tbl_graph(nodes, edges, directed = TRUE)
}

to_matrix <- function(x, ...) {
  UseMethod("to_matrix")
}
to_matrix.tbl_graph <- function(x, regulators = FALSE) {
  reg <- get_regulators(x)
  gen <- get_genes(x, regulators = TRUE)

  y <- Matrix::Matrix(0, nrow = length(reg), ncol = length(gen), dimnames = list(reg, gen))

  edges <- to_edges(x)

  for (k in seq_len(nrow(edges))) {
    y[edges[[k, "from"]], edges[[k, "to"]]] <- edges[[k, "activity"]]
  }
  if (!regulators)
    y <- y[, get_genes(x)]
  y
}

to_trn <- function(x) {
  UseMethod("to_trn")
}

to_trn.dgCMatrix <- function(x) {
  edges <- summary(x) %>% as_tibble() %>% rename(from = i, to = j, activity = x)
  nodes <- tibble(name = c(rownames(x), colnames(x))) %>% mutate(regulator = name %in% rownames(x))
  tbl_graph(nodes, edges)
}

plot_heatmap <- function(x, guide = TRUE) {
  UseMethod("plot_heatmap")
}

plot_heatmap.tbl_graph <- function(x) {
  plot_heatmap(to_matrix(x))
}

plot_heatmap.matrix <- function(x, guide = TRUE) {
  d <- x %>% as.data.frame() %>% rownames_to_column("regulator") %>%
    gather("gene", "activity", -regulator) %>%
    mutate(regulator = factor(regulator, rownames(x))) %>%
    mutate(gene = factor(gene, levels = colnames(x)))

  plot_heatmap(d)
}

plot_heatmap.data.frame <- function(x, guide = TRUE, axis = c(TRUE, TRUE)) {
  p <- ggplot(x, aes(x = gene, y = regulator, fill = activity)) + geom_tile() +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_gradient2(low = "steelblue1", mid = "white", high = "palevioletred1", midpoint = 0, guide = guide_legend(reverse = TRUE) ) +
    labs(x = NULL, y = NULL, title = NULL)

  if (!guide)
    p <- p + guides(fill = FALSE)

  if (axis[1])
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

  if (axis[2])
    p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

  p
}


plot_network <- function(x) {

  p <- ggraph(x, "circle") + theme_graph(base_family = "Arial", strip_text_size = 20)

  edges <- to_edges(x)
  ne <- nrow(edges)
  if (ne > 0) {
    p <- p +
      geom_edge_fan(
        aes(color = role, width = abs(activity)),
        arrow = arrow(length = unit(5, "mm")),
        start_cap = circle(4, "mm"),
        end_cap = circle(4, "mm")) +
      scale_edge_width("strength", range = c(1, 3), breaks = seq(1, 3, 1)) +
      scale_edge_color_manual(values = c(activator = "violetred", repressor = "slateblue3")) +
      facet_edges(~group, drop = FALSE)
  }
  p + geom_node_text(aes(label = name, color = regulator)) +
    scale_color_manual(values = c("TRUE" = "black", "FALSE" = "limegreen"))
}


# foo <- trn_empty()
# foo <- foo %>% add_edge(from = c(1, 2), to = c(5, 10), activity = c(1, -1), group = "A")
# foo <- foo %>% add_edge(from = c(1, 3), to = c(5, 10), activity = c(1, -1), group = "B")
#
# foo <- foo %>% activate(edges) %>% mutate(group = factor(group, levels = c("A", "B", "C")))
#
# plot_network(foo)
# plot_network(foo) + facet_edges(~group, drop = FALSE, ncol = 2) + theme(legend.position = "none")
#
#
# foo %>% to_matrix()
# foo2 <- foo %>% to_matrix() %>% to_trn()
# identical(foo, foo2)
