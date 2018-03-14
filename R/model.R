trn_empty <- function(N = 10, nr = 3, group = NULL) {
  ng <- N - nr
  gn <- paste0("G", seq_len(ng))
  rn <- paste0("R", seq_len(nr))

  create_empty(N) %>% activate(nodes) %>%
    mutate(name = c(rn, gn), regulator = rep(c(TRUE, FALSE), times = c(nr, ng))) %>%
    activate(edges) %>%
    mutate(from = numeric(), to = numeric(), activity = numeric(), role = character(), group = character())
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

to_role <- function(x) {
  c("activator", "repressor")[as.numeric(x > 0) + 1]
}

add_edge <- function(x, from, to, activity, group = NULL) {
  nodes <- to_nodes(x)
  edges <- to_edges(x)

  if (is.null(group)) group <- "group"

  new <- tibble(from = from, to = to, activity = activity, role = to_role(activity), group = group)

  edges <- edges %>% bind_rows(new)

  tbl_graph(nodes, edges, directed = TRUE)
}

to_matrix <- function(x) {
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


foo <- trn_empty()
foo <- foo %>% add_edge(from = c(1, 2), to = c(5, 10), activity = c(1, -1), group = "A")
foo <- foo %>% add_edge(from = c(1, 3), to = c(5, 10), activity = c(1, -1), group = "B")

foo <- foo %>% activate(edges) %>% mutate(group = factor(group, levels = c("A", "B", "C")))

plot_network(foo)
plot_network(foo) + facet_edges(~group, drop = FALSE, ncol = 2)
