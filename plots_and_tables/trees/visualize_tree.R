
plot_features_on_tree_with_complasso <- function(tree, 
                                                 gamma, # from trac 
                                                 A,
                                                 beta_complasso,
                                                 thickness = 0.5,
                                                 scaling = 1) {
  # get the aggregation nodes:
  agg_nodes <- which(abs(gamma) > 1e-8)
  agg_names <- rownames(gamma)[agg_nodes]
  agg_short_names <- str_split(agg_names, "::") %>% map_chr(~ rev(.x)[1])
  # partition the leaves based on this:
  leaf_assign <- apply(A[, agg_nodes], 1, function(row) which(row == 1)[1])
  leaf_assign <- agg_short_names[leaf_assign]
  #agg_short_names <- c(agg_short_names, "Other")
  #leaf_assign[is.na(leaf_assign)] <- "Other"
  # group the leaves so that the tree will be colored accordingly:
  grp <- agg_short_names %>% 
    map(~ which(leaf_assign == .x)) %>% 
    set_names(agg_short_names)
  tree <- groupOTU(tree, grp, "Selected Taxa")
  plt <- ggtree(tree, aes(color = `Selected Taxa`,
                          alpha = 1 - 0.01 * (`Selected Taxa` == "Other" |
                                               str_detect(`Selected Taxa`, "Bacteria"))
                          ),
                layout = "circular",
                branch.length = "none", 
                size = thickness
                ) + 
    theme(legend.position="bottom") +
    scale_alpha(guide = 'none')
  if ("0" %in% levels(plt$data$`Selected Taxa`))
    levels(plt$data$`Selected Taxa`)[levels(plt$data$`Selected Taxa`) == "0"] <- "Other"
  # put a heatmap along the perimeter giving beta's value:
  beta <- data.frame(beta = as.numeric(A %*% gamma))
  beta[beta == 0] <- NA
  rownames(beta) <- tree$tip.label
  plt <- gheatmap_modified(plt, beta, width = 0.1, color = NA, 
                           legend_title = "beta", scaling = scaling)

  plt$data <- plt$data %>%
    mutate(otu = str_extract(label,"::[^:]+$") %>% 
             str_remove("::")) %>% # extract leaf-level part of taxon label
    left_join(tibble(beta_cl = beta_complasso[,1],
                     otu = rownames(beta_complasso)),
              by = "otu")
  plt <- plt + geom_tippoint(data = plt$data %>% filter(beta_cl != 0),
                             color = "black", size = 0.5, alpha = 1)
  return(plt)
}
