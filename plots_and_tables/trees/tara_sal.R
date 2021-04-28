# https://bioconductor.org/install/
# BiocManager::install("ggtree")
library(Matrix)
library(tidyverse)
library(ggtree)
source("gheatmap_modified.R")
source("visualize_tree.R")

tara_sal_trac_figure <- function() {
  dat <- readRDS("../../Tara/tara_sal_processed.RDS")
  load("../../Tara/tara_sal_trac_fixed_level_multi_splits.Rdata")
  i1se <- cvfit[[1]]$OTU$cv[[1]]$i1se
  
  load("../../Tara/tara_sal_complasso_fixed_level_multi_splits.Rdata",
       envir = cl <- new.env())
  i1se_cl <- cl$cvfit[[1]]$OTU$cv$i1se
  
  gam <- fit[[1]]$OTU[[1]]$gamma[, i1se, drop = FALSE]  
  prefix <- c("", "k__", "p__", "c__", "o__", "f__", "g__", "s__", "")
  rownames(gam) <- rownames(gam) %>% 
    str_split("::") %>% 
    map(~ paste0(prefix[1:length(.x)], .x, collapse = "::"))
  
  g <- plot_features_on_tree_with_complasso(dat$tree, gam,
                                            A = dat$A,
                                            beta_complasso = cl$fit[[1]]$OTU$beta[, i1se_cl, drop = FALSE])
  # make lines thicker in "Selected Taxa" legend
  g <- g + guides(color = guide_legend(override.aes = list(size = 2),
                                       label.theme = element_text(size = 10)))
  g + 
    labs(title = "A: trac (a = 1)")
  ggsave("tara_sal_trac.pdf",
         width = 6.5,
         height = 4.58)
}

tara_sal_trac_weights_figure <- function() {
  dat <- readRDS("../../Tara/tara_sal_processed.RDS")
  load("../../Tara/tara_sal_trac_weights_fixed_level_multi_splits.Rdata")
  i1se <- cvfit[[1]]$OTU$cv[[1]]$i1se
  
  load("../../Tara/tara_sal_complasso_fixed_level_multi_splits.Rdata",
       envir = cl <- new.env())
  i1se_cl <- cl$cvfit[[1]]$OTU$cv$i1se
  
  gam <- fit[[1]]$OTU[[1]]$gamma[, i1se, drop = FALSE]  
  prefix <- c("", "k__", "p__", "c__", "o__", "f__", "g__", "s__", "")
  rownames(gam) <- rownames(gam) %>% 
    str_split("::") %>% 
    map(~ paste0(prefix[1:length(.x)], .x, collapse = "::"))
  
  g <- plot_features_on_tree_with_complasso(dat$tree, gam,
                                            A = dat$A,
                                            beta_complasso = cl$fit[[1]]$OTU$beta[, i1se_cl, drop = FALSE],
                                            scaling = max(abs(gam))/100)
  # make lines thicker in "Selected Taxa" legend
  g <- g + guides(color = guide_legend(override.aes = list(size = 2),
                                       label.theme = element_text(size = 8),
                                       keyheight = unit(0.05, "lines"))
  ) 
  g + 
    labs(title = "B: trac (a = 1/2)") +
    ggsave("tara_sal_trac_weights.pdf",
           width = 6.5,
           height = 4.58)
}