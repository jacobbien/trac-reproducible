# https://bioconductor.org/install/
# BiocManager::install("ggtree")
library(Matrix)
library(tidyverse)
library(ggtree)
source("gheatmap_modified.R")
source("visualize_tree.R")

cps_pH_trac_figure <- function() {
  dat <- readRDS("../../CentralParkSoil/cps_pH_processed.RDS")
  load("../../CentralParkSoil/cps_pH_trac_fixed_level_multi_splits.Rdata")
  i1se <- cvfit[[1]]$OTU$cv[[1]]$i1se
  
  load("../../CentralParkSoil/cps_pH_complasso_fixed_level_multi_splits.Rdata",
       envir = cl <- new.env())
  i1se_cl <- cl$cvfit[[1]]$OTU$cv$i1se
  gam <- fit[[1]]$OTU[[1]]$gamma[, i1se, drop = FALSE]
  g <- plot_features_on_tree_with_complasso(dat$tree, gam, A = dat$A,
                                            beta_complasso = cl$fit[[1]]$OTU$beta[, i1se_cl, drop = FALSE],
                                            scaling = max(abs(gam))/20)
  # make lines thicker in "Selected Taxa" legend
  g <- g + guides(color = guide_legend(override.aes = list(size = 2),
                                       label.theme = element_text(size = 10)))
  g + 
    labs(title = "A: trac (a = 1)")
  ggsave("cps_pH_trac.pdf",
         width = 6.5,
         height = 4.58)
  g + 
    labs(title = "A: pH")
  ggsave("cps_pH_trac_appendix.pdf",
         width = 6.5,
         height = 4.58)
  
}

cps_pH_trac_weights_figure <- function() {
  dat <- readRDS("../../CentralParkSoil/cps_pH_processed.RDS")
  load("../../CentralParkSoil/cps_pH_trac_weights_fixed_level_multi_splits.Rdata")
  i1se <- cvfit[[1]]$OTU$cv[[1]]$i1se
  
  load("../../CentralParkSoil/cps_pH_complasso_fixed_level_multi_splits.Rdata",
       envir = cl <- new.env())
  i1se_cl <- cl$cvfit[[1]]$OTU$cv$i1se
  
  gam <- fit[[1]]$OTU[[1]]$gamma[, i1se, drop = FALSE]
  g <- plot_features_on_tree_with_complasso(dat$tree, gam, A = dat$A,
                                            beta_complasso = cl$fit[[1]]$OTU$beta[, i1se_cl, drop = FALSE],
                                            scaling = max(abs(gam))/100)
  # make lines thicker in "Selected Taxa" legend
  g <- g + guides(color = guide_legend(override.aes = list(size = 2),
                                       label.theme = element_text(size = 5),
                                       keyheight = unit(0.05, "lines")
  )) 
  g + 
    labs(title = "B: trac (a = 1/2)") +
    ggsave("cps_pH_trac_weights.pdf",
           width = 6.5,
           height = 4.58)
}

cps_Mois_trac_figure <- function() {
  # for appendix: to compare Moisture with pH tree
  
  dat <- readRDS("../../CentralParkSoil/cps_Mois_processed.RDS")
  load("../../CentralParkSoil/cps_Mois_trac_fixed_level_multi_splits.Rdata")
  i1se <- cvfit[[1]]$OTU$cv[[1]]$i1se
  
  load("../../CentralParkSoil/cps_Mois_complasso_fixed_level_multi_splits.Rdata",
       envir = cl <- new.env())
  i1se_cl <- cl$cvfit[[1]]$OTU$cv$i1se
  gam <- fit[[1]]$OTU[[1]]$gamma[, i1se, drop = FALSE]
  g <- plot_features_on_tree_with_complasso(dat$tree, gam, A = dat$A,
                                            beta_complasso = cl$fit[[1]]$OTU$beta[, i1se_cl, drop = FALSE],
                                            scaling = max(abs(gam))/20)
  # make lines thicker in "Selected Taxa" legend
  g <- g + guides(color = guide_legend(override.aes = list(size = 2),
                                       label.theme = element_text(size = 4),
                                       keyheight = unit(0.07, "lines"),
                                       keywidth = unit(0.25, "cm"),
                                       title.theme = element_text(size = 8)))
  
  
  g + 
    labs(title = "B: Moisture") +
  ggsave("cps_Mois_trac_appendix.pdf",
         width = 6.5,
         height = 4.58)
}