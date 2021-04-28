# https://bioconductor.org/install/
# BiocManager::install("ggtree")
library(Matrix)
library(tidyverse)
library(ggtree)
source("gheatmap_modified.R")
source("visualize_tree.R")

sCD14_trac_figure <- function() {
  dat <- readRDS("../../sCD14_HIV/sCD14.RDS")
  load("../../sCD14_HIV/trac_fixed_level_multi_splits.Rdata")
  i1se <- cvfit[[1]]$OTU$cv[[1]]$i1se
  
  load("../../sCD14_HIV/complasso_fixed_level_multi_splits.Rdata",
       envir = cl <- new.env())
  i1se_cl <- cl$cvfit[[1]]$OTU$cv$i1se
  gam <- fit[[1]]$OTU[[1]]$gamma[, i1se, drop = FALSE]  
  prefix <- c("", "k__", "p__", "c__", "o__", "f__", "g__", "")
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
    ggsave("sCD14_HIV_trac.pdf",
             width = 6.5,
             height = 4.58)
}

sCD14_trac_weights_figure <- function() {
  dat <- readRDS("../../sCD14_HIV/sCD14.RDS")
  load("../../sCD14_HIV/trac_weights_fixed_level_multi_splits.Rdata")
  i1se <- cvfit[[1]]$OTU$cv[[1]]$i1se
  
  load("../../sCD14_HIV/complasso_fixed_level_multi_splits.Rdata",
       envir = cl <- new.env())
  i1se_cl <- cl$cvfit[[1]]$OTU$cv$i1se

  gam <- fit[[1]]$OTU[[1]]$gamma[, i1se, drop = FALSE]  
  prefix <- c("", "k__", "p__", "c__", "o__", "f__", "g__", "")
  rownames(gam) <- rownames(gam) %>% 
    str_split("::") %>% 
    map(~ paste0(prefix[1:length(.x)], .x, collapse = "::"))
  g <- plot_features_on_tree_with_complasso(dat$tree, gam,
                                            A = dat$A,
                                            beta_complasso = cl$fit[[1]]$OTU$beta[, i1se_cl, drop = FALSE])
  # make lines thicker in "Selected Taxa" legend
  g <- g + guides(color = guide_legend(override.aes = list(size = 2),
                                       label.theme = element_text(size = 8),
                                       keyheight = unit(0.1, "lines")
                                       )) 
  g + 
    labs(title = "B: trac (a = 1/2)") +
  ggsave("sCD14_HIV_trac_weights.pdf",
             width = 6.5,
             height = 4.58)
}

sCD14_trac_family_on_otu_tree_figure <- function() {
  dat_list <- readRDS("../../sCD14_HIV/sCD14_aggregated.RDS")
  load("../../sCD14_HIV/trac_fixed_level_multi_splits.Rdata")
  i1se <- cvfit[[1]]$Family$cv[[1]]$i1se
  
  load("../../sCD14_HIV/complasso_fixed_level_multi_splits.Rdata",
       envir = cl <- new.env())
  i1se_cl <- cl$cvfit[[1]]$Family$cv$i1se
  
  gamma <- fit[[1]]$Family[[1]]$gamma[, i1se, drop = FALSE]
  
  names_gamma_template <- rownames(fit[[1]]$OTU[[1]]$gamma)
  gamma_template <- matrix(0, nrow = length(names_gamma_template), ncol = 1)
  rownames(gamma_template) <- names_gamma_template
  
  ii_gamma_template <- match(rownames(gamma),rownames(gamma_template))
  ii_gamma <- which(!is.na(ii_gamma_template))
  gamma_template[ii_gamma_template[!is.na(ii_gamma_template)], ] <- gamma[ii_gamma]
  
  for (j in which(is.na(ii_gamma_template))) {
    # loop over families with single child (and therefore a node corresponding to
    # this family is not present in OTU-level tree)
    family_name <- rownames(gamma)[j]
    # get all descendants of this family:
    ind <- str_which(rownames(gamma_template), family_name)
    if (length(ind) > 1) {
      # if there are multiple descendants, let's find the child of the family
      # (which must be unique) so that we can assign gamma[j] to that taxon's 
      # gamma_template value:
      ii <- ind[which.min(str_count(names_gamma_template[ind], "::"))]
    } else {
      # if a unique descendant exists, then we will just assign gamma[j] to that descendant:
      ii <- ind
    }
    # do the assignment:
    gamma_template[ii] <- gamma[j]
    # let's remove the last part of the name so that it will be labeled
    # as a family rather than a genus:
    rownames(gamma_template)[ii] <- str_remove(rownames(gamma_template)[ii],
                                               "::[^:]+$")
  }
  
  otu_names <- rownames(dat_list$OTU$A)
  beta_cl <- matrix(0, nrow = length(otu_names), ncol = 1)
  rownames(beta_cl) <- otu_names
  
  beta_family_cl <- cl$fit[[1]]$Family$beta[, i1se_cl, drop = FALSE]
  cl_selected_families <- names(beta_family_cl[beta_family_cl!=0,])
  for (nam in cl_selected_families) {
    ii <- str_detect(otu_names, nam)
    beta_cl[ii,] <- 1
  }
  rownames(beta_cl) <- rownames(beta_cl) %>% 
    str_extract("::[^:]+$") %>% 
    str_remove("::")
  
  prefix <- c("", "k__", "p__", "c__", "o__", "f__", "g__", "")
  rownames(gamma_template) <- rownames(gamma_template) %>% 
    str_split("::") %>% 
    map(~ paste0(prefix[1:length(.x)], .x, collapse = "::"))

  g <- plot_features_on_tree_with_complasso(dat_list$OTU$tree, 
                                            gamma_template,
                                            A = dat_list$OTU$A,
                                            beta_complasso = beta_cl)
  # make lines thicker in "Selected Taxa" legend
  g <- g + guides(color = guide_legend(override.aes = list(size = 2)))
  g + 
    labs(title = "C: trac on family level (a = 1)") +
    ggsave("sCD14_HIV_trac_family_on_otu_tree.pdf",
           width = 6.5,
           height = 4.58)
}

sCD14_trac_weights_family_on_otu_tree_figure <- function() {
  dat_list <- readRDS("../../sCD14_HIV/sCD14_aggregated.RDS")
  load("../../sCD14_HIV/trac_weights_fixed_level_multi_splits.Rdata")
  i1se <- cvfit[[1]]$Family$cv[[1]]$i1se
  
  load("../../sCD14_HIV/complasso_fixed_level_multi_splits.Rdata",
       envir = cl <- new.env())
  i1se_cl <- cl$cvfit[[1]]$Family$cv$i1se
  
  gamma <- fit[[1]]$Family[[1]]$gamma[, i1se, drop = FALSE]
  
  names_gamma_template <- rownames(fit[[1]]$OTU[[1]]$gamma)
  gamma_template <- matrix(0, nrow = length(names_gamma_template), ncol = 1)
  rownames(gamma_template) <- names_gamma_template
  
  ii_gamma_template <- match(rownames(gamma),rownames(gamma_template))
  ii_gamma <- which(!is.na(ii_gamma_template))
  gamma_template[ii_gamma_template[!is.na(ii_gamma_template)], ] <- gamma[ii_gamma]
  
  for (j in which(is.na(ii_gamma_template))) {
    # loop over families with single child (and therefore a node corresponding to
    # this family is not present in OTU-level tree)
    family_name <- rownames(gamma)[j]
    # get all descendants of this family:
    ind <- str_which(rownames(gamma_template), family_name)
    if (length(ind) > 1) {
      # if there are multiple descendants, let's find the child of the family
      # (which must be unique) so that we can assign gamma[j] to that taxon's 
      # gamma_template value:
      ii <- ind[which.min(str_count(names_gamma_template[ind], "::"))]
    } else {
      # if a unique descendant exists, then we will just assign gamma[j] to that descendant:
      ii <- ind
    }
    # do the assignment:
    gamma_template[ii] <- gamma[j]
    # let's remove the last part of the name so that it will be labeled
    # as a family rather than a genus:
    rownames(gamma_template)[ii] <- str_remove(rownames(gamma_template)[ii],
                                               "::[^:]+$")
  }

  otu_names <- rownames(dat_list$OTU$A)
  beta_cl <- matrix(0, nrow = length(otu_names), ncol = 1)
  rownames(beta_cl) <- otu_names
  
  beta_family_cl <- cl$fit[[1]]$Family$beta[, i1se_cl, drop = FALSE]
  cl_selected_families <- names(beta_family_cl[beta_family_cl!=0,])
  for (nam in cl_selected_families) {
    ii <- str_detect(otu_names, nam)
    beta_cl[ii,] <- 1
  }
  rownames(beta_cl) <- rownames(beta_cl) %>% 
    str_extract("::[^:]+$") %>% 
    str_remove("::")

  prefix <- c("", "k__", "p__", "c__", "o__", "f__", "g__", "")
  rownames(gamma_template) <- rownames(gamma_template) %>% 
    str_split("::") %>% 
    map(~ paste0(prefix[1:length(.x)], .x, collapse = "::"))
  
  
  g <- plot_features_on_tree_with_complasso(dat_list$OTU$tree, 
                                            gamma_template,
                                            A = dat_list$OTU$A,
                                            beta_complasso = beta_cl)
  # make lines thicker in "Selected Taxa" legend
  g <- g + guides(color = guide_legend(override.aes = list(size = 2)))
  g + 
    labs(title = "D: trac on family level (a = 1/2)") +
    ggsave("sCD14_HIV_trac_weights_family_on_otu_tree.pdf",
           width = 6.5,
           height = 4.58)
}

