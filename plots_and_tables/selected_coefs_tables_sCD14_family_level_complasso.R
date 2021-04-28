library(tidyverse)
library(Matrix)

method_name <- "complasso"
method_long_name <- "the sparse log-contrast method on family level"

data_name <- "Gut (HIV): sCD14"
short_data_name <- "sCD14"

load("../sCD14_HIV/complasso_fixed_level_multi_splits.Rdata",
     envir = cl <- new.env())

dat_list <- readRDS("../sCD14_HIV/sCD14_aggregated.RDS")

cvfit <- cl$cvfit[[1]]$Family
cvfit$cv <- list(cvfit$cv)
out <- list(fit = list(cl$fit[[1]]$Family), 
            cvfit = cvfit, 
            A = dat_list$Family$A, 
            x = dat_list$Family$x)


tax_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "OTU")

ntop <- 15

family_names <- tibble(full = rownames(out$A),
                    short = rownames(out$A) %>% str_extract("[[:alnum:]_]+$"))
beta_cv1se <- round(out$fit[[1]]$beta[, out$cvfit$cv[[1]]$i1se], 6)
stopifnot(names(beta_cv1se) == family_names$full) # make sure order matches
names(beta_cv1se) <- family_names$full %>% str_remove_all("[a-z]__")
into <- tax_levels
selected <- enframe(beta_cv1se) %>%
  filter(value != 0) %>% 
  separate(name,
           into = into,
           sep = "::",
           remove = FALSE,
           fill = "right") %>% 
  mutate(across(1:9, ~ str_replace(.x, "^\\d+$", "[Unclassified]")))
caption <- ifelse(nrow(selected) > ntop,
                  sprintf("Top %s coefficients selected by %s for %s",
                          ntop, method_long_name, data_name),
                  sprintf("Coefficients selected by %s for %s", 
                          method_long_name, data_name))
# add latex label:
caption <- sprintf("\\label{fig:%s-family-level-%s} %s", 
                   short_data_name,
                   "complasso",
                   caption)
selected %>%
  #      arrange(!!!rlang::syms(tax_levels)) %>%
  arrange(-abs(value)) %>% 
  select(-name) %>%
  mutate(value = formatC(value, digits = 2, format = "f")) %>% 
  mutate_all(~ replace_na(., "")) %>% 
  mutate_at(vars(-OTU), ~ str_replace(., "^\\d+$", "-")) %>% 
  slice(1:ntop) %>% # top 15 only
  knitr::kable(format = "latex",
               align = c(rep("c", ncol(selected) - 2), "r"),
               booktabs = TRUE, 
               caption = caption
  ) %>% 
  #kableExtra::kable_styling() %>% 
  kableExtra::kable_styling(latex_options = c("striped", "scale_down")) %>% 
  kableExtra::row_spec(0, bold = TRUE) %>%
  # str_sub(29,-13) %>% 
  str_replace("value", "$\\\\beta$") %>% 
  cat(file = paste0("selected_tables/", short_data_name, "_family-level_", method_name, ".tex"))


