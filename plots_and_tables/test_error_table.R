library(tidyverse)

lower_levels <- c("OTU", "Genus", "Family")
methods_data <- list()

# load sCD14:
load("../sCD14_HIV/complasso_fixed_level_multi_splits.Rdata", e = cl <- new.env())
load("../sCD14_HIV/trac_fixed_level_multi_splits.Rdata", e = trac <- new.env())
load("../sCD14_HIV/trac_weights_fixed_level_multi_splits.Rdata", e = trac_weights <- new.env())
cl$cvfit <- cl$cvfit %>% map(~ .x %>% 
                               map(~ {
                                 .x$cv <- list(.x$cv)
                                 return(.x)
                               }))
methods_data[["Gut (HIV): sCD14"]] <- list(
  "Sparse Log-Contrast" = cl,
  "trac (a = 1)" = trac, 
  "trac (a = 1/2)" = trac_weights
)

# load Tara:
load("../Tara/tara_sal_complasso_fixed_level_multi_splits.Rdata", e = cl <- new.env())
load("../Tara/tara_sal_trac_fixed_level_multi_splits.Rdata", e = trac <- new.env())
load("../Tara/tara_sal_trac_weights_fixed_level_multi_splits.Rdata", e = trac_weights <- new.env())
cl$cvfit <- cl$cvfit %>% map(~ .x %>% 
                               map(~ {
                                 .x$cv <- list(.x$cv)
                                 return(.x)
                               }))
methods_data[["Ocean (TARA): Salinity"]] <- list(
  "Sparse Log-Contrast" = cl,
  "trac (a = 1)" = trac, 
  "trac (a = 1/2)" = trac_weights
)

# load cps_pH:
load("../CentralParkSoil/cps_pH_complasso_fixed_level_multi_splits.Rdata", e = cl <- new.env())
load("../CentralParkSoil/cps_pH_trac_fixed_level_multi_splits.Rdata", e = trac <- new.env())
load("../CentralParkSoil/cps_pH_trac_weights_fixed_level_multi_splits.Rdata", e = trac_weights <- new.env())
cl$cvfit <- cl$cvfit %>% map(~ .x %>% 
                               map(~ {
                                 .x$cv <- list(.x$cv)
                                 return(.x)
                               }))
methods_data[["Central Park Soil: pH"]] <- list(
  "Sparse Log-Contrast" = cl,
  "trac (a = 1)" = trac, 
  "trac (a = 1/2)" = trac_weights
)

## Additional data sets for supplement:

# load cps_Mois:
load("../CentralParkSoil/cps_Mois_complasso_fixed_level_multi_splits.Rdata", e = cl <- new.env())
load("../CentralParkSoil/cps_Mois_trac_fixed_level_multi_splits.Rdata", e = trac <- new.env())
load("../CentralParkSoil/cps_Mois_trac_weights_fixed_level_multi_splits.Rdata", e = trac_weights <- new.env())
cl$cvfit <- cl$cvfit %>% map(~ .x %>% 
                               map(~ {
                                 .x$cv <- list(.x$cv)
                                 return(.x)
                               }))
methods_data[["Central Park Soil: Mois"]] <- list(
  "Sparse Log-Contrast" = cl,
  "trac (a = 1)" = trac, 
  "trac (a = 1/2)" = trac_weights
)

# load AGP:
load("../AmericanGut/complasso_fixed_level_multi_splits.Rdata", e = cl <- new.env())
load("../AmericanGut/trac_fixed_level_multi_splits.Rdata", e = trac <- new.env())
load("../AmericanGut/trac_weights_fixed_level_multi_splits.Rdata", e = trac_weights <- new.env())
cl$cvfit <- cl$cvfit %>% map(~ .x %>% 
                               map(~ {
                                 .x$cv <- list(.x$cv)
                                 return(.x)
                               }))
methods_data[["Gut (AGP): BMI"]] <- list(
  "Sparse Log-Contrast" = cl,
  "trac (a = 1)" = trac, 
  "trac (a = 1/2)" = trac_weights
)

# load Fram Strait (PA): Leucine
load("../Marine/marine_leucine_large_complasso_fixed_level_multi_splits.Rdata", e = cl <- new.env())
load("../Marine/marine_leucine_large_trac_fixed_level_multi_splits.Rdata", e = trac <- new.env())
load("../Marine/marine_leucine_large_trac_weights_fixed_level_multi_splits.Rdata", e = trac_weights <- new.env())
cl$cvfit <- cl$cvfit %>% map(~ .x %>% 
                               map(~ {
                                 .x$cv <- list(.x$cv)
                                 return(.x)
                               }))
methods_data[["Fram Strait (PA): Leucine"]] <- list(
  "Sparse Log-Contrast" = cl,
  "trac (a = 1)" = trac, 
  "trac (a = 1/2)" = trac_weights
)

# load Fram Strait (FL): Leucine
load("../Marine/marine_leucine_small_complasso_fixed_level_multi_splits.Rdata", e = cl <- new.env())
load("../Marine/marine_leucine_small_trac_fixed_level_multi_splits.Rdata", e = trac <- new.env())
load("../Marine/marine_leucine_small_trac_weights_fixed_level_multi_splits.Rdata", e = trac_weights <- new.env())
cl$cvfit <- cl$cvfit %>% map(~ .x %>% 
                               map(~ {
                                 .x$cv <- list(.x$cv)
                                 return(.x)
                               }))
methods_data[["Fram Strait (FL): Leucine"]] <- list(
  "Sparse Log-Contrast" = cl,
  "trac (a = 1)" = trac, 
  "trac (a = 1/2)" = trac_weights
)

both_1se <- list()
for (d in names(methods_data)) {
  both_1se_list <- list()
  for (i in names(methods_data[[d]])) {
    testerr <- methods_data[[d]][[i]]$testerr %>%
      map_dfr( ~ .x %>% 
                 bind_cols() %>% 
                 mutate(ilam = row_number()) %>% 
                 pivot_longer(-ilam, names_to = "Level", values_to = "test_error"),
               .id = "Split")
    i1se <- methods_data[[d]][[i]]$cvfit %>% 
      map_dfr(~.x %>% map_dbl(~ .x$cv[[1]]$i1se)) %>%
      mutate(Split = as.character(row_number())) %>% 
      pivot_longer(-Split, names_to="Level", values_to="i1se")
    testerr_1se <- testerr %>% left_join(i1se, by = c("Level", "Split")) %>% 
      filter(ilam == i1se)
    nnz <- methods_data[[d]][[i]]$nnz %>%
      map_dfr( ~ .x %>% 
                 bind_cols() %>% 
                 mutate(ilam = row_number()) %>% 
                 pivot_longer(-ilam, names_to = "Level", values_to = "nnz"),
               .id = "Split")
    nnz_1se <- nnz %>% left_join(i1se, by = c("Level", "Split")) %>% 
      filter(ilam == i1se)
    both_1se_list[[i]] <- testerr_1se %>% left_join(nnz_1se, by = c("ilam", "Level", "Split"))
  }
  both_1se[[d]] <- bind_rows(both_1se_list, .id = "Method")
}

all <- bind_rows(both_1se, .id = "Data") %>% select(-ilam, -i1se.x, -i1se.y)

tbl <- all %>% 
  filter(Level %in% lower_levels) %>% 
  group_by(Data, Method, Level) %>% 
  summarize(mean_1se_err = mean(test_error),
            mean_1se_nnz = mean(nnz)) %>% 
  ungroup() %>% 
  mutate(value = sprintf("%.2g (%s)", mean_1se_err, round(mean_1se_nnz))) %>% 
  mutate(Data = fct_relevel(factor(Data), names(methods_data)[c(1,3,2)])) %>% 
  arrange(Data) %>% 
  select(-mean_1se_err, -mean_1se_nnz) %>% 
  pivot_wider(all_of(c("Data", "Level")), names_from = Method, values_from = value) %>% 
  relocate(`Sparse Log-Contrast`, .after = last_col())

# get dimension of each data-level pair
dims <- map_dfr(methods_data, 
        ~ map_dfr(.x[[1]]$fit[[1]][lower_levels], ~ nrow(.x$beta)), .id = "Data") %>% 
  pivot_longer(-Data, names_to = "Level", values_to = "$p$")

tbl <- tbl %>% 
  left_join(dims, by = c("Data", "Level")) %>% 
  relocate(`$p$`, .after = Level)

dat_labs <- unique(tbl$Data)
dat_shortlabs <- c("sCD14", # ordered to match dat_labs
                   "cps_pH",
                   "tara_sal",
                   "cps_Mois",
                   "marine_leucine_small",
                   "marine_leucine_large",
                   "AGP")
for (i in seq_along(dat_labs)) {
  tbl %>%
    filter(Data == dat_labs[i]) %>% 
    select(-Data) %>% 
    mutate(Level = fct_relevel(factor(Level), c("OTU", "Genus", "Family"))) %>% 
    arrange(Level) %>% 
    rename(`Base Level` = Level) %>% 
    knitr::kable(format = "latex",
                 booktabs = TRUE,
                 escape = FALSE,
                 linesep = "") %>%
    kableExtra::kable_styling(latex_options = c("striped", "scale_down")) %>%
    kableExtra::row_spec(0, bold = TRUE) %>%
    str_remove("\\\\begin\\{table\\}") %>% 
    str_remove("\\\\end\\{table\\}") %>% 
    cat(file = sprintf("test_error_table/%s.tex", dat_shortlabs[i]))
}
