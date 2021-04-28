library(tidyverse)
library(Matrix)
dat_files <- paste0("../",
                    c("sCD14_HIV/sCD14",
                      "AmericanGut/AGP_processed",
                      "CentralParkSoil/cps_pH_processed",
                      "CentralParkSoil/cps_Mois_processed",
                      "Marine/marine_leucine_large",
                      "Marine/marine_leucine_small",
                      "Tara/tara_sal_processed"), 
                    ".RDS")
trac_files <- paste0("../",
                     c("sCD14_HIV/",
                       "AmericanGut/",
                       "CentralParkSoil/cps_pH_",
                       "CentralParkSoil/cps_Mois_",
                       "Marine/marine_leucine_large_",
                       "Marine/marine_leucine_small_",
                       "Tara/tara_sal_"),
                     "trac_fixed_level_multi_splits.Rdata")
names(trac_files) <- c("Gut (HIV): sCD14",
                       "Gut (AGP): BMI",
                       "Central Park Soil: pH",
                       "Central Park Soil: Mois",
                       "Fram Strait (PA): Leucine",
                       "Fram Strait (FL): Leucine",
                       "Ocean (TARA): Salinity")

out <- map2(trac_files, 
            dat_files, 
            ~ {
              cat(.x, fill = TRUE)
              load(.x, envir = e <- new.env())
              dat <- readRDS(.y)
              return(list(fit = map(e$fit, ~ .x$OTU),
                          cvfit = map(e$cvfit, ~ .x$OTU),
                          n = nrow(dat$x),
                          p = ncol(dat$x)))
              })

tax_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "OTU")
tax_levels <- factor(tax_levels, levels = tax_levels)

tbl_1se_all <- out %>%
  map_dfr(~ { # for each data set...
    map2_dfr(.x$fit, .x$cvfit, ~ { # for each split...
      which(round(.x[[1]]$gamma[, .y$cv[[1]]$i1se], 6) != 0) %>% 
        names() %>% 
        str_remove("^Life::") %>% 
        str_count("::") %>% 
        {table(tax_levels[. + 1])}      
    }, .id = "Split")
  }, .id = "Data")

np_df <- out %>% map_dfr(~ tibble(n = .x$n, p = .x$p), .id = "Data")

tbl_1se <- np_df %>% 
  left_join(tbl_1se_all %>%
              group_by(Data) %>%
              summarize(across(Kingdom:OTU, mean)),
            by = "Data")
names(tbl_1se)[2:3] <- c("$n$", "$p$")
tbl_1se %>% 
  knitr::kable(format = "latex", booktabs = TRUE, escape = FALSE) %>% 
  kableExtra::kable_styling(latex_options = c("striped", "scale_down")) %>% 
  kableExtra::row_spec(0, bold = TRUE) %>% 
  cat(file = "meta_table.tex")

read_file("meta_table.tex") %>% 
  str_remove("\\\\begin\\{table\\}") %>% 
  str_remove("\\\\end\\{table\\}") %>% 
  write_file(file = "meta_table2.tex")
