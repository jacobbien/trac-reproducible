library(tidyverse)
library(Matrix)
method_name <- "complasso"
method_long_name <- "the sparse log-contrast method"

trac_files <- paste0("../",
                     c(
                       "sCD14_HIV/",
                       "AmericanGut/",
                       "CentralParkSoil/cps_pH_",
                       "CentralParkSoil/cps_Mois_",
                       "Marine/marine_leucine_large_",
                       "Marine/marine_leucine_small_",
                       "Tara/tara_sal_"
                        ), 
                     method_name,
                     "_fixed_level_multi_splits.Rdata")
names(trac_files) <- c("Gut (HIV): sCD14",
                       "Gut (AGP): BMI",
                       "Central Park Soil: pH",
                       "Central Park Soil: Mois",
                       "Fram Strait (PA): Leucine",
                       "Fram Strait (FL): Leucine",
                       "Ocean (TARA): Salinity")
short_names <- c("sCD14", 
                 "agp",
                 "cps_pH",
                 "cps_Mois", 
                 "marine_large", 
                 "marine_small",
                 "tara")

dat_files <- c(
  "../sCD14_HIV/sCD14_aggregated.RDS",
  "../AmericanGut/AGP_aggregated.RDS",
  "../CentralParkSoil/cps_pH_aggregated.RDS",
  "../CentralParkSoil/cps_Mois_aggregated.RDS",
  "../Marine/marine_leucine_large_aggregated.RDS",
  "../Marine/marine_leucine_small_aggregated.RDS",
  "../Tara/tara_sal_processed_aggregated.RDS"
)

out <- map2(trac_files, dat_files,
           ~ {
             load(.x, envir = e <- new.env())
             dat <- readRDS(.y)
             cvfit <- e$cvfit[[1]]$OTU
             cvfit$cv <- list(cvfit$cv)
             list(fit = list(e$fit[[1]]$OTU), 
                  cvfit = cvfit,
                  A = dat$OTU$A,
                  x = dat$OTU$x)
           })

tax_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "OTU")

ntop <- 15

latex_tbls <- map2(
  out,
  names(out),
  ~ {
    cat(.y, fill = TRUE)
    otu_names <- tibble(full = rownames(.x$A),
                        short = rownames(.x$A) %>% str_extract("[[:alnum:]_]+$"))
    beta_cv1se <- round(.x$fit[[1]]$beta[, .x$cvfit$cv[[1]]$i1se], 6)
    stopifnot(names(beta_cv1se) == otu_names$short) # make sure order matches
    names(beta_cv1se) <- otu_names$full %>% str_remove_all("[a-z]__")
    into <- tax_levels
    if (str_starts(names(beta_cv1se)[1], "Life")) {
      into <- c(NA, into) # some start with Life others with Kingdom
      if (str_count(names(beta_cv1se)[1], "::") == 7) {
        # missing Species
        names(beta_cv1se) <- names(beta_cv1se) %>% 
          str_replace("(::[^:]+)$", "::999\\1")
      }
    } else {
      if (str_count(names(beta_cv1se)[1], "::") == 6) {
        # missing Species
        names(beta_cv1se) <- names(beta_cv1se) %>% 
          str_replace("(::[^:]+)$", "::999\\1")
      }
    }
    selected <- enframe(beta_cv1se) %>%
      filter(value != 0) %>% 
      separate(name,
               into = into,
               sep = "::",
               remove = FALSE,
               fill = "right")
    caption <- ifelse(nrow(selected) > ntop,
                      sprintf("Top %s coefficients selected by %s for %s",
                              ntop, method_long_name, .y),
                      sprintf("Coefficients selected by %s for %s", 
                         method_long_name,
                         .y))
    # add latex label:
    caption <- sprintf("\\label{fig:%s-%s} %s", 
                       str_remove_all(.y,"[^[:alnum:]]"),
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
      str_replace("value", "$\\\\beta$")
  })

walk2(latex_tbls, short_names, ~ cat(.x, file = paste0("selected_tables/", method_name, "_", .y, ".tex")))
