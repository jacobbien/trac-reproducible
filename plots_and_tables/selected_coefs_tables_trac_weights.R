library(tidyverse)
library(Matrix)
method_name <- "trac_weights"
method_long_name <- "{\\tt trac} (a = 1/2)"

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
out <- map(trac_files,
           ~ {
             load(.x, envir = e <- new.env())
             return(list(fit = e$fit[[1]]$OTU, cvfit = e$cvfit[[1]]$OTU))
           })

tax_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "OTU")

ntop <- 15

latex_tbls <- map2(
  out,
  names(out),
  ~ {
    cat(.y, fill = TRUE)
    alpha_cv1se <- round(.x$fit[[1]]$alpha[, .x$cvfit$cv[[1]]$i1se], 6)
    names(alpha_cv1se) <- names(alpha_cv1se) %>% str_remove_all("[a-z]__")
    into <- tax_levels
    if (str_starts(names(alpha_cv1se)[1], "Life")) {
      into <- c(NA, into) # some start with Life others with Kingdom
      names(alpha_cv1se) <- if_else(str_count(names(alpha_cv1se), "::") == 7,
                                    str_replace(names(alpha_cv1se), # add missing Species
                                                "(::[^:]+)$", "::999\\1"),
                                    names(alpha_cv1se)
                                    )
    } else {
      names(alpha_cv1se) <- if_else(str_count(names(alpha_cv1se), "::") == 6,
                                    str_replace(names(alpha_cv1se), # add missing Species
                                                "(::[^:]+)$", "::999\\1"),
                                    names(alpha_cv1se)
                                    )
    }
    selected <- enframe(alpha_cv1se) %>%
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
                       str_remove_all(.y, "[^[:alnum:]]"),
                       "trac-weights",
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
      str_replace("value", "$\\\\alpha$")
  })

walk2(latex_tbls, short_names, ~ cat(.x, file = paste0("selected_tables/", method_name, "_", .y, ".tex")))
