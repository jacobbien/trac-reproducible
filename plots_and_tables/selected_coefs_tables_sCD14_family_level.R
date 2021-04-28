library(tidyverse)
library(Matrix)

method_names <- c("trac", "trac_weights")
method_long_names <- paste("{\\tt trac} on family level", 
                           c("(a = 1)", "(a = 1/2)"))
trac_files <- paste0("../sCD14_HIV/", 
                     method_names,
                     "_fixed_level_multi_splits.Rdata")
data_name <- "Gut (HIV): sCD14"
short_data_name <- "sCD14"

out <- map(trac_files,
           ~ {
             load(.x, envir = e <- new.env())
             return(list(fit = e$fit[[1]]$Family, cvfit = e$cvfit[[1]]$Family))
           })

tax_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "OTU")

ntop <- 15

for (i in 1:2) {
  cat(method_long_names[i], fill = TRUE)
  alpha_cv1se <- round(out[[i]]$fit[[1]]$alpha[, out[[i]]$cvfit$cv[[1]]$i1se], 6)
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
                            ntop, method_long_names[i], data_name),
                    sprintf("Coefficients selected by %s for %s", 
                            method_long_names[i], data_name))
  # add latex label:
  caption <- sprintf("\\label{fig:%s-family-level-%s} %s",
                     short_data_name,
                     method_names[i],
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
    str_replace("value", "$\\\\alpha$") %>% 
    cat(file = paste0("selected_tables/", short_data_name, "_family-level_", method_names[i], ".tex"))
}

