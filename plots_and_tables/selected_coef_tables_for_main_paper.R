# For the main paper (as opposed to the appendix), we want to embed the tabular
# for the selected coefficient tables within a figure.  Therefore, this script
# removes the \begin{table} \end{table} and caption.

library(tidyverse)

files <- c("selected_tables/trac_sCD14.tex",
          "selected_tables/trac_cps_pH.tex",
          "selected_tables/trac_tara.tex")

for (file in files) {
  tab <- read_file(file)
  tab %>% 
    str_remove("^.+\\n\\n\\\\caption.+\\n") %>%
    str_remove("\\\\end\\{table\\}") %>% 
    write_file(file = file %>%
                 str_replace("/", "/main/") %>% 
                 str_replace(".tex$", "_main.tex"))
}
