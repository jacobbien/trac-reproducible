library(tidyverse)
library(trac)

dat <- readRDS("../../AmericanGut/AGP_processed.RDS")
load("../../AmericanGut/trac_fixed_level_multi_splits.Rdata")

# show the y vs yhat plot for the first split:
tibble(yte = dat$y[-tr[[1]]],
       yhat_te = yhat_te[[1]]$OTU[, cvfit[[1]]$OTU$cv[[1]]$i1se]
) %>% 
  ggplot(aes(x = yhat_te, y = yte)) +
  geom_point(alpha = 0.1) + 
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(title = "Test Set Performance", x = "Predicted BMI", y = "Actual BMI") +
  coord_fixed() + 
  ggsave(filename = "bmi.pdf")

cor(dat$y[-tr[[1]]], yhat_te[[1]]$OTU[, cvfit[[1]]$OTU$cv[[1]]$i1se])

