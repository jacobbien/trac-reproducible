# make the yhat vs y plot for tara, tax (a=1)

library(tidyverse)
library(trac)

dat_list <- readRDS("../../Tara/tara_sal_processed_aggregated.RDS")
dat <- dat_list$OTU
rm(dat_list)

load("../../Tara/tara_sal_trac_fixed_level_multi_splits.Rdata")
ytr <- dat$y[tr[[1]]]
yte <- dat$y[-tr[[1]]]
# verify that the ordering of samples in train/test
# is what it should be:
stopifnot(dat$sample_data$Mean_Salinity..PSU..[tr[[1]]] == ytr)
stopifnot(dat$sample_data$Mean_Salinity..PSU..[-tr[[1]]] == yte)
# add a column of our predictions across entire data set:
dat$sample_data$yhat <- NA
dat$sample_data$yhat[tr[[1]]] <- yhat_tr[[1]]$OTU[, cvfit[[1]]$OTU$cv[[1]]$i1se]
dat$sample_data$yhat[-tr[[1]]] <- yhat_te[[1]]$OTU[, cvfit[[1]]$OTU$cv[[1]]$i1se]

dat$sample_data %>%
  rename(`Longhurst Province` = Marine.pelagic.biomes..Longhurst.2007.) %>% 
  ggplot(aes(x = yhat,
             y = Mean_Salinity..PSU..,
             color = `Longhurst Province`)) +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) +
  labs(x = "Predicted Salinity",
       y = "Actual Salinity",
       title = "All Samples") + 
  coord_fixed() +
  ggsave("tara.pdf")
cor(yhat_te[[1]]$OTU[, cvfit[[1]]$OTU$cv[[1]]$i1se], yte)
