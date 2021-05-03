library(tidyverse)
library(trac)

small <- new.env()
large <- new.env()
load("../../Marine/marine_leucine_small_trac_fixed_level_multi_splits.Rdata", envir = small)
load("../../Marine/marine_leucine_large_trac_fixed_level_multi_splits.Rdata", envir = large)
small$dat <- readRDS("../../Marine/marine_leucine_small.RDS")
large$dat <- readRDS("../../Marine/marine_leucine_large.RDS")


# verify that the ordering of samples in train/test
# is what it should be:
stopifnot(small$dat$sample_data$leucine == small$dat$y)

# add a column of our predictions across entire data set:
i1se <- small$cvfit[[1]]$OTU$cv[[1]]$i1se
small$dat$sample_data$yhat <- NA
small$dat$sample_data$yhat[small$tr[[1]]] <- small$yhat_tr[[1]]$OTU[, i1se]
small$dat$sample_data$yhat[-small$tr[[1]]] <- small$yhat_te[[1]]$OTU[, i1se]

# verify that the ordering of samples in train/test
# is what it should be:
stopifnot(large$dat$sample_data$leucine == large$dat$y)

# add a column of our predictions across entire data set:
i1se <- large$cvfit[[1]]$OTU$cv[[1]]$i1se
large$dat$sample_data$yhat <- NA
large$dat$sample_data$yhat[large$tr[[1]]] <- large$yhat_tr[[1]]$OTU[, i1se]
large$dat$sample_data$yhat[-large$tr[[1]]] <- large$yhat_te[[1]]$OTU[, i1se]


both <- list(small$dat$sample_data,
             large$dat$sample_data) %>% 
  set_names(c("Free living", "Particle associated")) %>% 
  bind_rows(.id = "Type")
both %>%
  ggplot(aes(x = yhat,
             y = leucine,
             color = Region)) +
  geom_point(size = 3) + 
  geom_abline(slope = 1, intercept = 0) +
  labs(x = "Predicted Leucine",
       y = "Actual Leucine",
       title = "All Samples") +
  facet_wrap(~ Type) +
  theme(strip.text = element_text(size = 8)) + 
  coord_fixed() + 
  ggsave("marine.pdf", width = 7.5, height = 3)


cor(small$yhat_te[[1]]$OTU[, small$cvfit[[1]]$OTU$cv[[1]]$i1se], 
    small$dat$y[-small$tr[[1]]])
cor(large$yhat_te[[1]]$OTU[, large$cvfit[[1]]$OTU$cv[[1]]$i1se], 
    large$dat$y[-large$tr[[1]]])


