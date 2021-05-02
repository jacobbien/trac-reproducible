library(tidyverse)
library(trac)
library(gridExtra)

dat <- readRDS("../sCD14_HIV/sCD14.RDS")
load("../sCD14_HIV/trac_fixed_level_multi_splits.Rdata")
yhat_tr <- list(yhat_tr[[1]]$OTU)
yhat_te <- list(yhat_te[[1]]$OTU)
testerr <- testerr[[1]]$OTU
cvfit <- cvfit[[1]]$OTU
fit <- fit[[1]]$OTU
yte <- dat$y[-tr[[1]]]
ytr <- dat$y[tr[[1]]]

g <- list()

# cv curve
g[[1]] <- tibble(mean_cv = cvfit$cv[[1]]$m,
                 lower_cv = cvfit$cv[[1]]$m - cvfit$cv[[1]]$se,
                 upper_cv = cvfit$cv[[1]]$m + cvfit$cv[[1]]$se,
                 non_zero = cvfit$cv[[1]]$nonzeros
) %>% 
  ggplot(aes(x = non_zero, y = mean_cv, ymin = lower_cv, ymax = upper_cv)) +
  geom_point() +
  geom_line() +
  geom_errorbar() +
  geom_vline(xintercept = cvfit$cv[[1]]$nonzeros[c(cvfit$cv[[1]]$i1se,
                                                   cvfit$cv[[1]]$ibest)],
             lty = 2:3) +
  labs(title = "B: Cross Validation", y = "CV Error", x = "Number of nonzero alpha")

# coefficient paths
g[[2]] <- fit[[1]]$alpha %>%
  t() %>% 
  as_tibble() %>% 
  mutate(frac = fit[[1]]$fraclist) %>% 
  pivot_longer(cols = -frac) %>% 
  ggplot(aes(x = frac, y = value, group = name, color = name)) + 
  geom_line() + 
  scale_x_log10() + 
  theme(legend.position = "none") + 
  geom_vline(xintercept = cvfit$cv[[1]]$fraclist[c(cvfit$cv[[1]]$i1se,
                                                   cvfit$cv[[1]]$ibest)],
             lty = 2:3) +
  labs(title = "A: Solution Path", y = "alpha", x = "Fraction of lambda_max")

# test error curve
g[[3]] <- tibble(testerr = testerr,
                 non_zero = cvfit$cv[[1]]$nonzeros
) %>% 
  ggplot(aes(x = non_zero, y = testerr)) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = cvfit$cv[[1]]$nonzeros[c(cvfit$cv[[1]]$i1se,
                                                   cvfit$cv[[1]]$ibest)],
             lty = 2:3) + 
  labs(title = "D: Error on Test Set", y = "Test Error", x = "Number of nonzero alpha")

g[[4]] <- tibble(yte = yte,
                 yhat_te = yhat_te[[1]][, cvfit$cv[[1]]$i1se]
) %>% 
  ggplot(aes(x = yhat_te, y = yte)) +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) +
  labs(title = "Test Set Performance", x = "Predicted sCD14", y = "Actual sCD14")

g[[5]] <- tibble(yte = yte,
                 `CV 1SE rule` = yhat_te[[1]][, cvfit$cv[[1]]$i1se],
                 `CV Best` = yhat_te[[1]][, cvfit$cv[[1]]$ibest],
) %>% 
  pivot_longer(cols = 2:3, names_to = "lambda") %>% 
  ggplot(aes(x = value, y = yte, color = lambda)) +
  geom_point(size = 2) + 
  geom_abline(slope = 1, intercept = 0) +
  labs(title = "C: Test Set Performance", x = "Predicted sCD14", y = "Actual sCD14") +
  # theme(legend.position = "bottom") +
  coord_equal()

cor(yte, yhat_te[[1]][, cvfit$cv[[1]]$ibest])
cor(ytr, yhat_tr[[1]][, cvfit$cv[[1]]$ibest])
cor(yte, yhat_te[[1]][, cvfit$cv[[1]]$i1se])
cor(ytr, yhat_tr[[1]][, cvfit$cv[[1]]$i1se])

gr <- g[c(2,1,5,3)] %>% 
map(ggplotGrob) %>% 
  arrangeGrob(grobs = ., ncol = 2)
ggsave(gr, filename = "sCD14.pdf", width = 10, height = 6)

