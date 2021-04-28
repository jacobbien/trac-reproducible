library(Matrix)
library(tidyverse)
library(trac)
library(reticulate)

dat_list <- readRDS("AGP_aggregated.RDS")
level_names <- rev(names(dat_list))
log_pseudo <- function(x, pseudo_count = 1) log(x + pseudo_count)

# apply trac weights on 10 random train-test splits...
set.seed(123)
nsplit <- 10
ntot <- length(dat_list$OTU$y)
n <- round(2/3 * ntot)

tr <- list()
fit <- list()
cvfit <- list()
yhat_tr <- list()
yhat_te <- list()
trainerr <- list()
testerr <- list()
nnz <- list()
for (j in seq(nsplit)) {
  cat("split", j, fill = TRUE)
  tr[[j]] <- sample(ntot, n)
  fit[[j]] <- list()
  cvfit[[j]] <- list()
  yhat_tr[[j]] <- list()
  yhat_te[[j]] <- list()
  trainerr[[j]] <- list()
  testerr[[j]] <- list()
  nnz[[j]] <- list()
  for (i in level_names) {
    cat(i, fill = TRUE)
    ytr <- dat_list[[i]]$y[tr[[j]]]
    yte <- dat_list[[i]]$y[-tr[[j]]]
    ztr <- log_pseudo(dat_list[[i]]$x[tr[[j]], ])
    zte <- log_pseudo(dat_list[[i]]$x[-tr[[j]], ])
    fit[[j]][[i]] <- trac(ztr, ytr, A = dat_list[[i]]$A, min_frac = 1e-3,
                          nlam = 30, w = colSums(dat_list[[i]]$A)^0.5)
    cvfit[[j]][[i]] <- cv_trac(fit[[j]][[i]], Z = ztr, y = ytr, A = dat_list[[i]]$A)
    yhat_tr[[j]][[i]] <- predict_trac(fit[[j]][[i]], new_Z = ztr)[[1]]
    yhat_te[[j]][[i]] <- predict_trac(fit[[j]][[i]], new_Z = zte)[[1]]
    trainerr[[j]][[i]] <- colMeans((yhat_tr[[j]][[i]] - ytr)^2)
    testerr[[j]][[i]] <- colMeans((yhat_te[[j]][[i]] - yte)^2)
    nnz[[j]][[i]] <- colSums(fit[[j]][[i]][[1]]$gamma != 0)
  }
}
# verify that these are same train-test splits as for trac:
load("trac_fixed_level_multi_splits.Rdata", envir = trac <- new.env())
stopifnot(identical(tr, trac$tr))
save(tr, fit, cvfit, yhat_tr, yhat_te, trainerr, testerr, nnz,
     file = "trac_weights_fixed_level_multi_splits.Rdata")
