# Go from phyloseq object to something ready to input into trac

library(phyloseq)
library(Matrix)
library(tidyverse)
library(trac)

cps <- readRDS("cps.RDS")

y <- sample_data(cps)$pH
summary(y)

tax <- cps@tax_table@.Data

# replace "unclassified" with the appropriate blank tag
blank <- paste0(c("k", "p", "c", "o", "f", "g", "s"), "__")
for (i in 1:7) tax[tax[, i] == "unclassified", i] <- blank[i]

tax <- cbind("Life", tax); colnames(tax)[1] <- "Rank0"
# add an OTU column
tax <- cbind(tax, rownames(tax))
colnames(tax)[ncol(tax)] <- "OTU"

# make it so labels are unique
for (i in seq(2, 8)) {
  # add a number when the type is unknown... e.g. "g__"
  ii <- nchar(tax[, i]) == 3
  if (sum(ii) > 0)
    tax[ii, i] <- paste0(tax[ii, i], 1:sum(ii))
}
# cumulative labels are harder to read but easier to work with:
for (i in 2:9) {
  tax[, i] <- paste(tax[, i-1], tax[, i], sep = "::")
}
tax <- data.frame(tax, stringsAsFactors = TRUE)

# form phylo object:
tree1 <- tax_table_to_phylo(~Rank0/Kingdom/Phylum/Class/Order/Family/Genus/Species/OTU,
                            data = tax, collapse = TRUE)

# convert this to an A matrix to be used for aggregation:
A <- phylo_to_A(tree1)

dat <- list(y = y,
            x = t(cps@otu_table@.Data),
            tree = tree1,
            tax = tax,
            A = A,
            sample_data = as_tibble(sample_data(cps)))
# rows of A correspond to OTUs as do columns of x
# rearrange columns of x to be in the order of rows of A:
dat$x <- dat$x[, match(str_match(rownames(A), "::([^:]+)$")[, 2],
                       colnames(dat$x))]
identical(str_match(rownames(A), "::([^:]+)$")[,2],
          colnames(dat$x))
saveRDS(dat, file = "cps_pH_processed.RDS")

# For Moisture
keep <- which(!is.na(dat$sample_data$Moisture))
dat$y <- dat$sample_data$Moisture
dat$y <- dat$y[keep]
dat$x <- dat$x[keep, ]
dat$sample_data <- dat$sample_data[keep, ]
saveRDS(dat, file = "cps_Mois_processed.RDS")



# aggregate to higher levels (for comparison to log-contrast regression at
# fixed aggregation levels)
dat <- readRDS("cps_pH_processed.RDS")
level_names <- c("Phylum",
                 "Class",
                 "Order",
                 "Family",
                 "Genus",
                 "Species",
                 "OTU")
dat_agg <- list()
for (i in seq(1,6)) {
  dat_agg[[level_names[i]]] <- aggregate_to_level(x = dat$x,
                                                  y = dat$y, 
                                                  A = dat$A,
                                                  tax = dat$tax, 
                                                  level = i + 2,
                                                  collapse = TRUE)
}
dat_agg[["OTU"]] <- dat
saveRDS(dat_agg, file = "cps_pH_aggregated.RDS")

# Moisture
dat <- readRDS("cps_Mois_processed.RDS")
dat_agg <- list()
for (i in seq(1,6)) {
  dat_agg[[level_names[i]]] <- aggregate_to_level(x = dat$x,
                                                  y = dat$y, 
                                                  A = dat$A,
                                                  tax = dat$tax, 
                                                  level = i + 2,
                                                  collapse = TRUE)
}
dat_agg[["OTU"]] <- dat
saveRDS(dat_agg, file = "cps_Mois_aggregated.RDS")

