library(tidyverse)
library(trac)
library(phyloseq)

# Data provided as phyloseq object by Rivera-Pinto (personal communication)

load("original/Data_HIV.RData")
# Less extreme filtering compared to J. Rivera-Pinto (see procSCD14.R)

# Build a filter (at OTU level)
filter.OTU <- genefilter_sample(x, filterfun_sample(function(x) x >= 1),
                                A = 0.1 * nsamples(x))
# Apply filter
x.filter <- prune_taxa(filter.OTU, x)
tax <- x.filter@tax_table@.Data

# replace "unclassified" with the appropriate blank tag
blank <- paste0(c("k", "p", "c", "o", "f", "g"), "__")
for (i in 1:6) tax[tax[, i] == "unclassified", i] <- blank[i]

# add an OTU column
tax <- cbind(tax, rownames(tax))
colnames(tax)[7] <- "OTU"

# make it so labels are unique
for (i in seq(2, 6)) {
  # add a number when the type is unknown... e.g. "g__"
  ii <- nchar(tax[, i]) == 3
  if (sum(ii) > 0)
    tax[ii, i] <- paste0(tax[ii, i], 1:sum(ii))
}

# add a root node that combines the three kingdoms into a tree:
tax <- as.data.frame(tax)
tax$Rank0 <- rep("Life", nrow(tax))
tax <- tax[, c(8, 1:7)]

# cumulative labels are harder to read but easier to work with:
for (i in 2:ncol(tax)) {
  tax[, i] <- paste(tax[, i-1], tax[, i], sep = "::")
}
# convert all columns from character to factors for tax_table_to_phylo
for (i in seq_along(tax)) tax[, i] <- factor(tax[, i])


# form phylo object:
tree1 <- tax_table_to_phylo(~Rank0/Kingdom/Phylum/Class/Order/Family/Genus/OTU,
                            data = tax, collapse = TRUE)

# convert this to an A matrix to be used for aggregation:
A <- phylo_to_A(tree1)
y <- sample_data(x.filter)$sCD14
yy <- as.numeric(levels(y))[y]
imiss <- which(is.na(yy)) # four samples have missing sCD14
# we have n = 152 samples whereas paper has only 151 for some reason.
dat <- list(y = yy[-imiss],
            x = t(x.filter@otu_table@.Data[, -imiss]),
            tree = tree1,
            tax = tax,
            A = A)
dat$x <- dat$x[, match(str_match(rownames(A), "::([^:]+)$")[,2],
                       colnames(dat$x))]
identical(str_match(rownames(A), "::([^:]+)$")[,2],
          colnames(dat$x))
saveRDS(dat, file = "sCD14.RDS")


# aggregate to higher levels using conventional fixed-level aggregation:

dat <- readRDS("sCD14.RDS")
level_names <- c("Phylum",
                 "Class",
                 "Order",
                 "Family",
                 "Genus",
                 "OTU")
dat_agg <- list()
for (i in seq(1,5)) {
  dat_agg[[level_names[i]]] <- aggregate_to_level(x = dat$x,
                                                  y = dat$y, 
                                                  A = dat$A,
                                                  tax = dat$tax, 
                                                  level = i + 2,
                                                  collapse = TRUE)
}
dat_agg[["OTU"]] <- dat
saveRDS(dat_agg, file = "sCD14_aggregated.RDS")


