# we predict Leucine (based on Figure 7 of Fadeev et al) 
# ... which measures (according to them) bacterial production

library(phyloseq)
library(tidyverse)
library(trac)

BAC_pruned <- readRDS("original/BAC_pruned.rds")
# See https://github.com/edfadeev/Bact-comm-PS85/blob/master/dataset_preprocess.R
# for code to generate

# subset to large size fractions:
ii_large <- sample_data(BAC_pruned) %>% 
  as_tibble() %>% 
  mutate(row = row_number()) %>% 
  filter(Fraction == 3) %>% 
  pull(row)

sample_data(BAC_pruned) <- sample_data(BAC_pruned)[ii_large, ]

filter.OTU <- genefilter_sample(BAC_pruned,
                                filterfun_sample(function(x) x >= 1),
                                A = 0.1 * nsamples(BAC_pruned))
# Apply filter  
BAC_pruned.filter <- prune_taxa(filter.OTU, BAC_pruned)

y <- sample_data(BAC_pruned.filter)$leucine

tax <- BAC_pruned.filter@tax_table@.Data
colnames(tax)[1] <- "Kingdom"
# replace "unclassified" with the appropriate blank tag
blank <- paste0(c("k", "p", "c", "o", "f", "g"), "__")
tax[str_detect(tax, "_unclassified")] <- "Unclassified"
for (i in 1:6) tax[tax[, i] == "Unclassified", i] <- blank[i]

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
# cumulative labels are harder to read but easier to work with:
for (i in 2:7) {
  tax[, i] <- paste(tax[, i-1], tax[, i], sep = "::")
}
tax <- as.data.frame(tax, stringsAsFactors = TRUE)

# form phylo object:
tree1 <- tax_table_to_phylo(~Kingdom/Phylum/Class/Order/Family/Genus/OTU, 
                            data = tax, collapse = TRUE)

# convert this to an A matrix to be used for aggregation:
A <- phylo_to_A(tree1)
imiss <- which(is.na(y))

dat <- list(y = y[-imiss],
            x = t(BAC_pruned.filter@otu_table@.Data[, -imiss]),
            tree = tree1, 
            tax = tax,
            A = A,
            sample_data = as_tibble(sample_data(BAC_pruned.filter))[-imiss, ])
# rows of A correspond to OTUs as do columns of x
# rearrange columns of x to be in the order of rows of A:
dat$x <- dat$x[, match(str_match(rownames(A), "::([^:]+)$")[, 2],
                       colnames(dat$x))]
identical(str_match(rownames(A), "::([^:]+)$")[,2], 
          colnames(dat$x))
saveRDS(dat, file = "marine_leucine_large.RDS")


# aggregate to higher levels (for comparison to log-contrast regression at
# fixed aggregation levels)
dat <- readRDS("marine_leucine_large.RDS")
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
                                                  level = i + 1,
                                                  collapse = TRUE)
}
dat_agg[["OTU"]] <- dat
saveRDS(dat_agg, file = "marine_leucine_large_aggregated.RDS")
