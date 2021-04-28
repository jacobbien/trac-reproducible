# Go from phyloseq object to something ready to input into trac
library(stringr)
library(dplyr)
library(Matrix)
library(tidyverse)
library(trac)
library(phyloseq)

agp <- readRDS("AGP_10.RDS")

y_factor <- sample_data(agp)$bmi
y <- as.numeric(levels(y_factor))[y_factor]
summary(y)
# filter 16 - 40 BMI (which excludes the most extreme categories
# of thinness and obesity. Note that those categories do not have
# lower/upper bounds)
# https://web.archive.org/web/20090418181049/http://www.who.int/bmi/index.jsp?introPage=intro_3.html
keep <- which(!is.na(y) & y <= 40 & y >= 16)

tax <- agp@tax_table@.Data
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
tax <- as.data.frame(tax, stringsAsFactors = TRUE)

# form phylo object:
tree1 <- tax_table_to_phylo(~Rank0/Rank1/Rank2/Rank3/Rank4/Rank5/Rank6/Rank7/OTU,
                            data = tax, collapse = TRUE)

# convert this to an A matrix to be used for aggregation:
A <- phylo_to_A(tree1)

dat <- list(y = y[keep],
            x = t(agp@otu_table@.Data[, keep]),
            tree = tree1,
            tax = tax,
            A = A,
            sample_data = as_tibble(sample_data(agp))[keep, ])
# rows of A correspond to OTUs as do columns of x
# rearrange columns of x to be in the order of rows of A:
dat$x <- dat$x[, match(str_match(rownames(A), "::([^:]+)$")[, 2],
                       colnames(dat$x))]
identical(str_match(rownames(A), "::([^:]+)$")[,2],
          colnames(dat$x))
saveRDS(dat, file = "AGP_processed.RDS")


# aggregate to higher levels (for comparison to log-contrast regression at
# fixed aggregation levels)
dat <- readRDS("AGP_processed.RDS")
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
saveRDS(dat_agg, file = "AGP_aggregated.RDS")
