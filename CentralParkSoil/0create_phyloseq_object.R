# This is a modification of the file "original/Data Analysis Pipeline for Soil PH Dataset.R"
# from Washburne, in which we instead create a phyloseq object
library(ape)
library(tidyverse)
library(phyloseq)

# OTU Table
Data <- read.csv('original/otu_table_wTax_40000_filt2.txt',
                 sep = '\t',
                 row.names = 1,
                 quote = '#')
taxonomy <- Data[, dim(Data)[2]]
Data <- Data[, 1:(dim(Data)[2] - 1)] # OTU by sample

otus <- rownames(Data) %>% str_replace('\"', '') %>% str_replace('\"', '')
rownames(Data) <- otus
colnames(Data) <- colnames(Data) %>% str_remove("^X.") %>% str_remove(".$")

taxonomy <- taxonomy %>% str_remove_all('\"')
Tax <- tibble(taxonomy = taxonomy) %>% 
  separate(col = taxonomy,
                 into = c("Kingdom", 
                          "Phylum",
                          "Class",
                          "Order",
                          "Family",
                          "Genus",
                          "Species"),
                 sep = "; ") %>% 
  as.matrix()
rownames(Tax) <- otus

# Sample data
MAP <- read.csv('original/CP_map.txt', sep = '\t', header = TRUE)
sample_names <- as.character(MAP$X.SampleID)
rownames(MAP) <- sample_names

# reorder columns of Data to match order in MAP:
Data <- Data[, sample_names]

# Phylogeny
tree <- read.tree('original/rep_set_aligned_gapfilt20_filt2.tre')
otus <- tree$tip.label

# reorder rows of Data and Tax to match order in tree tip labels:
Data <- Data[tree$tip.label, ]
Tax <- Tax[tree$tip.label, ]

# create phyloseq object:
sampledata <- sample_data(MAP)
otu <- otu_table(Data, taxa_are_rows = TRUE)
tax <- tax_table(Tax)
cps <- phyloseq(sampledata, otu, tax, tree)
saveRDS(cps, file = "cps.RDS")

