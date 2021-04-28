# Import and prune American Gut Project (AGP) data    
# McDonald et al. (2018)

library(phyloseq)

# import biom file (can take a few minutes due to large file size):
ag <- import_biom("original/8237_analysis.biom") 
# import metadata from mapping file:
map <- read.delim("original/8237_analysis_mapping.txt", 
                  sep = "\t",
                  header = TRUE, 
                  row.names = 1)
# assign metadata to phyloseq object:
sample_data(ag) <- map

# All fecal data
ag.fe <- subset_samples(ag, body_habitat == "UBERON:feces") ## only fecal samples

ag.fe
## Full fecal phyloseq object
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 27116 taxa and 8440 samples ]
# sample_data() Sample Data:       [ 8440 samples by 517 sample variables ]
# tax_table()   Taxonomy Table:    [ 27116 taxa by 7 taxonomic ranks ]

## Prune samples
depths <- colSums(ag.fe@otu_table@.Data) ## calculate sequencing depths

## Pruning (Minimum sequencing depth: at least 10000 reads per sample)
ag.filt1 <- prune_samples(depths > 10000, ag.fe) 
ag.filt1
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 27116 taxa and 7203 samples ]
# sample_data() Sample Data:       [ 7203 samples by 517 sample variables ]
# tax_table()   Taxonomy Table:    [ 27116 taxa by 7 taxonomic ranks ]

## Pruning (taxa present in at least 10% of samples)
freq <- rowSums(sign(ag.filt1@otu_table@.Data))
ag.filt2 <- prune_taxa(freq > 0.1 * nsamples(ag.filt1), ag.filt1) 
ag.filt2
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 1387 taxa and 7203 samples ]
# sample_data() Sample Data:       [ 7203 samples by 517 sample variables ]
# tax_table()   Taxonomy Table:    [ 1387 taxa by 7 taxonomic ranks ]

## Save data to RDS file (taxa present in at least 10% of samples)
saveRDS(ag.filt2, file = "AGP_10.RDS")
