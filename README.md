# All code for trac paper

This directory contains the code needed to fully reproduce all the data examples in the manuscript [Tree-Aggregated Predictive Modeling of Microbiome Data](https://www.biorxiv.org/content/10.1101/2020.09.01.277632v1.full).  It depends on the R package `trac`, whose latest version is available [here](https://github.com/jacobbien/trac).

The data set directories are the following:

- AmericanGut

- CentralParkSoil

- Marine

- sCD14_HIV

- Tara

Each data set directory has a common structure and an identical workflow (with step 0 sometimes not being necessary).

**Directory structure:**

- original directory: contains the files from others that is our starting point

- numbered .R files: this is our workflow, detailed below. 

- .RDS and .Rdata files: these are created by running the .R files.

**Workflow:**

*Step 0:* Create phyloseq object (see file *0create_phyloseq_object.R*)

- in some cases 

- some (generally minimal) amount of filtering of samples  and/or OTUs

*Step 1:* Prepare Data (see file *1prep_data_all_levels.R*)

 - starts by reading in phyloseq object

 - gets the relevant subset of the samples (e.g. response should not be `NA`)

 - goes from the tax table to a phylo object to the `A` matrix used in `trac` 

 - saves a .RDS file that has a list of objects to be used by `trac`

 - in particular, this list contains `y`, `x`, `A` (and also the tree as a `phylo` object and the tax table)

 - performs fixed-level aggregation (by simple summation) so that `trac` can be applied at different base levels

 - saves an "aggregated" .RDS file that has the data ready for `trac` at each base level

*Step 2:* Fit `trac` model (see file *2trac_multi_splits.R*)

 - repeats for 10 random splits into train / test sets and for each base level...

 - fits `trac` (a=1)

 - performs 5-fold CV

 - saves `trac` and `cvtrac` objects as well as other information to file

*Step 3:* Identical to Step 2 but for trac (a=1/2) (see file *3trac_weights_multi_splits.R*)

*Step 4:* Identical to Step 2 but for sparse log contrast (*4complasso_multi_splits.R*)

In addition to the data set directories described above there is also a directory called `plots_and_tables`.  The .R files in this directory create all of the plots and tables in the paper and supplement based on the outputs of the data workflows described above.
