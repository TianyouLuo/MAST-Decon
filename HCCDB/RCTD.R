library(spacexr)
library(Matrix)
library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

slide = args[1]
reffile = args[2]

stdir = "/proj/yunligrp/users/tianyou/data_ST/HCCDB/ST/expression/cleaned/"
scRNAdir = "/proj/yunligrp/users/tianyou/data_ST/HCCDB/sc/processed/"
RCTDdir = paste0("/proj/yunligrp/users/tianyou/MASTDecon/HCCDB/results_", reffile, "/RCTD/")

print(RCTDdir)
print(paste0("Slide:", slide))
if (!file.exists(RCTDdir)){
  dir.create(RCTDdir)
}

options(check.name = F) ## to avoid changing numerical gene names
##### read in single cell and ST data ####
ref = fread(file.path(scRNAdir, paste0(reffile, "_filtered.txt")))
cluster = fread(file.path(scRNAdir, paste0(reffile, "_meta.csv")))
spots <- fread(file.path(stdir, paste0(slide,"_count.tsv"))) # load in counts matrix


## convert ST data to the required format: gene on rows and barcodes on columns
spots = as.data.frame(spots, check.names = FALSE)
rownames(spots) = spots[,1]
spots[,1] = NULL

##### convert count df into matrix #####
ref_mat = as.data.frame(ref, check.names = FALSE)
rownames(ref_mat) <- ref_mat[,1]
ref_mat[,1] <- NULL # Move first column to rownames

# keep only overlap genes with minimum count > 0 
spots = spots[rowSums(spots) > 0,]
overlap_gene = intersect(rownames(spots), rownames(ref_mat))
ref_overlap = ref_mat[overlap_gene,]
spots = spots[overlap_gene,]


##### clean the cluster labels #####
# In the step creating the reference count matrix, 
# we already filtered on celltypes and only kept celltypes with enough cells
cell_types = cluster$Celltype
names(cell_types) <- cluster$barcode # create cell_types named list
cell_types <- as.factor(cell_types) # convert to factor data type


##### convert count df into matrix #####
nUMI <- colSums(ref_overlap)
names(nUMI) <- colnames(ref_overlap) # create nUMI named list

## create reference data object
reference <- Reference(ref_overlap, cell_types, nUMI)

##### Create Spatial RNA object #####
coords <- fread(file.path(stdir, paste0(slide, "_meta.csv")))
coords = coords[coords$barcode %in% colnames(spots), c("barcode", "imagerow", "imagecol")]
coords = as.data.frame(coords)
colnames(coords) = c("barcodes", "xcoord", "ycoord")
rownames(coords) <- coords$barcodes
coords$barcodes <- NULL # Move barcodes to rownames
coords$xcoord = as.numeric(coords$xcoord)
coords$ycoord = as.numeric(coords$ycoord)
nUMI_spot <- colSums(spots) # In this case, total counts per pixel is nUMI

### Create SpatialRNA object
puck <- SpatialRNA(coords, spots, nUMI_spot)


## create RCTD object
# myRCTD <- create.RCTD(puck, reference, max_cores = 1,
#                       gene_cutoff = 1e-4, fc_cutoff = 0.5,
#                       gene_cutoff_reg = 2e-4, fc_cutoff_reg = 0.7)
myRCTD <- create.RCTD(puck, reference, max_cores = 1)

#### run RCTD ####
myRCTD@config$N_epoch = 20 ## the number of iterations that choose_sigma_c takes
myRCTD_full <- run.RCTD(myRCTD, doublet_mode = 'full')
saveRDS(myRCTD_full, file.path(RCTDdir, paste0(slide,"_RCTD.rds")))


#### summarize results ####
RCTD_results_full <- myRCTD_full@results
RCTD_norm_weights = as.matrix(sweep(RCTD_results_full$weights, 1, rowSums(RCTD_results_full$weights), '/'))
#colnames(RCTD_norm_weights) = numeric_celltype$assigned_cluster
RCTD_norm_weights_tibb = as_tibble(RCTD_norm_weights, rownames = "barcodes")
RCTD_results_loc = coords %>%
  as_tibble(rownames = "barcodes") %>%
  right_join(RCTD_norm_weights_tibb, by="barcodes")

write_csv(RCTD_results_loc, file.path(RCTDdir, paste0(slide,"_RCTD.csv")))


