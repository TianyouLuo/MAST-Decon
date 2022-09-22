library(spacexr)
library(Matrix)
library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

#slice = "ST8059048"
#subsample = "5000"
slice = args[1]
subsample = args[2]
stdir = file.path("/proj/yunligrp/users/tianyou/MASTDecon/mouse_brain_cell2location/data/ST/", slice)
scRNAdir = "/proj/yunligrp/users/tianyou/MASTDecon/mouse_brain_cell2location/data/sc/ssp_processed/1m/"
RCTDdir = file.path("/proj/yunligrp/users/tianyou/MASTDecon/mouse_brain_cell2location/RCTD/", slice)

print(RCTDdir)
print(paste0("Subsampling:", subsample))
if (!file.exists(RCTDdir)){
  dir.create(RCTDdir)
}

options(check.name = F) ## to avoid changing numerical gene names
##### read in single cell and ST data ####
ref = fread(file.path(scRNAdir,"sc_cnt.raw.txt"))
cluster = fread(file.path(scRNAdir,"sc_cnt.subclass.tsv"))
if (subsample == "orig"){
  spots <- fread(file.path(stdir,paste0(slice,"_processed.tsv"))) # load in counts matrix
} else {
  spots <- fread(file.path(stdir,paste0(slice,"_sub",subsample,"_processed.tsv"))) # load in counts matrix
}

## convert ST data to the required format: gene on rows and barcodes on columns
spots$gene = as.character(spots$gene)
spots = as.data.frame(spots, check.names = FALSE)
rownames(spots) = spots[,1]
spots[,1] = NULL

# keep only overlap genes with minimum count > 0 
spots = spots[rowSums(spots) > 0,]
overlap_gene = intersect(rownames(spots), colnames(ref)[2:dim(ref)[2]])
ref_overlap = ref[,c("sample_name", overlap_gene), with=FALSE]
spots = spots[overlap_gene,]


##### clean the cluster labels #####
# In the step creating the reference count matrix, 
# we already filtered on celltypes and only kept celltypes with enough cells
# cluster_sel = cluster %>%
#   filter(subclass_label %in% c("L2/3 IT CTX-1", "L4/5 IT CTX", "L5 IT CTX",
#                                "L5 NP CTX", "L5 PT CTX", "L6 CT CTX", "L6 IT CTX",
#                                "L6b CTX", "Lamp5", "Vip", "Sst"))
# 
cell_types <- cluster$subclass_label
names(cell_types) <- cluster$sample_name # create cell_types named list
cell_types[cell_types=="L2/3 IT CTX-1"]="L2-3 IT CTX-1"
cell_types[cell_types=="L4/5 IT CTX"]="L4-5 IT CTX"
cell_types <- as.factor(cell_types) # convert to factor data type

#ref_overlap_sel = ref_overlap[ref_overlap$sample_name %in% cluster_sel$sample_name, ]

##### convert count df into matrix #####
ref_mat = as.data.frame(ref_overlap, check.names = FALSE)
rownames(ref_mat) <- ref_mat[,1]
ref_mat[,1] <- NULL # Move first column to rownames
ref_mat_t = t(ref_mat)

nUMI <- colSums(ref_mat_t)
names(nUMI) <- colnames(ref_mat_t) # create nUMI named list

## create reference data object
reference <- Reference(ref_mat_t, cell_types, nUMI)

##### Create Spatial RNA object #####
coords <- fread(file.path(stdir, paste0(slice, "_tissue_positions.csv")))
coords = coords[coords$V1 %in% colnames(spots), c(1,5,6)]
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
if (subsample == "orig") {
  saveRDS(myRCTD_full, file.path(RCTDdir, "RCTD_full.rds"))
} else {
  saveRDS(myRCTD_full, file.path(RCTDdir, paste0("RCTD_full_","sub",subsample,".rds")))
}

#### Visualize results ####
RCTD_results_full <- myRCTD_full@results
RCTD_norm_weights = as.matrix(sweep(RCTD_results_full$weights, 1, rowSums(RCTD_results_full$weights), '/'))
#colnames(RCTD_norm_weights) = numeric_celltype$assigned_cluster
RCTD_norm_weights_tibb = as_tibble(RCTD_norm_weights, rownames = "barcodes")
RCTD_results_loc = coords %>%
  as_tibble(rownames = "barcodes") %>%
  right_join(RCTD_norm_weights_tibb, by="barcodes") %>%
  select(barcodes, xcoord, ycoord, `L2-3 IT CTX-1`, `L4-5 IT CTX`,
         `L5 IT CTX`, `L5 NP CTX`, `L5 PT CTX`, `L6 CT CTX`, `L6 IT CTX`,
         `L6b CTX`, Lamp5, Pvalb, Sst, Vip, Astro, Oligo)
if (subsample == "orig"){
  write_csv(RCTD_results_loc, file.path(RCTDdir, "RCTD_full_norm.csv"))
} else {
  write_csv(RCTD_results_loc, file.path(RCTDdir, paste0("RCTD_full_norm_", "sub",subsample, ".csv")))
}

library(scatterpie)
cell_type_names = colnames(RCTD_results_loc)[4:dim(RCTD_results_loc)[2]]

STpie = ggplot() +
  geom_scatterpie(aes(x = ycoord, y = 10000 - xcoord), data = RCTD_results_loc,
                  cols=cell_type_names, pie_scale = 0.5) +
  coord_equal()

if (subsample == "orig") {
  ggsave(file.path(RCTDdir, "RCTD_full.png"), STpie, width = 12, height = 8)
} else {
  ggsave(file.path(RCTDdir, paste0("RCTD_full_", "sub", subsample,".png")), STpie, width = 12, height = 8)
}
