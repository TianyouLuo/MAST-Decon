#devtools::install_github("dmcable/RCTD", build_vignettes = TRUE)
#browseVignettes('RCTD')

library(RCTD)
library(Matrix)
library(tidyverse)
library(data.table)
#library(jcolors)

slice = "ST8059050"
subsample = "3500"
stdir = paste0("/pine/scr/t/i/tianyou/ST/mouse_brain_cell2location/data/ST/", slice)
scRNAdir = "/pine/scr/t/i/tianyou/ST/mouse_brain_cell2location/data/scRNA/ssp_processed/"
RCTDdir = paste0("/pine/scr/t/i/tianyou/ST/mouse_brain_cell2location/RCTD/", slice)

##### read in single cell and ST data ####
ref = fread(file.path(scRNAdir,"sc_cnt.raw.txt"))
cluster = fread(file.path(scRNAdir,"sc_cnt.subclass.tsv"))
if (is.na(subsample)){
  spots <- fread(file.path(stdir,paste0(slice,"_processed.tsv"))) # load in counts matrix
} else {
  spots <- fread(file.path(stdir,paste0(slice,"_sub",subsample,"_processed.tsv"))) # load in counts matrix
}
spots$gene = as.character(spots$gene)
overlap_gene = intersect(spots$gene, colnames(ref)[2:dim(ref)[2]])
ref_overlap = ref[,c("sample_name", overlap_gene), with=FALSE]
spots = spots[spots$gene %in% overlap_gene]

## clean the cluster labels to keep only clusters with enough cells ##
cluster_sel = cluster %>%
  filter(subclass_label %in% c("L2/3 IT CTX-1", "L4/5 IT CTX", "L5 IT CTX",
                               "L5 NP CTX", "L5 PT CTX", "L6 CT CTX", "L6 IT CTX",
                               "L6b CTX", "Lamp5", "Vip", "Sst"))

cell_types <- cluster_sel$subclass_label
names(cell_types) <- cluster_sel$sample_name # create cell_types named list
cell_types[cell_types=="L2/3 IT CTX-1"]="L2-3 IT CTX-1"
cell_types[cell_types=="L4/5 IT CTX"]="L4-5 IT CTX"
cell_types <- as.factor(cell_types) # convert to factor data type

ref_overlap_sel = ref_overlap[ref_overlap$sample_name %in% cluster_sel$sample_name, ]

nUMI <- rowSums(ref_overlap_sel[,2:dim(ref_overlap_sel)[2]])
names(nUMI) <- ref_overlap_sel$sample_name # create nUMI named list

## convert count df into matrix ###
ref_mat = as.data.frame(ref_overlap_sel)
rownames(ref_mat) <- ref_mat[,1]
ref_mat[,1] <- NULL # Move first column to rownames
ref_mat_t = t(ref_mat)


## create reference data object
reference <- Reference(ref_mat_t, cell_types, nUMI)
#saveRDS(reference, file.path(RCTDdir,'SCRef.rds'))

##reference=readRDS(file.path(RCTDdir,'SCRef.rds'))

#### Spatial RNA files ####
spots_mat = as.data.frame(spots)
coords <- fread(file.path(stdir,"/spatial/tissue_positions_list.csv"))
rownames(spots_mat) <- spots_mat[,1]
spots_mat[,1] <- NULL # Move first column to rownames

coords = coords[coords$V1 %in% colnames(spots), c(1,5,6)]
coords = as.data.frame(coords)
colnames(coords) = c("barcodes", "xcoord", "ycoord")
rownames(coords) <- coords$barcodes
coords$barcodes <- NULL # Move barcodes to rownames
nUMI_spot <- colSums(spots_mat) # In this case, total counts per pixel is nUMI

### Create SpatialRNA object
puck <- SpatialRNA(coords, spots_mat, nUMI_spot)


## create RCTD object
myRCTD <- create.RCTD(puck, reference, max_cores = 1,
                      gene_cutoff = 1e-4, fc_cutoff = 0.5,
                      gene_cutoff_reg = 2e-4, fc_cutoff_reg = 0.7)

#### run RCTD ####
myRCTD@config$N_epoch = 20 ## the number of iterations that choose_sigma_c takes
myRCTD_full <- run.RCTD(myRCTD, doublet_mode = 'full')
if (is.na(subsample)) {
  saveRDS(myRCTD_full, file.path(RCTDdir, "RCTD_full.rds"))
} else {
  saveRDS(myRCTD_full, file.path(RCTDdir, paste0("sub",subsample),"RCTD_full.rds"))
}

#### Visualize results ####
RCTD_results_full <- myRCTD_full@results
RCTD_norm_weights = as.matrix(sweep(RCTD_results_full$weights, 1, rowSums(RCTD_results_full$weights), '/'))
#colnames(RCTD_norm_weights) = numeric_celltype$assigned_cluster
RCTD_norm_weights_tibb = as_tibble(RCTD_norm_weights, rownames = "barcodes")
RCTD_results_loc = coords %>%
  as_tibble(rownames = "barcodes") %>%
  right_join(RCTD_norm_weights_tibb, by="barcodes")
if (is.na(subsample)){
  write_csv(RCTD_results_loc, file.path(RCTDdir, "RCTD_full_norm.csv"))
} else {
  write_csv(RCTD_results_loc, file.path(RCTDdir, paste0("sub",subsample), "RCTD_full_norm.csv"))
}

library(scatterpie)
cell_type_names = colnames(RCTD_results_loc)[4:dim(RCTD_results_loc)[2]]

STpie = ggplot() +
 geom_scatterpie(aes(x = ycoord, y = 10000 - xcoord), data = RCTD_results_loc,
                 cols=cell_type_names, pie_scale = 0.5) +
 coord_equal()

if (is.na(subsample)) {
  ggsave(file.path(RCTDdir, "RCTD_full.png"), STpie, width = 12, height = 8)
} else {
  ggsave(file.path(RCTDdir, paste0("sub", subsample),"RCTD_full.png"), STpie, width = 12, height = 8)
}
