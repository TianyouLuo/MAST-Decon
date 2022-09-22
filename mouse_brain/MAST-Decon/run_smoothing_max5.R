library(spacexr)
library(Matrix)
library(tidyverse)
library(optimx)
library(scatterpie)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

slice = args[1]
subsample = args[2]
smoothing = args[3]

RCTDdir = file.path("/proj/yunligrp/users/tianyou/MASTDecon/mouse_brain_cell2location/RCTD/",slice)
resultdir = file.path("/proj/yunligrp/users/tianyou/MASTDecon/mouse_brain_cell2location/MASTDecon/",slice)

print(paste0("RCTD dirctory: ", RCTDdir))
print(paste0("Subsampling: ", subsample))
print(paste0("Smoothing level: ", smoothing))

if (!file.exists(resultdir)){
  dir.create(resultdir)
}

source("/proj/yunligrp/users/tianyou/MASTDecon/mouse_brain_cell2location/MASTDecon/smoothing.R")
## Read in RCTD results as initial starting point
if (subsample == "orig"){
  RCTD = readRDS(file.path(RCTDdir, "RCTD_full.rds"))
} else {
  RCTD = readRDS(file.path(RCTDdir, paste0("RCTD_full_","sub",subsample,".rds")))
}


w = as.matrix(RCTD@results$weights)
UMI = RCTD@spatialRNA@nUMI
spatialRNA = as.matrix(RCTD@spatialRNA@counts)
cell_type_mean_norm = as.matrix(RCTD@cell_type_info$renorm[[1]])

loc = RCTD@spatialRNA@coords
distmat = dist.mat(loc)


#### smoothing parameter ####
if (smoothing == "med") {
  sm = 1
} else if (smoothing == "high") {
  sm = 2
} else if (smoothing == "low") {
  sm = 0.5
}
results = train.smooth(spatialRNA = spatialRNA, theta0 = w, UMI = UMI, 
                       cell_type_mean_norm = cell_type_mean_norm, distmat = distmat,
                       c = sm, radius_seq = c(221, 383, 444), tol = 10^-6, trace = 4,
                       max_CT = 5, threshold = 0.01)


## normalize results
sum_celltype = rowSums(results)
results_norm = as.matrix(sweep(results, 1, sum_celltype, "/"))
results_tibb = as_tibble(results_norm, rownames = "barcodes")
# results_norm_loc = results_norm %>%
#   left_join(as_tibble(loc, rownames="barcodes"), by = "barcodes") %>%
#   rename(xcoord = x, ycoord = y)
results_norm_loc = loc %>%
  as_tibble(rownames = "barcodes") %>%
  rename(xcoord = x, ycoord = y) %>%
  right_join(results_tibb, by="barcodes") %>%
  select(barcodes, xcoord, ycoord, `L2-3 IT CTX-1`, `L4-5 IT CTX`,
         `L5 IT CTX`, `L5 NP CTX`, `L5 PT CTX`, `L6 CT CTX`, `L6 IT CTX`,
         `L6b CTX`, Lamp5, Pvalb, Sst, Vip, Astro, Oligo)

if (subsample == "orig"){
  write_csv(results_norm_loc, 
            file.path(resultdir, 
                      paste0("MASTDecon_",smoothing, "C_max5_norm.csv")))
} else {
  write_csv(results_norm_loc, 
            file.path(resultdir, 
                      paste0("MASTDecon_",smoothing, "C_max5_norm_sub", subsample, ".csv")))
}

## generate pie plots
celltype_names = colnames(results_norm_loc)[4:dim(results_norm_loc)[2]]

STpie = ggplot() +
  geom_scatterpie(aes(x = ycoord, y = 10000-xcoord), data = results_norm_loc,
                  cols=celltype_names, pie_scale = 0.5) +
  coord_equal()

if (subsample == "orig") {
  ggsave(file.path(resultdir, paste0("MASTDecon_",smoothing,"C_max5.png")), 
         STpie, width = 12, height = 8)
} else {
  ggsave(file.path(resultdir, paste0("MASTDecon_",smoothing,"C_max5_sub", subsample,".png")),
         STpie, width = 12, height = 8)
}
