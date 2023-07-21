library(spacexr)
library(Matrix)
library(tidyverse)
library(optimx)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

slide = args[1]
smoothing = args[2]
reffile = args[3]

RCTDdir = paste0("/proj/yunligrp/users/tianyou/MASTDecon/HCCDB/results_", reffile, "/RCTD/")
resultdir = paste0("/proj/yunligrp/users/tianyou/MASTDecon/HCCDB/results_", reffile, "/MASTDecon/")

print(paste0("RCTD dirctory: ", RCTDdir))
print(paste0("Slide:", slide))
print(paste0("Smoothing level: ", smoothing))
print(paste0("Reference file: ", reffile))

if (!file.exists(resultdir)){
  dir.create(resultdir)
}

source("/proj/yunligrp/users/tianyou/MASTDecon/HCCDB/codes/smoothing.R")
## Read in RCTD results as initial starting point
RCTD = readRDS(file.path(RCTDdir, paste0(slide, "_RCTD.rds")))


w = as.matrix(RCTD@results$weights)
UMI = RCTD@spatialRNA@nUMI
spatialRNA = as.matrix(RCTD@spatialRNA@counts)
cell_type_mean_norm = as.matrix(RCTD@cell_type_info$renorm[[1]])

loc = RCTD@spatialRNA@coords
distmat = dist.mat(loc)
#sort(unique(as.vector(distmat)))[1:100]

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
                       c = sm, radius_seq = c(7.5, 12.5, 14.5), tol = 10^-6, trace = 0)


## normalize results
sum_celltype = rowSums(results)
results_norm = as.matrix(sweep(results, 1, sum_celltype, "/"))
results_tibb = as_tibble(results_norm, rownames = "barcodes")
# results_norm_loc = results_norm %>%
#   left_join(as_tibble(loc, rownames="barcodes"), by = "barcodes") %>%
#   rename(xcoord = x, ycoord = y)
results_norm_loc = loc %>%
  as_tibble(rownames = "barcodes") %>%
  right_join(results_tibb, by="barcodes")

write_csv(results_norm_loc, file.path(resultdir, paste0(slide, "_MASTDecon_",smoothing, "C.csv")))


