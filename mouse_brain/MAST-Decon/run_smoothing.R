library(RCTD)
library(Matrix)
library(tidyverse)
library(optimx)
library(scatterpie)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

slice = args[1]
subsample = args[2]
smoothing = args[3]

RCTDdir = paste0("/pine/scr/t/i/tianyou/ST/mouse_brain_cell2location/RCTD/",slice,"/sub",subsample)
resultdir = paste0("/pine/scr/t/i/tianyou/ST/mouse_brain_cell2location/MAST-Decon/",slice,"/sub",subsample)

print(paste0("RCTD dirctory: ", RCTDdir))
print(paste0("Smoothing level: ", smoothing))

source("/pine/scr/t/i/tianyou/ST/mouse_brain_cell2location/MAST-Decon/smoothing.R")
RCTD = readRDS(file.path(RCTDdir, "RCTD_full.rds"))

w = as.matrix(RCTD@results$weights)
UMI = RCTD@spatialRNA@nUMI
spatialRNA = as.matrix(RCTD@spatialRNA@counts)
cell_type_mean_norm = as.matrix(RCTD@cell_type_info$renorm[[1]])

loc = RCTD@spatialRNA@coords
distmat = dist.mat(loc)


#### medium smoothing ####
if (smoothing == "med") {
  sm = 1
} else if (smoothing == "high") {
  sm = 2
} else if (smoothing == "low") {
  sm = 0.5
}
results = train.smooth(spatialRNA = spatialRNA, theta0 = w, UMI = UMI, 
                       cell_type_mean_norm = cell_type_mean_norm, distmat = distmat,
                       c = sm, radius_seq = c(221, 382, 442), tol = 10^-6, trace = 0)

results_tibb = as_tibble(results, rownames = "barcodes")
write_csv(results_tibb, file.path(resultdir, paste0("smoothing_", smoothing, "C.csv")))

results_norm = results_tibb
sum_celltype = apply(results_tibb[,-1], 1, sum)
results_norm[,-1] = sweep(results_norm[,-1], 1, sum_celltype, "/")
results_norm_loc = results_norm %>%
  left_join(as_tibble(loc, rownames="barcodes"), by = "barcodes") %>%
  rename(xcoord = x, ycoord = y)
celltype_names = colnames(results_tibb)[-1]

ggplot() +
  geom_scatterpie(aes(x = ycoord, y = 10000-xcoord), data = results_norm_loc,
                  cols=celltype_names, pie_scale = 0.5) +
  coord_equal() +
  ggsave(file.path(resultdir, paste0("smoothing_",smoothing, "C.png")), width = 12, height = 8)
