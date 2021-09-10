library(RCTD)
library(Matrix)
library(tidyverse)
library(optimx)
library(scatterpie)

RCTDdir = "/pine/scr/t/i/tianyou/ST/mouse_brain_cell2location/RCTD/ST8059050/"
resultdir = "/pine/scr/t/i/tianyou/ST/mouse_brain_cell2location/MAST-Decon/ST8059050/"

source("/pine/scr/t/i/tianyou/ST/mouse_brain_cell2location/MAST-Decon/smoothing.R")
RCTD = readRDS(file.path(RCTDdir, "RCTD_full.rds"))

w = as.matrix(RCTD@results$weights)
UMI = RCTD@spatialRNA@nUMI
spatialRNA = as.matrix(RCTD@spatialRNA@counts)
cell_type_mean_norm = as.matrix(RCTD@cell_type_info$renorm[[1]])

loc = RCTD@spatialRNA@coords
distmat = dist.mat(loc)


#### medium smoothing #####
results = train.smooth(spatialRNA = spatialRNA, theta0 = w, UMI = UMI, 
                       cell_type_mean_norm = cell_type_mean_norm, distmat = distmat,
                       c = 1, radius_seq = c(221, 382, 442), tol = 10^-6, trace = 0)

results_tibb = as_tibble(results, rownames = "barcodes")
#write_csv(results_tibb, "smoothing_medC.csv")

results_norm = results_tibb
sum_celltype = apply(results_tibb[,-1], 1, sum)
results_norm[,-1] = sweep(results_norm[,-1], 1, sum_celltype, "/")
#write_csv(results_norm_tibb, "smoothing_norm_lowC.csv")
results_norm_loc = results_norm %>%
  left_join(coords, by = "barcodes") %>%
  select(barcodes, xcoord, ycoord, eval(numeric_celltype$assigned_cluster))
write_csv(results_norm_loc, "smoothing_medC_norm.csv")

ggplot() +
  geom_scatterpie(aes(x = ycoord, y = xcoord), data = results_norm_loc, 
                  cols=numeric_celltype$assigned_cluster, pie_scale = 0.8) +
  coord_equal() +
  ggsave("smooth_medC.png")


##### low smoothing ###
results = train.smooth(spatialRNA = spatialRNA, theta0 = w, UMI = UMI, 
                       cell_type_mean_norm = cell_type_mean_norm, distmat = distmat,
                       r0 = 62, c = 0.5, s = 1.42, iter = 3, tol = 10^-6, trace = 0)

results_tibb = as_tibble(results, rownames = "barcodes")

results_norm = results_tibb
sum_celltype = apply(results_tibb[,-1], 1, sum)
results_norm[,-1] = sweep(results_norm[,-1], 1, sum_celltype, "/")
results_norm_loc = results_norm %>%
  left_join(coords, by = "barcodes") %>%
  select(barcodes, xcoord, ycoord, eval(numeric_celltype$assigned_cluster))
write_csv(results_norm_loc, "smoothing_lowC_norm.csv")

ggplot() +
  geom_scatterpie(aes(x = ycoord, y = xcoord), data = results_norm_loc, 
                  cols=numeric_celltype$assigned_cluster, pie_scale = 0.8) +
  coord_equal() +
  ggsave("smooth_lowC.png")


##### high smoothing ###
results = train.smooth(spatialRNA = spatialRNA, theta0 = w, UMI = UMI, 
                       cell_type_mean_norm = cell_type_mean_norm, distmat = distmat,
                       r0 = 62, c = 2, s = 1.42, iter = 3, tol = 10^-6, trace = 0)

results_tibb = as_tibble(results, rownames = "barcodes")

results_norm = results_tibb
sum_celltype = apply(results_tibb[,-1], 1, sum)
results_norm[,-1] = sweep(results_norm[,-1], 1, sum_celltype, "/")
results_norm_loc = results_norm %>%
  left_join(coords, by = "barcodes") %>%
  select(barcodes, xcoord, ycoord, eval(numeric_celltype$assigned_cluster))
write_csv(results_norm_loc, "smoothing_highC_norm.csv")

ggplot() +
  geom_scatterpie(aes(x = ycoord, y = xcoord), data = results_norm_loc, 
                  cols=numeric_celltype$assigned_cluster, pie_scale = 0.8) +
  coord_equal() +
  ggsave("smooth_highC.png")






########### draw original RCTD results ########
RCTD_result = w

sum_celltype = apply(RCTD_result, 1, sum)
RCTD_results_norm = sweep(RCTD_result, 1, sum_celltype, "/")
RCTD_results_norm_tibb = as_tibble(RCTD_results_norm)
RCTD_results_norm_tibb$loc = rownames(RCTD_results_norm)
RCTD_results_norm_tibb = RCTD_results_norm_tibb %>% select(loc, everything())
#write_csv(RCTD_results_norm_tibb, "RCTD_results.csv")

RCTD_results_full = coords %>%
  right_join(RCTD_results_norm_tibb, by=c("barcodes"="loc"))
#cell_type_names = colnames(RCTD_results_full)[5:13]

ggplot() +
  geom_scatterpie(aes(x = ycoord, y = xcoord), data = RCTD_results_full, 
                  cols=celltypes, pie_scale = 0.8) +
  coord_equal() +
  ggsave("RCTD_result.png")

