library(data.table)
library(tidyverse)
library(aricode)

result_dir = "/proj/yunligrp/users/tianyou/MASTDecon/HCCDB/results_Gu2022/"
savedir = "/proj/yunligrp/users/tianyou/MASTDecon/HCCDB/comparison"

slide_list = c("HCC-5A", "HCC-5B", "HCC-5C", "HCC-5D", "HCC-2L")

renormalize = function(weight){
  return(sweep(weight, 1, rowSums(weight), '/'))
}

cal.ARI = function(annotation, weight, n_center = 6, seed = 330){
  weight[,2:dim(weight)[2]] = renormalize(weight[,2:dim(weight)[2]])
  result_merge = left_join(weight, annotation, by = "barcodes") %>%
    select(barcodes, layer, everything())
  
  set.seed(seed)
  k.res = kmeans(result_merge[,-c(1:2)], centers = n_center, iter.max = 100, nstart = 1)
  result_merge$cluster = k.res$cluster
  
  #adjusted Rand index (ARI), 
  return(ARI(result_merge$cluster, result_merge$layer))
}

ARI = matrix(NA, nrow = length(slide_list), ncol = 2)
colnames(ARI) = c("RCTD", "MASTDecon")
rownames(ARI) = slide_list

for (i in 1:length(slide_list)){
  slide = slide_list[i]
  print(paste0("Processing slide ", slide))
  
  annotation_path = paste0("/proj/yunligrp/users/tianyou/data_ST/HCCDB/ST/expression/cleaned/",slide,"_meta.csv")
  RCTD = fread(file.path(result_dir, "RCTD", paste0(slide, "_RCTD.csv")))[,-c(2:3)]
  MASTDecon = fread(file.path(result_dir, "MASTDecon", paste0(slide, "_MASTDecon_medC.csv")))[,-c(2:3)]
  annotation = fread(annotation_path) %>%
    select(barcode, Type)
  colnames(annotation) = c("barcodes", "layer")
  n_cluster = length(table(annotation$layer))
  
  if (length(unique((c(dim(RCTD)[1], dim(MASTDecon)[1]))))!=1){
    warning("Different numbers of spots for different methods' results. Please check.")
  }
  
  if (length(unique((c(dim(RCTD)[2], dim(MASTDecon)[2]))))!=1){
    warning("Cell types for different methods' results are not aligned.")
  }
  
  ARI[i,1] = cal.ARI(annotation, RCTD, n_cluster)
  ARI[i,2] = cal.ARI(annotation, MASTDecon, n_cluster)
}

ARI_table = as_tibble(ARI, rownames = "slide")
write_csv(ARI_table, file.path(savedir, "ARI_comparison.csv"))

