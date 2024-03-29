library(data.table)
library(tidyverse)
library(aricode)

slice = "ST8059048"

ARI = matrix(NA, nrow = 4, ncol = 5)
rownames(ARI) = c("1000", "1500", "2500", "3500")
colnames(ARI) = c("RCTD", "low-smoothing", "medium-smoothing", "high-smoothing", "Stereoscope")

renormalize = function(weight){
  return(sweep(weight, 1, rowSums(weight), '/'))
}

cal.ARI = function(annotation, weight, n_center = 6, seed = 330){
  weight[,2:12] = renormalize(weight[,2:12])
  result_merge = left_join(weight, annotation, by = "barcodes")
  
  set.seed(seed)
  k.res = kmeans(result_merge[,c(2:12)], centers = n_center, iter.max = 100, nstart = 1)
  result_merge$cluster = k.res$cluster
  
  #adjusted Rand index (ARI), 
  return(ARI(result_merge$cluster, result_merge$layer))
}

for (i in 1:4){
  subsample = rownames(ARI)[i]
  annotation_path = paste0("/pine/scr/t/i/tianyou/ST/mouse_brain_cell2location/data/ST/",slice,"/annotation/",slice,"_markers_layers.csv")
  RCTD_path = paste0("/pine/scr/t/i/tianyou/ST/mouse_brain_cell2location/RCTD/",slice,"/sub",subsample,"/RCTD_full_norm.csv")
  low_path = paste0("/pine/scr/t/i/tianyou/ST/mouse_brain_cell2location/MAST-Decon/",slice,"/sub",subsample,"/smoothing_lowC.csv")
  med_path = paste0("/pine/scr/t/i/tianyou/ST/mouse_brain_cell2location/MAST-Decon/",slice,"/sub",subsample,"/smoothing_medC.csv")
  high_path = paste0("/pine/scr/t/i/tianyou/ST/mouse_brain_cell2location/MAST-Decon/",slice,"/sub",subsample,"/smoothing_highC.csv")
  stereoscope_folder = paste0("/pine/scr/t/i/tianyou/ST/mouse_brain_cell2location/stereoscope/result/",slice,"_sub",subsample,"_cnt.mapped/")
  stereoscope_path = list.files(stereoscope_folder, pattern = "W*.tsv")
  stereoscope_map_path = paste0("/pine/scr/t/i/tianyou/ST/mouse_brain_cell2location/stereoscope/result/",slice,"_mapping.txt")
  
  annotation = fread(annotation_path)
  colnames(annotation) = c("barcodes", "layer")
  RCTD = fread(RCTD_path)[,c(1,4:14)]
  lowsmooth = fread(low_path)
  medsmooth = fread(med_path)
  highsmooth = fread(high_path)
  stereoscope = fread(file.path(stereoscope_folder,stereoscope_path))
  stereoscope_map = fread(stereoscope_map_path)
  stereoscope_merge = left_join(stereoscope, 
                                stereoscope_map %>% select(V1, V7), by = c("V1"="V7")) %>%
    select(V1.y, everything()) %>% select(-V1) %>%
    rename(barcodes = V1.y)
  
  ARI[i,1] = cal.ARI(annotation, RCTD)
  ARI[i,2] = cal.ARI(annotation, lowsmooth)
  ARI[i,3] = cal.ARI(annotation, medsmooth)
  ARI[i,4] = cal.ARI(annotation, highsmooth)
  ARI[i,5] = cal.ARI(annotation, stereoscope_merge)
}

ARI_table = as_tibble(ARI, rownames = "nUMI")
write_csv(ARI_table, "comparison.csv")
