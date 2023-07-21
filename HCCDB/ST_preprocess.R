library(tidyverse)
library(png)
library(ggforce)

gene_dir = "/proj/yunligrp/users/tianyou/data_ST/HCCDB/ST/expression"
out_dir = file.path(gene_dir, "cleaned")
patient = "1"
img_dir = file.path("/proj/yunligrp/users/tianyou/data_ST/HCCDB/ST/histology", paste0("HCC", patient))

dat = readRDS(file.path(gene_dir, paste0("HCC-", patient, "-expr.RDS")))
count = dat@assays$RNA@counts
count = count[!duplicated(rownames(count)),]
meta = dat@meta.data

samp_list = names(table(meta$sample.ident))
for (samp in samp_list){
  print(samp)
  meta_sel = meta %>% filter(sample.ident == samp)
  cnt_sel = count[,rownames(meta_sel)]
  meta_sel = as_tibble(meta_sel, rownames = "barcode")
  cnt_sel = as.matrix(cnt_sel)

  ## gene filtering: remove genes appearing in less than 1% spots
  gene_ncell = rowSums(cnt_sel > 0)
  row_to_keep = names(which(gene_ncell >= dim(cnt_sel)[2]*0.01))
  cnt_nonzero = cnt_sel[row_to_keep,]
  
  write_csv(meta_sel, file.path(out_dir, paste0(samp, "_meta.csv")))
  write_tsv(as_tibble(cnt_nonzero,rownames = "gene"),file.path(out_dir, paste0(samp, "_count.tsv")))
  
  ## overlay clusters onto the histology
  image <- readPNG(file.path(img_dir, paste0(samp, ".png")))
  imgh = unique(meta_sel$height)
  imgw = unique(meta_sel$width)
  
  fig = meta_sel %>%
    mutate(Type = as.factor(Type)) %>%
    ggplot() +
    annotation_raster(image, xmin = 0, xmax = imgw, 
                      ymin = 0, ymax = imgh) +
    geom_point(aes(y = imgh-imagerow, x = imagecol, color = Type), size = 2) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 18),
          legend.position = "bottom")
  
  ggsave(file.path(out_dir, paste0(samp, "_annotation.png")), fig, width = 8, height = 8)
}


