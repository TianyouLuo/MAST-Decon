library(tidyverse)
library(viridis)
library(scatterpie)
library(data.table)

slide = "HCC-5C"
gene_list = c("CD4", "GATA3", "CD8A", "CD8B", "CCR7", "CD14", "CD68", "MSR1", "CD11b", "ITGAM", "CCR2")

data_dir = "/proj/yunligrp/users/tianyou/data_ST/HCCDB/ST/expression/cleaned"
figures_dir = "/proj/yunligrp/users/tianyou/MASTDecon/HCCDB/figures/gene_map"

dat = fread(file.path(data_dir, paste0(slide, "_count.tsv")))
dat = as.data.frame(dat, check.names = FALSE)
rownames(dat) = dat[,1]
dat[,1] = NULL
meta = read_csv(file.path(data_dir, paste0(slide, "_meta.csv"))) %>%
  select(barcode, imagerow, imagecol) %>%
  rename(xcoord = imagerow, ycoord = imagecol)

print(slide)
for (gene in gene_list){
  dat_tmp = dat[rownames(dat) == gene,]
  if (dim(dat_tmp)[1] == 0){
    warning(paste0("Provided gene ", gene, " is not present in the ST data."))
  } else{
    dat_tmp_t = as_tibble(t(dat_tmp), rownames = "barcode") %>%
      left_join(meta, by = "barcode")
    fig = ggplot(dat_tmp_t, aes(x = ycoord, y = 600-xcoord, fill = get(gene)))+
      geom_point(shape = 21, size = 2) +
      coord_cartesian() +
      scale_fill_viridis(alpha=0.6) +
      labs(fill = "", x = "", y = "", title = gene) +
      theme_bw() +
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(), axis.title = element_blank(),
            plot.title = element_text(size = 18),
            legend.text=element_text(size = 18),
            legend.title = element_text(size=16))
    ggsave(file.path(figures_dir, paste0(slide, "_",gene,".png")), fig,
           width = 7.5, height = 6)
  }
}


