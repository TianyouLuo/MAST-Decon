library(tidyverse)
library(viridis)
library(scatterpie)

method = "MASTDecon"

slide_list = c("HCC-2L") #, "HCC-5A", "HCC-5B", "HCC-5C", "HCC-5D"
reffile = "GuilliamsMa"

results_dir = file.path(paste0("/proj/yunligrp/users/tianyou/MASTDecon/HCCDB/results_", reffile), method)
figures_dir = file.path(paste0("/proj/yunligrp/users/tianyou/MASTDecon/HCCDB/figures/sc_", reffile), method)

## Added basophils for Gu2022
## Added pDCs, Cholangiocytes and Hepatocytes for GuilliamsMa
celltypes = c("Malignant cells", "Stromal cells", "Endothelial cells",
              "CD4T cell", "CD8T cell", "Tregs",
              "B cells", "Plasma B cells", "Macrophages",
              "cDC1s", "cDC2s", "pDCs",
              "Basophils", "Cholangiocytes", "Hepatocytes",
              "Monocytes&Monocyte-derived cells", "Neutrophils",
              "Circulating NK_NKT", "Resident NK cell"
              )

for (slide in slide_list){
  result = read_csv(file.path(results_dir, paste0(slide, "_", method, "_medC.csv"))) %>%
    rename(xcoord = x, ycoord = y)
  print(slide)
  for (celltype in celltypes){
    celltype_ano = str_replace_all(celltype, " ", "_")
    celltype_ano = str_replace_all(celltype_ano, "&", "_")
    fig = ggplot(result, aes(x = ycoord, y = 600-xcoord, fill = get(celltype)))+
      geom_point(shape = 21, size = 2) +
      coord_cartesian() +
      scale_fill_viridis(alpha=0.6) +
      labs(fill = "Proportion", x = "", y = "", title = celltype) +
      theme_bw() +
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(), axis.title = element_blank(),
            plot.title = element_text(size = 18),
            legend.text=element_text(size = 18),
            legend.title = element_text(size=16))
    ggsave(file.path(figures_dir, paste0(method, "_", slide, "_",celltype_ano,".png")), fig,
           width = 7.5, height = 6)
  }
  
  pie = ggplot() +
    geom_scatterpie(aes(x = ycoord, y = 600-xcoord), data = result,
                    cols=celltypes, pie_scale = 0.36) +
    coord_equal() +
    labs(x = "", y = "", fill = "Celltypes", title = paste0(method, "(", slide, ")")) +
    theme(legend.position = "bottom",
          axis.line = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(), axis.title = element_blank(),
          plot.title = element_text(size = 28),
          legend.text=element_text(size = 20),
          legend.title = element_text(size=16)) +
    guides(fill = guide_legend(ncol = 3))
  ggsave(file.path(figures_dir, paste0(method,"_",slide,"_pie.png")), pie, width = 12, height = 13.5)
}

