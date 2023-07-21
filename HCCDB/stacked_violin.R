library(tidyverse)
library(viridis)
library(cowplot)

slide_list = c("HCC-2L")
reffile = "GuilliamsMa"

result_dir = file.path(paste0("/proj/yunligrp/users/tianyou/MASTDecon/HCCDB/results_", reffile))
fig_dir = file.path("/proj/yunligrp/users/tianyou/MASTDecon/HCCDB/figures/stacked_violin")

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
annot_dir = "/proj/yunligrp/users/tianyou/data_ST/HCCDB/ST/expression/cleaned/"

for (slide in slide_list){
  MASTDecon = read_csv(file.path(result_dir, "MASTDecon", paste0(slide, "_MASTDecon_medC.csv"))) %>%
    rename("Circulating NK/NKT" = "Circulating NK_NKT")
  RCTD = read_csv(file.path(result_dir, "RCTD", paste0(slide, "_RCTD.csv"))) %>%
    rename("Circulating NK/NKT" = "Circulating NK_NKT")
  truth = read_csv(file.path(annot_dir, paste0(slide, "_meta.csv"))) %>%
    select(barcode, Type) %>%
    rename(barcodes = barcode)
  
  result_m = MASTDecon %>%
    pivot_longer(`B cells`:Tregs, names_to = "celltype", values_to = "MASTDecon") %>%
    left_join(truth, by = "barcodes")
  
  fig_m = ggplot(result_m, aes(factor(Type), MASTDecon, fill = celltype)) +
    geom_violin(scale = "width", adjust = 1, trim = TRUE) +
    scale_y_continuous(expand = c(0, 0), position="right", labels = function(x)
      c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
    facet_grid(rows = vars(celltype), scales = "free", switch = "y") +
    # theme_cowplot(font_size = 12) +
    theme_bw() +
    theme(text=element_text(size=14),
          legend.position = "none", panel.spacing = unit(0, "lines"),
          plot.title = element_text(hjust = 0.5),
          # panel.background = element_rect(fill = NA, color = "black"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          strip.text.y.left = element_text(angle = 0)) +
    xlab("Annotated regions") + ylab("Celltype proportions")
  # ggtitle(method)
  
  ggsave(file.path(fig_dir, paste0(slide, "_", reffile, "_MASTDecon_stacked_violin.png")), fig_m, width = 8, height = 8)
  
  result_r = RCTD %>%
    pivot_longer(`B cells`:Tregs, names_to = "celltype", values_to = "RCTD") %>%
    left_join(truth, by = "barcodes")
  
  fig_r = ggplot(result_r, aes(factor(Type), RCTD, fill = celltype)) +
    geom_violin(scale = "width", adjust = 1, trim = TRUE) +
    scale_y_continuous(expand = c(0, 0), position="right", labels = function(x)
      c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
    facet_grid(rows = vars(celltype), scales = "free", switch = "y") +
    # theme_cowplot(font_size = 12) +
    theme_bw() +
    theme(text=element_text(size=14),
          legend.position = "none", panel.spacing = unit(0, "lines"),
          plot.title = element_text(hjust = 0.5),
          # panel.background = element_rect(fill = NA, color = "black"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          strip.text.y.left = element_text(angle = 0)) +
    xlab("Annotated regions") + ylab("Celltype proportions")
  # ggtitle(method)
  
  ggsave(file.path(fig_dir, paste0(slide, "_", reffile, "_RCTD_stacked_violin.png")), fig_r, width = 8, height = 8)
}
