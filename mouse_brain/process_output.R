library(tidyverse)
library(data.table)
library(scatterpie)

RCTD_dir = "/proj/yunligrp/users/tianyou/MASTDecon/mouse_brain_cell2location/RCTD"
RCTD_savedir = "/proj/yunligrp/users/tianyou/MASTDecon/mouse_brain_cell2location/RCTD_cleaned"
card_dir = "/proj/yunligrp/users/jiawen/st/tianyou/card/"
card_savedir = "/proj/yunligrp/users/tianyou/MASTDecon/mouse_brain_cell2location/CARD"
cell2location_dir = "/proj/yunligrp/users/jiawen/st/tianyou/cell2location/result"
cell2location_savedir = "/proj/yunligrp/users/tianyou/MASTDecon/mouse_brain_cell2location/cell2location"
stereoscope_dir = "/proj/yunligrp/users/jiawen/st/tianyou/stereoscope"
stereoscope_savedir = "/proj/yunligrp/users/tianyou/MASTDecon/mouse_brain_cell2location/stereoscope"
figs_dir = "/proj/yunligrp/users/tianyou/MASTDecon/mouse_brain_cell2location/figures/piechart"

slide_list = c("ST8059048", "ST8059049", "ST8059050", "ST8059051", "ST8059052")
sub_list = c("1000", "2000", "3000", "4000", "5000", "7500", "10000", "full")
cols = c("L2/3 IT CTX-1", "L4/5 IT CTX", "L5 IT CTX", "L5 NP CTX", "L5 PT CTX", 
         "L6 CT CTX", "L6 IT CTX", "L6b CTX", "Lamp5", "Pvalb", "Sst", 
         "Vip", "Astro", "Oligo")

keep.maxCT = function(x, keep = 5, threshold = 0.01){
  x[rank(-x, ties.method = "random") > keep] = 0
  x = x/sum(x)
  x[x <= threshold] = 0
  x = x/sum(x)
  return(x)
}

for (s in 1:length(slide_list)){
  slide = slide_list[s]
  if (!dir.exists(file.path(RCTD_savedir, slide))){
    dir.create(file.path(RCTD_savedir, slide))
  }
  if (!dir.exists(file.path(card_savedir, slide))){
    dir.create(file.path(card_savedir, slide))
  }
  if (!dir.exists(file.path(cell2location_savedir, slide))){
    dir.create(file.path(cell2location_savedir, slide))
  }
  if (!dir.exists(file.path(stereoscope_savedir, slide))){
    dir.create(file.path(stereoscope_savedir, slide))
  }
  for (i in 1:length(sub_list)){
    subsample = sub_list[i]
    if (subsample == "full") {
      RCTD = fread(file.path(RCTD_dir, slide, "RCTD_full_norm.csv"))
      card = fread(file.path(card_dir, paste0(slide, "_full.csv")))
      c2l = fread(file.path(cell2location_dir, paste0(slide, "_full.csv")))
      stereo_f = list.files(file.path(stereoscope_dir, slide, "full",
                                      paste0(slide, "_processed")),
                                      "W*.tsv", full.names = T)
      stereoscope = fread(stereo_f)
    } else {
      RCTD = fread(file.path(RCTD_dir, slide, 
                             paste0("RCTD_full_norm_sub", subsample, ".csv")))
      card = fread(file.path(card_dir, paste0(slide, "_sub", subsample, ".csv")))
      c2l = fread(file.path(cell2location_dir, paste0(slide, "_sub", subsample, ".csv")))
      stereo_f = list.files(file.path(stereoscope_dir, slide, paste0("sub", subsample),
                                      paste0(slide, "_sub", subsample,"_processed")),
                            "W*.tsv", full.names = T)
      stereoscope = fread(stereo_f)
    }
    
    spot_info = RCTD[,1:3]
    ## process RCTD results
    RCTD = RCTD %>%
      rename(`L2/3 IT CTX-1` = `L2-3 IT CTX-1`,
             `L4/5 IT CTX` = `L4-5 IT CTX`)
    RCTD_new = data.table(t(apply(RCTD %>% select(all_of(cols)), 1, 
                                  keep.maxCT,keep = 5,threshold = 0.01)))
    RCTD_new$barcodes = RCTD$barcodes
    RCTD_new$xcoord = RCTD$xcoord
    RCTD_new$ycoord = RCTD$ycoord
    RCTD_new = RCTD_new %>% select(barcodes, xcoord, ycoord, everything())
    RCTDpie = ggplot() +
      geom_scatterpie(aes(x = ycoord, y = 10000-xcoord), data = RCTD_new,
                      cols=cols, pie_scale = 0.5) +
      coord_equal()
    
    ## process CARD results
    card = card %>% rename(barcodes = V1) %>%
      mutate(barcodes = str_replace(barcodes, "\\.", "-"))
    card_new = data.table(t(apply(card %>% select(all_of(cols)), 1, 
                                  keep.maxCT,keep = 5,threshold = 0.01)))
    card_new$barcodes = card$barcodes
    card_new = left_join(spot_info, card_new, by = "barcodes")
    cardpie = ggplot() +
      geom_scatterpie(aes(x = ycoord, y = 10000-xcoord), data = card_new,
                      cols=cols, pie_scale = 0.5) +
      coord_equal()
    
    ## process cell2location results
    c2l = c2l %>% rename(barcodes = V1)
    c2l_new = data.table(t(apply(c2l %>% select(all_of(cols)), 1, 
                                 keep.maxCT,keep = 5,threshold = 0.01)))
    c2l_new$barcodes = c2l$barcodes
    c2l_new = left_join(spot_info, c2l_new, by = "barcodes")
    c2lpie = ggplot() +
      geom_scatterpie(aes(x = ycoord, y = 10000-xcoord), data = c2l_new,
                      cols=cols, pie_scale = 0.5) +
      coord_equal()
    
    ## process stereoscope results
    stereoscope = stereoscope %>% rename(barcodes = V1)
    stereoscope_new = data.table(t(apply(stereoscope %>% select(all_of(cols)), 1, 
                                         keep.maxCT,keep = 5,threshold = 0.01)))
    stereoscope_new$barcodes = stereoscope$barcodes
    stereoscope_new = left_join(spot_info, stereoscope_new, by = "barcodes")
    stereopie = ggplot() +
      geom_scatterpie(aes(x = ycoord, y = 10000-xcoord), data = stereoscope_new,
                      cols=cols, pie_scale = 0.5) +
      coord_equal()
    
    if (subsample == "full"){
      write_csv(RCTD_new, file.path(RCTD_savedir, slide, "RCTD_max5_norm.csv"))
      write_csv(card_new, file.path(card_savedir, slide, 
                                    paste0("CARD_max5_norm.csv")))
      write_csv(c2l_new, file.path(cell2location_savedir, slide, 
                                   paste0("cell2location_max5_norm.csv")))
      write_csv(stereoscope_new, file.path(stereoscope_savedir, slide, 
                                           paste0("stereoscope_max5_norm.csv")))
      ggsave(file.path(figs_dir, paste0(slide, "_RCTD_max5_norm.png")), 
             RCTDpie, width = 12, height = 8)
      ggsave(file.path(figs_dir, paste0(slide, "_CARD_max5_norm.png")), 
             cardpie, width = 12, height = 8)
      ggsave(file.path(figs_dir, paste0(slide, "_cell2location_max5_norm.png")), 
             c2lpie, width = 12, height = 8)
      ggsave(file.path(figs_dir, paste0(slide, "_stereoscope_max5_norm.png")), 
             stereopie, width = 12, height = 8)
    } else {
      write_csv(RCTD_new, file.path(RCTD_savedir, slide, 
                                    paste0("RCTD_max5_norm_sub", subsample,".csv")))
      write_csv(card_new, file.path(card_savedir, slide, 
                                    paste0("CARD_max5_norm_sub", subsample,".csv")))
      write_csv(c2l_new, file.path(cell2location_savedir, slide, 
                                   paste0("cell2location_max5_norm_sub", subsample,".csv")))
      write_csv(stereoscope_new, file.path(stereoscope_savedir, slide, 
                                           paste0("stereoscope_max5_norm_sub", subsample,".csv")))
      ggsave(file.path(figs_dir, paste0(slide, "_RCTD_max5_norm_sub",
                                        subsample, ".png")), 
             RCTDpie, width = 12, height = 8)
      ggsave(file.path(figs_dir, paste0(slide, "_CARD_max5_norm_sub",
                                        subsample, ".png")), 
             cardpie, width = 12, height = 8)
      ggsave(file.path(figs_dir, paste0(slide, "_cell2location_max5_norm_sub",
                                        subsample, ".png")), 
             c2lpie, width = 12, height = 8)
      ggsave(file.path(figs_dir, paste0(slide, "_stereoscope_max5_norm_sub",
                                        subsample, ".png")), 
             stereopie, width = 12, height = 8)
    }
  }
}
