library(spacexr)
library(Matrix)
library(tidyverse)
library(optimx)
library(scatterpie)
library(data.table)


dir = "/proj/yunligrp/users/tianyou/MASTDecon/mouse_brain_cell2location/RCTD/"
slice = c("ST8059048", "ST8059049", "ST8059050", "ST8059051", "ST8059052")
print(paste0("RCTD dirctory: ", dir))

sc_dir = "/proj/yunligrp/users/tianyou/MASTDecon/mouse_brain_cell2location/data/sc/ssp_processed/1m"
sc = fread(file.path(sc_dir, "sc_cnt.raw.txt"))
outdir = "/proj/yunligrp/users/tianyou/MASTDecon/mouse_brain_cell2location/data/sc/ssp_processed/1m/RCTD_markers"

keep_lst = list()
sub_vec = c("1000", "2000", "3000", "4000", "5000", "7500", "10000", "orig")
for (sub in 1:8){
  subsample = sub_vec[sub]
  print(paste0("Subsampling: ", subsample))
  gene_list = list()
  for (i in 1:length(slice)){
    RCTDdir = file.path(dir, slice[i])
    if (subsample == "orig"){
      RCTD = readRDS(file.path(RCTDdir, "RCTD_full.rds"))
    } else {
      RCTD = readRDS(file.path(RCTDdir, paste0("RCTD_full_","sub",subsample,".rds")))
    }
    gene_list[[slice[i]]] = RCTD@internal_vars$gene_list_bulk
  }
  keep_lst[[subsample]] = Reduce(union, gene_list)
  col_to_keep = c("sample_name", Reduce(union, gene_list))
  sc_out = sc[, ..col_to_keep]
  if (subsample == "orig"){
    fwrite(sc_out,file.path(outdir, "sc_cnt_marker_full.txt"),
           col.names=T,row.names=F,sep="\t",quote=F)
  } else {
    fwrite(sc_out,file.path(outdir, paste0("sc_cnt_marker_sub", subsample, ".txt")),
           col.names=T,row.names=F,sep="\t",quote=F)
  }
}

