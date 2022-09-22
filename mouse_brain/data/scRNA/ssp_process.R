library(data.table)

input_dir = "/proj/yunligrp/users/jiawen/data/Yao2021/1m/"
out_dir = "/proj/yunligrp/users/tianyou/MASTDecon/mouse_brain_cell2location/data/sc/ssp_processed/1m/"

region = fread(file.path(input_dir, "metadata.csv"))
SSp_meta = region[region$region_label == "SSp",]
SSp_meta = SSp_meta[SSp_meta$subclass_label != ""]
dim(SSp_meta)
celltype_keep = names(table(SSp_meta$subclass_label)[table(SSp_meta$subclass_label)>=200])  ##40 for 7w
SSp_meta = SSp_meta[SSp_meta$subclass_label %in% celltype_keep]
dim(SSp_meta)

label = SSp_meta[,c("sample_name", "subclass_label")]
fwrite(label, file.path(out_dir, "sc_cnt.subclass.tsv"),col.names=T,row.names=F,
       sep="\t",quote=F)

count = fread(file.path(input_dir, "matrix.csv"))
SSp_count = count[count$sample_name %in% SSp_meta$sample_name,]
sumUMI = colSums(SSp_count[,2:dim(SSp_count)[2]])
col_to_keep = c("sample_name", names(which(sumUMI != 0)))
SSp_count_nonzero = SSp_count[, ..col_to_keep]
dim(SSp_count_nonzero)

cell_UMI = rowSums(SSp_count_nonzero[,2:dim(SSp_count_nonzero)[2]])
cell_nGene = rowSums(SSp_count_nonzero[,2:dim(SSp_count_nonzero)[2]] > 0)
SSp_count_filtered = SSp_count_nonzero[cell_UMI > 500 & cell_nGene > 200, ]
dim(SSp_count_filtered)

fwrite(SSp_count_filtered,file.path(out_dir, "sc_cnt.raw.txt"),
       col.names=T,row.names=F,sep="\t",quote=F)
fwrite(SSp_meta,file.path(out_dir, "sc_cnt.meta.tsv"),
       col.names=T,row.names=F,sep="\t",quote=F)

