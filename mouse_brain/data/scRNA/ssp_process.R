library(data.table)

region = fread("metadata.csv")
SSp_meta = region[region$region_label == "SSp",]
SSp_meta = SSp_meta[SSp_meta$subclass_label != ""]
dim(SSp_meta)

label = SSp_meta[,c("sample_name", "subclass_label")]
fwrite(label,"./ssp_processed/sc_cnt.subclass.tsv",col.names=T,row.names=F,
       sep="\t",quote=F)

count = fread("matrix.csv")
SSp_count = count[count$sample_name %in% SSp_meta$sample_name,]
sumUMI = colSums(SSp_count[,2:dim(SSp_count)[2]])
col_to_keep = c("sample_name", names(which(sumUMI != 0)))
SSp_count_nonzero = SSp_count[, ..col_to_keep]
dim(SSp_count_nonzero)

fwrite(SSp_count_nonzero,"./ssp_processed/sc_cnt.raw.txt",col.names=T,row.names=F,
       sep="\t",quote=F)
fwrite(SSp_meta,"./ssp_processed/sc_cnt.meta.tsv",col.names=T,row.names=F,
       sep="\t",quote=F)

