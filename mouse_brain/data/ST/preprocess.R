library(data.table)

STdir = "/pine/scr/t/i/tianyou/ST/mouse_brain_cell2location/data/ST/ST8059051/"
count = fread(file.path(STdir,"raw_feature_bc_matrix/matrix.mtx.gz"),header=F,fill=TRUE,skip=2)
cnt = matrix(rep(0,count$V1[1]*count$V2[1]),nrow=count$V2[1])
for(i in 2:nrow(count)){
  cnt[count$V2[i],count$V1[i]]=count$V3[i]
}
gene = fread(file.path(STdir, "raw_feature_bc_matrix/features.tsv.gz"),header=F)
barcode = fread(file.path(STdir, "raw_feature_bc_matrix/barcodes.tsv.gz"),header=F)
annotation = fread(file.path(STdir, "annotation/ST8059051_markers_layers.csv"))
annotation = annotation[annotation$ST8059051 != ""]
duplicate_gene = gene$V2[which(duplicated(gene$V2))]
sum(duplicated(gene$V2))

colnames(cnt) = gene$V1
rownames(cnt) = barcode$V1
cnt = cnt[annotation$Barcode,]

# for genes with more than one Entrez id, we do not include them here
gene = gene[!(gene$V2 %in% duplicate_gene),]
cnt = cnt[,gene$V1]
sum(colnames(cnt)==gene$V1)
colnames(cnt) = gene$V2

gene_sum = apply(cnt, 2, sum)
cnt = cnt[,gene_sum!=0]

gene_present = apply(cnt, 2, function(x){sum(x>0)})
cnt = cnt[,gene_present > 0.02 * dim(cnt)[1]]
cnt_t = t(cnt)

write.table(cnt_t, file.path(STdir, "ST8059051_processed.tsv"),row.names=T,col.names=T,sep="\t",quote=F)
