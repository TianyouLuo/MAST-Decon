library(tidyverse)

gene_dir = "/proj/yunligrp/users/tianyou/data_ST/HCCDB/sc"
out_dir = file.path(gene_dir, "processed")

## For Ma2021 data
source = "Ma2021"

dat = readRDS(file.path(gene_dir, paste0(source, ".download.Rds")))
count = dat@assays$RNA@counts
count = count[!duplicated(rownames(count)),]
meta = dat@meta.data
meta = meta %>%
  filter(!(Celltype %in% c("Basophils", "pDCs"))) %>%
  mutate(Celltype = ifelse(Celltype == "Circulating NK/NKT", "Circulating NK_NKT", Celltype)) %>%
  mutate(Celltype.major = ifelse(Celltype.major == "NK/T cells", "NK_T cells", Celltype.major))
count = count[, rownames(meta)]

cell_UMI = colSums(count)
cell_nGene = colSums(count > 0) ## There is no cells with less than 200 genes
cnt_filtered = count[, cell_UMI >= 500]
dim(cnt_filtered)

gene_ncell = rowSums(cnt_filtered > 0)
row_to_keep = names(which(gene_ncell >= dim(cnt_filtered)[2]*0.01))
cnt_nonzero = cnt_filtered[row_to_keep,]
dim(cnt_nonzero)

write_csv(as_tibble(meta, rownames = "barcode"), file.path(out_dir, paste0(source, "_meta.csv")))
write_tsv(as_tibble(as.matrix(cnt_nonzero),rownames = "gene"),
          file.path(out_dir, paste0(source, "_filtered.txt")))




## For Gu2022 data
source = "Gu2022"

dat = readRDS(file.path(gene_dir, paste0(source, ".download.Rds")))
count = dat@assays$RNA@counts
count = as.matrix(count)
count = count[!duplicated(rownames(count)),]
meta = dat@meta.data
meta = meta %>%
  filter(!(Celltype %in% c("pDCs"))) %>%
  mutate(Celltype = ifelse(Celltype == "Circulating NK/NKT", "Circulating NK_NKT", Celltype)) %>%
  mutate(Celltype.major = ifelse(Celltype.major == "NK/T cells", "NK_T cells", Celltype.major))
rownames(meta) = str_replace(rownames(meta), "data/??????processed_10x/20180416-HCC-P2-Tumor-New/filtered_feature_bc_matrix/_", "")
colnames(count) = str_replace(colnames(count), "data/??????processed_10x/20180416-HCC-P2-Tumor-New/filtered_feature_bc_matrix/_", "")
count = count[, rownames(meta)]

cell_UMI = colSums(count)
cell_nGene = colSums(count > 0) ## There is no cells with less than 200 genes
cnt_filtered = count[, cell_UMI >= 500]
dim(cnt_filtered)

gene_ncell = rowSums(cnt_filtered > 0)
row_to_keep = names(which(gene_ncell >= dim(cnt_filtered)[2]*0.01))
cnt_nonzero = cnt_filtered[row_to_keep,]
dim(cnt_nonzero)

write_csv(as_tibble(meta, rownames = "barcode"), file.path(out_dir, paste0(source, "_meta.csv")))
write_tsv(as_tibble(as.matrix(cnt_nonzero),rownames = "gene"),
          file.path(out_dir, paste0(source, "_filtered.txt")))



## For Guilliams2022 data
source = "Guilliams2022"

dat = readRDS(file.path(gene_dir, paste0(source, ".download.Rds")))
count = dat@assays$RNA@counts
count = as.matrix(count)
count = count[!duplicated(rownames(count)),]
meta = dat@meta.data
meta = meta %>%
  filter(!(Celltype %in% c("Malignant cells"))) %>%
  mutate(Celltype = ifelse(Celltype == "Circulating NK/NKT", "Circulating NK_NKT", Celltype)) %>%
  mutate(Celltype.major = ifelse(Celltype.major == "NK/T cells", "NK_T cells", Celltype.major))
rownames(meta) = str_replace(rownames(meta), "Guilliams2022/", "")
colnames(count) = str_replace(colnames(count), "Guilliams2022/", "")
count = count[, rownames(meta)]

# cell_UMI = colSums(count)
# cell_nGene = colSums(count > 0) ## There is no cells with less than 200 genes
# cnt_filtered = count[, cell_UMI >= 500]
# dim(cnt_filtered)  ## There is no cells needed to be filtered

gene_ncell = rowSums(count > 0)
row_to_keep = names(which(gene_ncell >= dim(count)[2]*0.01))
cnt_nonzero = count[row_to_keep,]
dim(cnt_nonzero)

write_csv(as_tibble(meta, rownames = "barcode"), file.path(out_dir, paste0(source, "_meta.csv")))
write_tsv(as_tibble(as.matrix(cnt_nonzero),rownames = "gene"),
          file.path(out_dir, paste0(source, "_filtered.txt")))

macnt = read_table(file.path(out_dir, "Ma2021_filtered.txt"))
mameta = read_csv(file.path(out_dir, "Ma2021_meta.csv"))

macnt = as.data.frame(macnt, check.names = FALSE)
rownames(macnt) <- macnt[,1]
macnt[,1] <- NULL # Move first column to rownames

overlap_gene = intersect(rownames(macnt), rownames(cnt_nonzero))
macnt_overlap = macnt[overlap_gene,]
cnt_overlap = cnt_nonzero[overlap_gene,]

col_keep = c("barcode", "nCount_RNA", "nFeature_RNA", "Patient_ID",
             "Celltype.major", "dataset", "Tissue", "Phase", "Celltype",
             "UMAP_1", "UMAP_2")
cnt_merge = cbind(cnt_overlap, macnt_overlap)
meta_merge = bind_rows(as_tibble(meta, rownames = "barcode") %>% select(all_of(col_keep)), 
                       mameta %>% select(all_of(col_keep)))
write_csv(meta_merge, file.path(out_dir, "GuilliamsMa_meta.csv"))
write_tsv(as_tibble(as.matrix(cnt_merge),rownames = "gene"),
          file.path(out_dir, "GuilliamsMa_filtered.txt"))

