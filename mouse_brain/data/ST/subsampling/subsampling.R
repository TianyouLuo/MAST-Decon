library(data.table)

slice = "ST8059051"
STdir = paste0("/pine/scr/t/i/tianyou/ST/mouse_brain_cell2location/data/ST/",slice)
ST_orig = fread(file.path(STdir, paste0(slice, "_processed.tsv")))

## function to subsample to a read depth for one spot
subsample.spot = function(ST, depth){
  prob_vec = ST / sum(ST)
  reads = rmultinom(1, depth, prob_vec)
  return(reads)
}

subsample.ST = function(STdat, depth = 100000, sd = 5000){
  subsampled_mat = STdat
  depth_vec = round(rnorm(dim(STdat)[2]-1, depth, sd))
  for (i in 2:dim(STdat)[2]){
    vec = STdat[,..i][[1]] ##..used for data.table
    subsample_vec = subsample.spot(vec, depth_vec[i-1])
    subsampled_mat[,i] = subsample_vec
  }
  return(subsampled_mat)
}

process.ST = function(cnt){
  gene_sum = apply(cnt[,2:dim(cnt)[2]], 1, sum)
  cnt = cnt[gene_sum!=0,]
  
  gene_present = apply(cnt[,2:dim(cnt)[2]], 1, function(x){sum(x>0)})
  cnt = cnt[gene_present > 0.02 * (dim(cnt)[2]-1),]
}

set.seed(2021)
ST_sub2500 = subsample.ST(ST_orig, depth = 2500, sd = 500)
ST_sub2500_process = process.ST(ST_sub2500)
write.table(ST_sub2500_process, file.path(STdir, paste0(slice, "_sub2500_processed.tsv")),
            row.names=F,col.names=T,sep="\t",quote=F)

set.seed(9050)
ST_sub3500 = subsample.ST(ST_orig, depth = 3500, sd = 500)
ST_sub3500_process = process.ST(ST_sub3500)
write.table(ST_sub3500_process, file.path(STdir, paste0(slice, "_sub3500_processed.tsv")),
            row.names=F,col.names=T,sep="\t",quote=F)

set.seed(805)
ST_sub1500 = subsample.ST(ST_orig, depth = 1500, sd = 200)
ST_sub1500_process = process.ST(ST_sub1500)
write.table(ST_sub1500_process, file.path(STdir, paste0(slice, "_sub1500_processed.tsv")),
            row.names=F,col.names=T,sep="\t",quote=F)


set.seed(1000)
ST_sub1000 = subsample.ST(ST_orig, depth = 1000, sd = 120)
ST_sub1000_process = process.ST(ST_sub1000)
write.table(ST_sub1000_process, file.path(STdir, paste0(slice, "_sub1000_processed.tsv")),
            row.names=F,col.names=T,sep="\t",quote=F)

set.seed(500)
ST_sub5000 = subsample.ST(ST_orig, depth = 5000, sd = 500)
ST_sub5000_process = process.ST(ST_sub5000)
write.table(ST_sub5000_process, file.path(STdir, paste0(slice, "_sub5000_processed.tsv")),
            row.names=F,col.names=T,sep="\t",quote=F)

set.seed(751)
ST_sub7500 = subsample.ST(ST_orig, depth = 7500, sd = 500)
ST_sub7500_process = process.ST(ST_sub7500)
write.table(ST_sub7500_process, file.path(STdir, paste0(slice, "_sub7500_processed.tsv")),
            row.names=F,col.names=T,sep="\t",quote=F)

set.seed(104)
ST_sub10000 = subsample.ST(ST_orig, depth = 10000, sd = 1000)
ST_sub10000_process = process.ST(ST_sub10000)
write.table(ST_sub10000_process, file.path(STdir, paste0(slice, "_sub10000_processed.tsv")),
            row.names=F,col.names=T,sep="\t",quote=F)
