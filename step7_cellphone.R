# conda activate cellphonedbpy36
# cd ~/AM/cellphonedb
# write.table(as.matrix(seurat_obj@assays$RNA@data), 'cellphonedb_count.txt', sep='\t', quote=F)
# meta_data <- cbind(rownames(seurat_obj@meta.data), pred.hesc$pruned.labels)  
# meta_data <- as.matrix(meta_data)
# meta_data[is.na(meta_data)] = "Unkown" #  细胞类型中不能有NA
# 
# write.table(meta_data, 'cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)
# cellphonedb method statistical_analysis  cellphonedb_meta.txt  cellphonedb_count.txt      --counts-data=gene_name
# #如果我们count的基因是基因名格式，需要添加参数--counts-data=gene_name，如果行名为ensemble名称的话，可以不添加这个参数，使用默认值即可。
