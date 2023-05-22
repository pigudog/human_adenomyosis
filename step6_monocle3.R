# BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
#                        'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
#                        'SummarizedExperiment', 'batchelor', 'HDF5Array',
#                        'terra', 'ggrastr'),force = TRUE)
# # speedglm包违反CRAN规定被删除了，所以要从rstudio下载
# install.packages("pak")
# install.packages("pak", repos = "https://r-lib.github.io/p/pak/dev/")
# pkgs <- c("hhoeflin/hdf5r","mojaveazure/loomR@develop",
#           "url::https://packagemanager.rstudio.com/cran/2023-03-31/src/contrib/speedglm_0.3-4.tar.gz")
# pak::pkg_install(pkgs)
# devtools::install_github('cole-trapnell-lab/leidenbase')
# devtools::install_github('cole-trapnell-lab/monocle3')



##1、重新创建seurat对象
rm(list=ls())
library(Seurat)
library(tidyverse)
library(Matrix)
library( magrittr)
library(sp)
library(monocle3)
load(file = 'AM_after_cluster.Rdata')

ps <- function(filename, plot = FALSE, w = 9, h = 6) {
  if (is.object(plot)) {
    print(plot)
  }
  plot <- recordPlot()
  pdf(file = paste0(filename,".pdf"), onefile = T, width = w, height = h)
  replayPlot(plot)
  dev.off()
}
n_jobs <- 48 # 使用的线程数

##Read1
AM@meta.data$cell_type <- AM@active.ident
AM_EC <- AM[,AM$orig.ident %in% "EC"]
# AM1, unkown1 to (Mac and Mast)
AM1 <- AM_EC[,AM_EC$cell_type %in% c("Mast","Mac","Unkown1")]
pbmc <-AM1
# AM2, 第二个monocle：fibro-like，unkonwn2，SMA
AM2 <- AM_EC[,AM_EC$cell_type %in% c("VSMC", "Pericyte","SMC","ESC")]
pbmc <-AM2
##2.拆分成3个文件
# data <- GetAssayData(scRNAsub, assay = 'RNA', slot = 'counts')
# cell_metadata <- scRNAsub@meta.data
# gene_annotation <- data.frame(gene_short_name = rownames(data))
# rownames(gene_annotation) <- rownames(data)

expression_matrix<-as(as.matrix(pbmc@assays$RNA@counts),"sparseMatrix")
cell_metadata<-pbmc@meta.data
gene_annotation<-data.frame(gene_short_name=row.names(expression_matrix), row.names=row.names(expression_matrix))


##本步为数据的标准化，以及预降维，为之后的降维做准备
#preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

hvg <- pbmc@assays$SCT@var.features

pbmc_embed <- Seurat::Embeddings(pbmc, reduction = "umap")

# rm(pbmc, cell_metadata)
#通过降维，可以再二维图中观察细胞之间的关系
# cds <- reduce_dimension(cds,preprocess_method="UMAP")
## Normalize and pre-process the data
cds = preprocess_cds(cds,norm_method="log",method="PCA",num_dim=100)
# cds %<>% preprocess_cds(num_dim = 30,method="PCA",norm_method="log",use_genes=hvg)
## Remove batch effects with cell alignment
# cds %<>% align_cds(alignment_group = "orig.ident",preprocess_method="PCA")
plot_pc_variance_explained(cds)
ps("./AM/PCA", w = 9, h = 6)
## Reduce the dimensions using UMAP
cds %<>% reduce_dimension(umap.fast_sgd = TRUE, cores = n_jobs, reduction_method = "UMAP", preprocess_method = "PCA")
## 从seurat导入整合过的umap坐标
cds_embed <- cds@int_colData$reducedDims$UMAP
cds@int_colData$reducedDims$UMAP <- pbmc_embed[rownames(cds_embed),]

plot_cells(cds, reduction_method = "UMAP", 
           color_cells_by = "cell_type",
           group_label_size=4)
ps("./AM/UMAP_CellType")








#为了找出哪些细胞在一个发育轨迹上
# 这里的cluster其实是做分区，不同分区的细胞会进行单独的轨迹分析
cds <- cluster_cells(cds,preprocess_method = "UMAP",cluster_method = "leiden",random_seed=1314)
plot_cells(cds,  reduction_method="UMAP",group_label_size=3)
ps("./AM/EC_AM1UMAP_monocle")
#在上一步分群的基础上，找出每个群内的细胞发育轨迹（无监督，没有基于生物学背景知识）
cds <- learn_graph(cds)
plot_cells(cds,trajectory_graph_segment_size = 1.5,
           group_label_size = 6,group_cells_by="cluster",
           color_cells_by = "cell_type", label_groups_by_cluster = FALSE,
           label_leaves = FALSE, label_branch_points = FALSE
)
ps("./AM/EC_AM1_UMAP_monocle_learn_graph_auto")
# 上面这个图将被用于许多下游分析，比如分支分析和差异表达分析
plot_cells(cds, color_cells_by = "cell_type", label_groups_by_cluster=FALSE,
           label_leaves=TRUE, label_branch_points=TRUE,graph_label_size=1.5)
#黑色的线显示的是graph的结构。
# 数字带白色圆圈表示不同的结局，也就是叶子。
# 数字带黑色圆圈代表分叉点，从这个点开始，细胞可以有多个结局。
# 这些数字可以通过label_leaves和label_branch_points参数设置。
save(cds,file = "monocle2_AM2")
load("monocle2_AM2")
cds <- order_cells(cds)
plot_cells(cds,
           group_label_size = 6,
           color_cells_by = "pseudotime",
           label_cell_groups = TRUE,
           label_leaves = TRUE,
           label_branch_points = TRUE,
           graph_label_size = 0
)
ps("./AM/UMAP_monocle_Trajectory_Pseudotime.pdf")
# 作者将Cluster1细胞、上皮细胞和内皮细胞提取出来做拟时间分析，
# 发现Cluster1细胞（亚群1）呈现由上皮细胞（亚群7、17）特征
# 向内皮细胞特征（亚群2、9、10、15）分化。
# 在拟时间轴上，marker基因也呈现出从上皮细胞（EPCAM、CDH1、KRT7）
# 向内皮细胞（PECAM1、VWF、CDH5）的转化。
# 作者进一步发现该轨迹的末端出现了血管生成拟态（VM）相关的基因（F2R、FLT1、KDR）上调表达。

dea_res <- graph_test(cds, neighbor_graph = "principal_graph", cores = n_jobs) %>%
  dplyr::filter(q_value < 0.05)
write_csv(dea_res, "./AM/EC_AM2_Trajectory_genes.csv")
degs <- dea_res$gene_short_name

top_genes <- dea_res %>%
  dplyr::top_n(n = 10, morans_I) %>%
  dplyr::pull(gene_short_name) %>%
  as.character()
plot_genes_in_pseudotime(cds[top_genes, ],
                         color_cells_by = "cell_type",
                         min_expr = 0.5, ncol = 2
)
ps("./AM/EC_AM2_Top_Genes_Jitterplot.pdf", w = 9, h = 6)

p <- plot_cells(cds,
                genes = top_genes, show_trajectory_graph = FALSE,
                label_cell_groups = FALSE, label_leaves = FALSE
)
p$facet$params$ncol <- 5
ps("./AM/Top_Genes_Featureplot.pdf", plot = p)


#查看特定基因的表达情况
ciliated_genes <- c("EPCAM","PECAM1","F2R","CDH1","VWF","FLT1","KRT7","CDH5","KDR")
plot_cells(cds,
           genes=ciliated_genes,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)
ps("./AM/ciliated_genes_cell_plot.pdf", w = 9, h = 6)

plot_genes_in_pseudotime(cds[ciliated_genes, ],
                         color_cells_by = "cell_type",
                         min_expr = 0.5, ncol = 2
)
ps("./AM/ciliated_genes_pseudotime.pdf", w = 9, h = 6)




gene_module_df <- find_gene_modules(cds[degs, ], resolution = 1e-2, cores = n_jobs)
write_csv(gene_module_df, "./AM/AM2_Genes_Module.csv")

cell_group_df <- tibble::tibble(
  cell = row.names(colData(cds)),
  cell_group = colData(cds)$cell_type
)

agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))

pheatmap::pheatmap(agg_mat,
                   cluster_rows = TRUE, cluster_cols = TRUE,
                   scale = "column", clustering_method = "ward.D2",
                   fontsize = 6
)
ps("./AM/heatmap.pdf", w = 5, h = 16)
plot_cells(cds,
           genes = gene_module_df %>% dplyr::filter(module %in% c(200, 169, 138, 22,17,91)),
           group_cells_by = "partition",
           color_cells_by = "partition",
           show_trajectory_graph = FALSE
)
ps("./AM/module.pdf")
saveRDS(cds,file="cds.rds")
