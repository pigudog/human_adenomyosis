library(reticulate)
library(patchwork)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(tidyverse)
library(Matrix)
library( magrittr)
library(sp)
library(monocle3)
options(stringsAsFactors = FALSE)
rm(list=ls())

load(file="AM2_HVG.Rdata")
# AM2, 第二个monocle：fibro-like，unkonwn2，SMA
library(ggsci)
library(Cairo)
AM2@meta.data$cell_type <- AM2@active.ident
color = c(pal_d3("category20")(20),
          pal_d3("category20b")(20),
          pal_d3("category20c")(20),
          pal_d3("category10")(10))

AM_SMM <- AM2[,AM2$cell_type %in% c("SMC","ESC")]
AM_SMM <- SCTransform(AM_SMM, vst.flavor = "v2")
AM_SMM <- RunPCA(AM_SMM, npcs=50, verbose=FALSE)
ElbowPlot(AM_SMM)
pc.num=1:20
AM_SMM <- RunUMAP(AM_SMM, dims=pc.num)
# 0.5
AM_SMM <- FindNeighbors(AM_SMM, dims = pc.num) 
AM_SMM <- FindClusters(AM_SMM,resolution=0.5)
p1 <- DimPlot(AM_SMM, reduction = 'umap', 
              group.by = "seurat_clusters",
              cols = color , #设置颜色
              repel = T,  #label不重叠
              label.size = 5, #调整label的大小
              label = TRUE,  #是否展示label
              pt.size = 0.5)
p2 <- DimPlot(AM_SMM, reduction = 'umap', group.by = "orig.ident",cols = color , 
              label = F, pt.size = 0.5)

ggsave(p1,filename=paste0("Figure5",'/picture_a.png'),width = 5,height = 4)
ggsave(p2,filename=paste0("Figure5",'/picture_b.png'),width = 5,height = 4)
AM_SMC.markers <- FindAllMarkers(object = AM_SMC,assay = "SCT", only.pos = TRUE, min.pct = 0.25,thresh.use = 0.25)
# wriet to placental_decidual/HVGmarkers_res_0.5
write.csv(AM_SMM.markers,file=paste0("Figure5",'/AM_SMMHVGmarkers_res_0.5.csv'))
library(dplyr) 
top10 <- AM_SMM.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
p3 <- DoHeatmap(AM_SMM,top10$gene,size=3)
ggsave(p3,filename=paste0("Figure5",'/AM_SMM_top10_heatmap.pdf'),width = 24,height = 18)

p_ADD <- DimPlot(AM_SMM, reduction = 'umap', 
              group.by = "cell_type",
              cols = color , 
              label = F, pt.size = 0.5)
ggsave(p_ADD,filename=paste0("Figure5",'/picture_c.png'),width = 5,height = 4)

pbmc <-AM_SMM
expression_matrix<-as(as.matrix(pbmc@assays$RNA@counts),"sparseMatrix")
cell_metadata<-pbmc@meta.data
gene_annotation<-data.frame(gene_short_name=row.names(expression_matrix), row.names=row.names(expression_matrix))
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

hvg <- pbmc@assays$SCT@var.features
pbmc_embed <- Seurat::Embeddings(pbmc, reduction = "umap")
cds = preprocess_cds(cds,norm_method="log",method="PCA",num_dim=100)
cds %<>% reduce_dimension(umap.fast_sgd = TRUE, cores = 48, reduction_method = "UMAP", preprocess_method = "PCA")
cds_embed <- cds@int_colData$reducedDims$UMAP
cds@int_colData$reducedDims$UMAP <- pbmc_embed[rownames(cds_embed),]
cds <- cluster_cells(cds,cores = 48,preprocess_method = "UMAP",cluster_method = "leiden",random_seed=1314)
cds <- learn_graph(cds)
cds <- order_cells(cds)
p4 <- plot_cells(cds,
                 group_label_size = 6,
                 color_cells_by = "pseudotime",
                 label_cell_groups = TRUE,
                 label_leaves = TRUE,
                 label_branch_points = TRUE,
                 graph_label_size = 0
)+labs(y="Percentage")+RotatedAxis()
ggsave(p4,filename=paste0("Figure5",'/picture_d.png'),width = 5,height = 4)


dea_res <- graph_test(cds, neighbor_graph = "principal_graph", cores = 48) %>%
  dplyr::filter(q_value < 0.05)
write_csv(dea_res, "./Figure5/Trajectory_genes.csv")
degs <- dea_res$gene_short_name

top_genes <- dea_res %>%
  dplyr::top_n(n = 10, morans_I) %>%
  dplyr::pull(gene_short_name) %>%
  as.character()
p5 <- plot_genes_in_pseudotime(cds[top_genes, ],
                         color_cells_by = "cell_type",
                         min_expr = 0.5, ncol = 2
)
# ps("./Figure5/Top_Genes_Jitterplot.pdf", w = 9, h = 6)
ggsave(p5,filename=paste0("Figure5",'/picture_e.png'),width = 9,height = 6)
p6 <- plot_cells(cds,
                genes = top_genes, show_trajectory_graph = FALSE,
                label_cell_groups = FALSE, label_leaves = FALSE
)
p6$facet$params$ncol <- 5
# ps("./Figure5/Top_Genes_Featureplot.pdf", plot = p)
ggsave(p6,filename=paste0("Figure5",'/picture_f.png'),width = 12,height = 6)
save(AM_SMM,AM_SMM.markers,cds,file = "./Figure5/AM_SMM.Rdata")



