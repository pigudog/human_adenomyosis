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
#设置group顺序
# AM@meta.data$cell_type <- factor(AM@active.ident )
df=AM2@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  cbind(cell_type=AM2@meta.data$cell_type)

#用stat_ellipse加置信区间
p_UMAP1=ggplot(df, aes(UMAP_1, UMAP_2, color=cell_type))+
  geom_point(size = 1) + #细胞亚群点的大小
  scale_color_manual(values=color)+
  theme(panel.border = element_blank(), #边框
        axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),
        # panel.grid = element_blank(),
        # panel.grid.major = element_blank(), #网格线
        # panel.grid.minor = element_blank(), #网格线
        panel.background = element_rect(fill = 'white'), #背景色
        plot.background=element_rect(fill="white"),
        legend.title = element_blank(), #去掉legend.title 
        legend.key=element_rect(fill='white'), 
        legend.text = element_text(size=14), #设置legend标签的大小
        legend.key.size=unit(0.6,'cm') )+
  stat_ellipse(aes(x = UMAP_1, y = UMAP_2, fill = cell_type),
               geom = "polygon",#多边形
               linetype=2, # 椭圆，虚线
               alpha = 0.05,  #椭圆不透明度
               show.legend = FALSE, #去掉椭圆对应的图例
               level = 0.93)+ #置信区间
  guides(fill= guide_legend(override.aes = list(size = 4)))+ #图例点的大小
  scale_color_manual(values=color) 

#theme_dr加左下角的坐标
#呈现效果如图，如果还想更加精确的加上细胞群外圈的框框的话，PPT和AI或许是更好的选择。
#theme_dr
p_UMAP2<-p_UMAP1+theme_dr(xlength = 0.22, ylength = 0.22,
                          arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")
)+theme(panel.grid = element_blank())
#先提取出tSNE_1的位置信息
celltype_position <- df %>%
  group_by(cell_type) %>%
  dplyr::summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2))

library(ggrepel)
CairoPDF("./Figure4/UMAP_label.pdf", width=8, height=7) 
p_UMAP3<-p_UMAP2+geom_point(aes(color=factor(cell_type)), size=3)+
  geom_label_repel(data = celltype_position,aes(label=cell_type), 
                   fontface="bold",point.padding=unit(0.1, "lines"))+theme(legend.position = "none");p_UMAP3

dev.off()


#分样本
sample_table <- as.data.frame(table(AM2@meta.data$orig.ident,AM2@meta.data$cell_type))
names(sample_table) <- c("Samples","celltype","CellNumber")

plot_sample<-ggplot(sample_table,aes(x=Samples,weight=CellNumber,fill=celltype))+
  geom_bar(position="fill")+
  scale_fill_manual(values=color) + 
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        axis.line.x = element_line(colour = "black") ,
        axis.line.y = element_line(colour = "black") ,
        axis.title = element_text(lineheight=2, face="bold", hjust=2, size =15),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15),
        plot.title = element_text(size = 20)
  )+labs(y="Percentage")+RotatedAxis()
ggsave(plot_sample,filename=paste0("Figure4",'/picture_b.png'),width = 6,height = 8)


## pictureD
p1 <- p_UMAP3
p2 <- DimPlot(AM2, reduction = 'umap', 
              group.by = "seurat_clusters",
              cols = color , #设置颜色
              repel = T,  #label不重叠
              label.size = 5, #调整label的大小
              label = TRUE,  #是否展示label
              pt.size = 0.5)
p3 <- DimPlot(AM2, reduction = 'umap', group.by = "orig.ident",cols = color , 
              label = F, pt.size = 0.5)

ggsave(p2,filename=paste0("Figure4",'/picture_c.png'),width = 5,height = 4)
ggsave(p3,filename=paste0("Figure4",'/picture_d.png'),width = 5,height = 4)
# saveRDS(cds,file="cds.rds")
pbmc <-AM2
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
ggsave(p4,filename=paste0("Figure4",'/picture_e.png'),width = 5,height = 4)

cds <- order_cells(cds)
p5 <- plot_cells(cds,
                 group_label_size = 6,
                 color_cells_by = "pseudotime",
                 label_cell_groups = TRUE,
                 label_leaves = TRUE,
                 label_branch_points = TRUE,
                 graph_label_size = 0
)+labs(y="Percentage")+RotatedAxis()
ggsave(p5,filename=paste0("Figure4",'/picture_f.png'),width = 5,height = 4)

############ 分步去做monocle3的基因查看？
# dea_res <- graph_test(cds, neighbor_graph = "principal_graph", cores = n_jobs) %>%
#   dplyr::filter(q_value < 0.05)
# write_csv(dea_res, "./AM2/Trajectory_genes.csv")
degs <- Trajectory_genes$gene_short_name

top_genes <- Trajectory_genes %>%
  dplyr::top_n(n = 10, morans_I) %>%
  dplyr::pull(gene_short_name) %>%
  as.character()
p4 <- plot_genes_in_pseudotime(cds[top_genes, ],
                         color_cells_by = "cell_type",
                         min_expr = 0.5, ncol = 2
)
# Top_Genes_Jitterplot.pdf
ggsave(p4,filename=paste0("Figure4",'/picture_d.png'),width = 8,height = 11)

p5 <- plot_cells(cds,
                genes = top_genes, show_trajectory_graph = FALSE,
                label_cell_groups = FALSE, label_leaves = FALSE
)
p5$facet$params$ncol <- 5
# Top_Genes_Featureplot.pdf
ggsave(p4,filename=paste0("Figure4",'/picture_e.png'),width = 10,height = 5)
