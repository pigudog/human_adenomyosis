library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(stringr)
# devtools::install_github("sajuukLyu/ggunchull", type = "source")
# install.packages("tidydr")
# devtools::install_github(repo = "samuel-marsh/scCustomize", ref = "develop")
# remotes::install_github(repo = "samuel-marsh/scCustomize", ref = "develop")
# install.packages("scCustomize")
# devtools::install_github("xmc811/Scillus", ref = "development")
# library(Scillus)
library(cowplot)
library(dplyr)
library(tidydr)
library(stringr)
library(viridis)
library(scCustomize)
rm(list=ls())
options(stringsAsFactors = F)
load(file =  'AM_HVG.Rdata')
# Set color palette
# viridis提供了五种色带：
# viridis：option D，为默认色带，翠绿色；
# magma：option A，岩溶色；
# inferno：option B，火焰色；
# plasma：option C，血色；
# cividis：option E;
# "#0D0887FF" "#47039FFF" "#7301A8FF" "#9C179EFF" "#BD3786FF" "#D8576BFF" "#ED7953FF" "#FA9E3BFF" "#FDC926FF" "#F0F921FF"
pal <- viridis(n = 15, option = "D", direction = -1)
#之后使用scCustomize：：FeaturePlot_scCustom()函数来画一下FeaturePlot。
#因为通常画featureplot都是多个gene，所以在这里使用for循环。

#Continuous Palettes颜色

#viridis_plasma_dark_high
#viridis_plasma_light_high
#viridis_magma_dark_high#备选
#viridis_magma_light_high
#viridis_inferno_dark_high
#viridis_inferno_light_high
#viridis_dark_high
#使用Seurat::FeaturePlot()做一个常规的可视化
FeaturePlot(object = sce.all,features = gene, cols = pal, order = T)
ggsave(file="featureplot1.pdf",width=8, height=9)

# featureplot
# for (i in 1:length(genes_to_check)){
#   plots[[i]]=FeaturePlot_scCustom(seurat_object = AM, 
#                                   colors_use = pal, 
#                                   features = genes_to_check[i])+NoAxes() #+
#   # NoLegend()+
#     # theme(panel.border = element_rect(fill = NA,color = "black",
#     #                                   size=1.5,linetype = "solid"))
# }
# p<-wrap_plots(plots, ncol = 3)

  #theme(panel.border = element_rect(fill = NA,color = "black",
  #size=1.5,linetype = "solid"))


library(ComplexHeatmap)
library(circlize)
# Specify genes  
genes_to_check = c("ACTA2","DES",
                   "PDGFRB","GJA4","PDGFRA","COL1A2")
genes_to_check = c("CD14","MS4A7", # Mc
                   "NCAM1","GZMB", # NK
                   "CD3G","CD3D",  # T cells
                   "TPSB2","TPSAB1", # Mast
                   "PDGFRA","COL1A2", # ESC
                   "KRT8","KRT18", # Epi 
                   "PECAM1", "CDH5",  # Endo
                   "DES", "ACTA2",# SMC
                   "PDGFRB","GJA4" #Pericyte
)
pal <- viridis_magma_dark_high
# #给图例定义范围和颜色
# col_fun = colorRamp2(c(0,2,4), pal)
# lgd3 = Legend(at = c(0,2,4), col_fun = col_fun, 
#               nrow = 1,#border = "black",direction = "horizontal",
#               legend_height = unit(4, "cm")
# )
plots <- list()
for (i in 1:length(genes_to_check)){
  plots[[i]]=FeaturePlot_scCustom(seurat_object = AM, 
                                  colors_use = pal, 
                                  features = genes_to_check[i])+NoAxes()
}
p<-wrap_plots(plots, ncol = 3)
  
  # draw(lgd3,x = unit(0.99, "npc"), y = unit(0.4, "npc"),
  #                                      just = c("right", "bottom"));p
#可以调整图例的方向，高度，还有加不加border。
#用draw把图例画在合适的位置。
ggsave(p,filename=paste0("AM",'/second _fearureplot_sepcify.pdf'),width = 16,height = 12)

# All on Dotplot 
pal2 <- pal[50:100]
p <- DotPlot_scCustom(AM, features = genes_to_check,colors_use = pal2) + coord_flip()
ggsave(p,filename=paste0("AM",'/second_dotplot_sepcify.pdf'),width = 20,height = 12)

#  Need to look at the picture, determine the cell subsets:
celltype=data.frame(ClusterID=0:19,
                    celltype='unkown')
celltype[celltype$ClusterID %in% c(14),2]='NK/T cells' 
celltype[celltype$ClusterID %in% c(16),2]='Mast'
celltype[celltype$ClusterID %in% c(19),2]='Mac'
# not immune
celltype[celltype$ClusterID %in% c(7,17,18),2]='Epi'
celltype[celltype$ClusterID %in% c(4,5,13),2]='Endo'
celltype[celltype$ClusterID %in% c(9,10,15),2]='SMC'

# celltype[celltype$ClusterID %in% c(2,11,12),2]='Fibro-like'
# celltype[celltype$ClusterID %in% c(0,1,8,10,15),2]='Actin-Fibro'
# celltype[celltype$ClusterID %in% c(3),2]='Infra-Fibro'
celltype[celltype$ClusterID %in% c(3),2]='ESC'
celltype[celltype$ClusterID %in% c(0,1),2]='Unkown1'
celltype[celltype$ClusterID %in% c(6),2]='Unkown2'
# celltype[celltype$ClusterID %in% c(8,11,12),2]='fibroblast'
celltype[celltype$ClusterID %in% c(2,11,12),2]='Pericyte'
celltype[celltype$ClusterID %in% c(8),2]='VSMC'
table(celltype)


new.cluster.ids <- celltype$celltype
names(new.cluster.ids) <- levels(AM)
AM <- RenameIdents(AM, new.cluster.ids)
# legend
p<-DimPlot(AM, reduction = "umap", label = TRUE, pt.size = 1)
ggsave(p,filename=paste0("AM",'/umap_secondcluster_legend.pdf'),width = 16,height = 12)

p <-DimPlot(AM, reduction = "umap", split.by = "level0",label=TRUE )
ggsave(p,filename=paste0("AM",'/umap_secondcluster_legend_split.pdf'),width = 16,height = 12)

genes_to_check = c("PDGFRB","GJA4","RGS5", #Pericyte
                   "PDGFRA","COL1A2", # ESC
                   "PECAM1", "CDH5",  # Endo
                   "KRT8","KRT18", # Epi 
                   "DES", "ACTA2",# SMC
                   
                   # "COL1A1", "COL3A1","MDK", # Fibro
                   "GZMB", # NK
                   "CD3G","CD3D",  # T cells
                   "TPSB2","TPSAB1", # Mast
                    "CD14","MS4A7" # Mc
)
pal <- viridis_magma_dark_high
plots <- list()
for (i in 1:length(genes_to_check)){
  plots[[i]]=FeaturePlot_scCustom(seurat_object = AM,
                                  colors_use = pal,
                                  features = genes_to_check[i])+NoAxes()
}
p<-wrap_plots(plots, ncol = 3)
ggsave(p,filename=paste0("AM",'/second _fearureplot_sepcify.pdf'),width = 16,height = 12)
pal2 <- pal[50:100]
p <- DotPlot_scCustom(AM, features = genes_to_check,colors_use = pal2,flip_axes = T,x_lab_rotate = TRUE,
                      remove_axis_titles = FALSE) + coord_flip()
ggsave(p,filename=paste0("AM",'/second_dotplot_sepcify.pdf'),width = 8,height = 8)


save(AM,file = 'AM_after_cluster.Rdata')

















# RECLUSTER
p <- DotPlot(AM, features = genes_to_check) + coord_flip()
# immune
# Annotate Immune vs Nonimmune clusters
# At this point we dont care for a more detailed annotation as we will annotate immune and non-immune separately later
dat=p$data 
cd45=dat[dat$features.plot=='PTPRC',]
fivenum(cd45$avg.exp.scaled)
imm=cd45[cd45$avg.exp.scaled > 0,]$id
imm
AM@meta.data$immune_annotation <-ifelse(AM@meta.data$SCT_snn_res.0.5  %in% imm ,'immune','non-immune')
# MAke a table 
table(AM@meta.data$immune_annotation)
# Make and save relevant plots 
p<-DimPlot(object=AM,group.by='immune_annotation')
ggsave(p,filename=paste0("AM",'/umap_immune-or-not.pdf'),width = 16,height = 12)
phe=AM@meta.data
save(phe,file = 'phe-of-immune-or-not.Rdata')

# fibroblast
p <- DotPlot(AM, features = genes_to_check) + coord_flip()
dat=p$data 
COL1A1=dat[dat$features.plot=='COL1A1',]
fivenum(COL1A1$avg.exp.scaled)
fibro=COL1A1[COL1A1$avg.exp.scaled > -0.3,]$id
fibro
AM@meta.data$fibroblast_annotation <-ifelse(AM@meta.data$SCT_snn_res.0.5  %in% fibro ,'fibroblast','non-fibroblast')
# MAke a table 
table(AM@meta.data$fibroblast_annotation)
# Make and save relevant plots 
p<-DimPlot(object=AM,group.by='fibroblast_annotation')
ggsave(p,filename=paste0("AM",'/umap_fibroblast_annotation.pdf'),width = 16,height = 12)
phe=AM@meta.data
save(phe,file = 'phe-of-fibroblast_annotation.Rdata')


#  Need to look at the picture, determine the cell subsets:
celltype=data.frame(ClusterID=0:19,
                    celltype='unkown')
celltype[celltype$ClusterID %in% c(14),2]='T cells' 
celltype[celltype$ClusterID %in% c(0,1,2,3,11,12,16),2]='Fibroblast'
celltype[celltype$ClusterID %in% c(7,17,18),2]='Ep'
celltype[celltype$ClusterID %in% c(4,5,13),2]='En'
celltype[celltype$ClusterID %in% c(8,15),2]='ILC'

celltype[celltype$ClusterID %in% c(6,9,10),2]='SMC'
celltype[celltype$ClusterID %in% c(19),2]='macrophage'
table(celltype)


new.cluster.ids <- celltype$celltype
names(new.cluster.ids) <- levels(AM)
AM <- RenameIdents(AM, new.cluster.ids)
# legend
p<-DimPlot(AM, reduction = "umap", label = TRUE, pt.size = 0.5)
ggsave(p,filename=paste0("AM",'/umap_secondcluster_legend.pdf'),width = 16,height = 12)

p <-DimPlot(AM, reduction = "umap", split.by = "level0",label=TRUE )
ggsave(p,filename=paste0("AM",'/umap_secondcluster_legend_split.pdf'),width = 16,height = 12)



# 设置可能用到的主题
library(ggpubr)
library(ggsignif)
my_comparisons <- list(c("CTRL","EC"),c("EC","EM"),c("CTRL","EM"))
theme.set2 = theme(axis.title.x=element_blank())
# 设置绘图元素
# PDGFRA和PDGFC
plot.featrures = c("MMP1","CSF3","CXCL8","IGF1")
group = "orig.ident"
# 质控前小提琴图
plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(AM, group.by=group, pt.size = 0,
                       features = plot.featrures[i])+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    stat_compare_means(comparisons = my_comparisons) +
    ylim(-1, 8)+
    geom_boxplot(position=position_dodge(0.6),
                 size=0.4,
                 width=0.1,
                 # fill="gold",
                 # color="blue",
                 # outlier.color = "blue",
                 # outlier.fill = "red",
                 outlier.shape = 19,
                 outlier.size = 1.5,
                 outlier.stroke = 0.5,
                 outlier.alpha = 45,
                 notch = F,
                 notchwidth = 0.5)+
  theme.set2 + NoLegend()}
p <- wrap_plots(plots = plots, nrow=2)     
ggsave(p,filename=paste0("AM",'/umap_MMP1.pdf'),width = 16,height = 12)


##############################
# THRID CLUSTER
# - B cells - CD79A, IGKC
# - Progenitor-like cells,PC - TOP2A, UBE2C
# - Uncil Ep - CD24, EPCAM
# - Act fibro - LEFTY2, SVIL, KCNMA1, DES, TAGLN,ACTA2
# - Macrophage - MS4A7，CD14
# - NK - GNLY, NKG7
# - Pericyte - NOTCH3, COX4I2,PDGFRB,MCAM,
# - SMC - CNN1
# - T cells - CCL5, CD3D, CD69
# - En - PECAM1, ICAM
# - Fibroblast-like cells - MMP7,MMP11,SFRP4, LUM, COL3A1,MDK
# - MAST - "TPSB2","TPSAB1"
# Specify genes  
# FBs共享ECM基因如COL1A1、COL1A2和COL3A1(图B)和平滑肌肌动蛋白ACTA2的表达。
# FB1和FB3通过表达DES、MFAP4、OGN和S100A4基因表现出与肌成纤维细胞的特征相似性。
# FB3还表达促炎基因，如IL6、PTGDS、CFD、CXCL2和BDKRB1。
# FB2s的特点是独特表达REN和AGTR1，这些基因与血压、钠和体液平衡的调节有关，以及IGFBP7和AREG。
genes_to_check = c("CD14","MS4A7", # Mc
                   "NCAM1","GZMB", # NK
                   "CD3G","CD3D",  # T cells
                   "TPSB2","TPSAB1", # Mast
                   "PDGFRA","PDPN", # CSC
                   "KRT8","KRT18", # Epi 
                   "PECAM1", "CDH5", "LYVE1", # Endo
                   "CNN1", "DES", "TAGLN","ACTA2",#SMC
                   "RPL37A","RPL27A","RPL37", #VMC
                   "SVIL","KCNMA1", #Actin Fibro
                   "COL1A1", "COL3A1","MDK", # Fibro
                   "IL6","PTGDS","CFD","CXCL2"# Infra Fibro
                   )
# featureplot
p <- FeaturePlot(AM, features = genes_to_check)
ggsave(p,filename=paste0("AM",'/third _fearureplot_sepcify.pdf'),width = 16,height = 12)
# All on Dotplot 
p <- DotPlot(AM, features = genes_to_check) + coord_flip()
ggsave(p,filename=paste0("AM",'/third_second_dotplot_sepcify.pdf'),width = 16,height = 12)

#  Need to look at the picture, determine the cell subsets:
celltype=data.frame(ClusterID=0:19,
                    celltype='unkown')
# immune
celltype[celltype$ClusterID %in% c(14),2]='NK/T cells' 
celltype[celltype$ClusterID %in% c(16),2]='Mast'
celltype[celltype$ClusterID %in% c(19),2]='Mac'
# not immune
celltype[celltype$ClusterID %in% c(7,17,18),2]='Epi'
celltype[celltype$ClusterID %in% c(4,5,13),2]='Endo'
celltype[celltype$ClusterID %in% c(6,9),2]='SMC'

# celltype[celltype$ClusterID %in% c(2,11,12),2]='Fibro-like'
# celltype[celltype$ClusterID %in% c(0,1,8,10,15),2]='Actin-Fibro'
# celltype[celltype$ClusterID %in% c(3),2]='Infra-Fibro'
celltype[celltype$ClusterID %in% c(0,1,2,3,8,10,11,12,15),2]='Fibro-like'
celltype[celltype$ClusterID %in% c(0,1),2]='Unkown1'
celltype[celltype$ClusterID %in% c(2,11),2]='Unkown2'
table(celltype)


new.cluster.ids <- celltype$celltype
names(new.cluster.ids) <- levels(AM)
AM <- RenameIdents(AM, new.cluster.ids)
# legend
p<-DimPlot(AM, reduction = "umap", label = TRUE, pt.size = 0.5)
ggsave(p,filename=paste0("AM",'/umap_thirdcluster_legend.pdf'),width = 12,height = 9)

p <-DimPlot(AM, reduction = "umap", split.by = "level0",label=TRUE )
ggsave(p,filename=paste0("AM",'/umap_thirdcluster_legend_split.pdf'),width = 16,height = 9)
genes_to_check = c(
                  # "SVIL","KCNMA1", #Actin Fibro
                   "COL1A1", "COL3A1", # Fibro
                   "APOD",  "FGF7", 
                   # "PTGDS","CFD",# Infra Fibro
                   "PECAM1", "CDH5", "LYVE1", # Endo
                   "CNN1", "DES", "TAGLN","ACTA2",#SMC
                   "KRT8","KRT18", # Epi 
                   "GZMB", # NK
                   "CD3G","CD3D",  # T cells
                   "TPSB2","TPSAB1","CXCL2", # Mast
                    "CD14","MS4A7" # Mc
)
# featureplot
# All on Dotplot 
p <- DotPlot(AM, features = genes_to_check,dot.scale = 7) + coord_flip()
ggsave(p,filename=paste0("AM",'/third_dotplot_sepcify.pdf'),width = 8,height = 6)


save(AM,file = 'AM_after_cluster.Rdata')

