rm(list=ls())
library(Seurat)
library(tidyverse)
library(Matrix)
library( magrittr)
library(sp)
library(monocle3)
load(file = 'AM_after_cluster.Rdata')
# devtools::install_github("TheHumphreysLab/plot1cell")
# ## or the development version, devtools::install_github("HaojiaWu/plot1cell")
# 
# ## You might need to install the dependencies below if they are not available in your R library.
# bioc.packages <- c("biomaRt","GenomeInfoDb","EnsDb.Hsapiens.v86","GEOquery","simplifyEnrichment","ComplexHeatmap")
# BiocManager::install(bioc.packages)
# 
# dev.packages <- c("chris-mcginnis-ucsf/DoubletFinder","Novartis/hdf5r","mojaveazure/loomR")
# devtools::install_github(dev.packages)
# ## If you can't get the hdf5r package installed, please see the fix here:
# ## https://github.com/hhoeflin/hdf5r/issues/94

library(plot1cell)
###Check and see the meta data info on your Seurat object
colnames(AM@meta.data)  

###Prepare data for ploting 准备圈图数据
circ_data <- prepare_circlize_data(AM, scale = 0.8 )
set.seed(1234)
# 设置细胞分群信息的颜色
cluster_colors<-rand_color(length(levels(AM)))
group_colors<-rand_color(length(names(table(AM$Group))))
rep_colors<-rand_color(length(names(table(AM$orig.ident))))

###plot and save figures
# 绘制细胞分群圈图
png(filename =  'circlize_plot.png', width = 6, height = 6,units = 'in', res = 300)
plot_circlize(circ_data,do.label = T, pt.size = 0.01, col.use = cluster_colors ,bg.color = 'white', kde2d.n = 200, repel = T, label.cex = 0.6)
# 添加细胞群注释信息
add_track(circ_data, group = "Group", colors = group_colors, track_num = 2) ## can change it to one of the columns in the meta data of your seurat object
add_track(circ_data, group = "orig.ident",colors = rep_colors, track_num = 3) ## can change it to one of the columns in the meta data of your seurat object
dev.off()
