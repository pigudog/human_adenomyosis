library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(stringr)
rm(list=ls())
options(stringsAsFactors = F)
load(file =  'AM_HVG.Rdata')

AM<-FindClusters(AM,resolution=0.8)
sce <- AM
sce_for_SingleR <- GetAssayData(sce, slot="data")
sce_for_SingleR
library(SingleR)
hpca.se <- HumanPrimaryCellAtlasData()
clusters=sce@meta.data$seurat_clusters
pred.hesc <- SingleR(test = sce_for_SingleR, ref = hpca.se, labels = hpca.se$label.main,
                     method = "cluster", clusters = clusters, 
                     assay.type.test = "logcounts", assay.type.ref = "logcounts")
table(pred.hesc$labels)
celltype = data.frame(ClusterID=rownames(pred.hesc), celltype=pred.hesc$labels, stringsAsFactors = F) 
sce@meta.data$singleR=celltype[match(clusters,celltype$ClusterID),'celltype']
p <-DimPlot(sce, reduction = "umap",label=TRUE ,group.by = "singleR")
ggsave(p,filename=paste0("AM",'/umap_singler.pdf'),width = 16,height = 12)
p <-DimPlot(sce, reduction = "umap",group.by = "singleR", split.by = "level0",label=TRUE )
ggsave(p,filename=paste0("AM",'/umap_singler_level0.pdf'),width = 16,height = 12)
