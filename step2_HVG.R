library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(stringr)
library(harmony)
rm(list=ls())
options(stringsAsFactors = F)
load('AM.Rdata')

AM <- SCTransform(AM, vst.flavor = "v2")
AM <- RunPCA(AM, npcs=50, verbose=FALSE)
ElbowPlot(AM)
pc.num=1:20
AM <- RunUMAP(AM, dims=pc.num)
# 0.5
AM <- FindNeighbors(AM, dims = pc.num) 
AM<-FindClusters(AM,resolution=0.5)
p <- DimPlot(AM2, group.by = "orig.ident")
ggsave("AM2/UMAP_Samples_before_singler.pdf", p, width = 8, height = 6)
p <- DimPlot(AM, group.by = "orig.ident", split.by = "orig.ident", ncol = 4)
ggsave("AM/UMAP_Samples_Split_before_singler.pdf", p, width = 18, height = 12)
p <- DimPlot(AM,  label = T, group.by = "seurat_clusters")
ggsave("AM/UMAP_cluster_0.5.pdf", p, width = 18, height = 12)

AM.markers <- FindAllMarkers(object = AM,assay = "SCT", only.pos = TRUE, min.pct = 0.25,thresh.use = 0.25)
# wriet to placental_decidual/HVGmarkers_res_0.5
write.csv(AM.markers,file=paste0("AM",'/HVGmarkers_res_0.5.csv'))


library(dplyr) 
top10 <- AM.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
p <- DoHeatmap(AM,top10$gene,size=3)
ggsave(p,filename=paste0("AM",'/top10_heatmap.pdf'),width = 24,height = 18)

save(AM,AM.markers,file = 'AM_HVG.Rdata')
