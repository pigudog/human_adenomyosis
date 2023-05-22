BiocManager::install("DESeq2")
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
BiocManager::install("DESeq2")
BiocManager::install("clusterProfiler")
library(DESeq2)
library(edgeR)
library(limma)
library(clusterProfiler)
gene_name <- rownames(AM.markers)
#将基因名转换为entrez_id
id = bitr(gene_name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
ego <- enrichGO(OrgDb="org.Hs.eg.db", gene = id$ENTREZID, ont = "MF", pvalueCutoff = 0.01, readable= TRUE) #GO富集分析
dotplot(ego,showCategory=10,title="Enrichment GO Top10") #泡泡图
barplot(ego, showCategory=20,title="EnrichmentGO")  #柱状图
library(clusterProfile)
plotGOgraph(ego) 	#GO图，看不清楚可以尝试左上角另存为pdf
ekk <- enrichKEGG(gene= id$ENTREZID,organism  = 'hsa', qvalueCutoff = 0.05)	 #KEGG富集分析
dotplot(ekk,font.size=8)	# 画气泡图
browseKEGG(ekk,'mmu01100')	
