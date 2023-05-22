library(clusterProfiler)
# library(org.Mm.eg.db) ##加载小鼠
library(org.Hs.eg.db) ##加载人类
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
load(file = 'AM_after_cluster.Rdata')
AM@meta.data$cell_type <- AM@active.ident

AM_SMC <- AM[,AM$cell_type %in% "SMC"]
# CTRL VS EM
# future::plan("multicore", workers = 8)
SMC_CTRLEM <- FindMarkers(AM_SMC, ident.1 = 'EM',ident.2 = 'CTRL', 
                            group.by = 'orig.ident',logfc.threshold = 0,min.pct = 0)
sig_dge.all <- subset(SMC_CTRLEM, p_val_adj<0.05&abs(avg_log2FC)>0.7) #所有差异基因
# # BP, CC和MF三种通路都一起富集
# ego_ALL <- enrichGO(gene          = row.names(sig_dge.all),
#                     #universe     = row.names(dge.celltype),
#                     OrgDb         = 'org.Hs.eg.db',
#                     keyType       = 'SYMBOL',
#                     ont           = "ALL",  #设置为ALL时BP, CC, MF都计算
#                     pAdjustMethod = "BH",
#                     pvalueCutoff  = 0.01,
#                     qvalueCutoff  = 0.05)
# ego_all <- data.frame(ego_ALL)
# 分别对BP, CC和MF进行富集
ego_CC <- enrichGO(gene          = row.names(sig_dge.all),
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)
# ego_cc <- data.frame(ego_CC)
ego_MF <- enrichGO(gene          = row.names(sig_dge.all),
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)
# ego_mf <- data.frame(ego_MF)
ego_BP <- enrichGO(gene          = row.names(sig_dge.all),
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05) 
# ego_bp <- data.frame(ego_BP)
# p_BP <- barplot(ego_BP,showCategory = 10) + ggtitle("barplot for Biological process")
# p_CC <- barplot(ego_CC,showCategory = 10) + ggtitle("barplot for Cellular component")
# p_MF <- barplot(ego_MF,showCategory = 10) + ggtitle("barplot for Molecular function")
# plotc <- p_BP/p_CC/p_MF
#但是有时候我们富集到的通路太多了，不方便绘制图片，可以选择部分绘制，这时候跳过第4步，来到第5步
display_number = c(10, 10, 10)#这三个数字分别代表选取的BP、CC、MF的通路条数，这个自己设置就行了
ego_result_BP <- as.data.frame(ego_BP)[1:display_number[1], ]
ego_result_CC <- as.data.frame(ego_CC)[1:display_number[2], ]
ego_result_MF <- as.data.frame(ego_MF)[1:display_number[3], ]

##将以上我们摘取的部分通路重新组合成数据框
go_enrich_df <- data.frame(
  ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),                         
  Description=c(ego_result_BP$Description,ego_result_CC$Description,ego_result_MF$Description),
  GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
  type=factor(c(rep("biological process", display_number[1]), 
                rep("cellular component", display_number[2]),
                rep("molecular function", display_number[3])), 
              levels=c("biological process", "cellular component","molecular function" )))

##通路的名字太长了，我选取了通路的前五个单词作为通路的名字
# for(i in 1:nrow(go_enrich_df)){
#   description_splite=strsplit(go_enrich_df$Description[i],split = " ")
#   description_collapse=paste(description_splite[[1]][1:5],collapse = " ") #这里的5就是指5个单词的意思，可以自己更改
#   go_enrich_df$Description[i]=description_collapse
#   go_enrich_df$Description=gsub(pattern = "NA","",go_enrich_df$Description)
# }

##开始绘制GO柱状图
###横着的柱状图
go_enrich_df$type_order=factor(rev(as.integer(rownames(go_enrich_df))),labels=rev(go_enrich_df$Description))#这一步是必须的，为了让柱子按顺序显示，不至于很乱
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")#设定颜色

plotc <- ggplot(data=go_enrich_df, aes(x=type_order,y=GeneNumber, fill=type)) + #横纵轴取值
  geom_bar(stat="identity", width=0.8) + #柱状图的宽度，可以自己设置
  scale_fill_manual(values = COLS) + ###颜色
  coord_flip() + ##这一步是让柱状图横过来，不加的话柱状图是竖着的
  xlab("GO term") +
  ylab("Gene_Number") +
  labs(title = "The Most Enriched GO Terms")+
  theme_bw()

ggsave('./Figure6/CTRLEM_SMC_enrichGO.png', plotc, width = 12,height = 14)

# KEGG
genelist <- bitr(row.names(sig_dge.all), fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
genelist <- pull(genelist,ENTREZID)               
ekegg <- enrichKEGG(gene = genelist, organism = 'hsa',pvalueCutoff = 0.8)
p1 <- barplot(ekegg, showCategory=15)
p2 <- dotplot(ekegg, showCategory=15)
plotc = p1/p2
ggsave("./Figure6/CTRLEM_SMC_enrichKEGG.png", plot = plotc, width = 12, height = 14)

#########################
# EC VS EM
# future::plan("multicore", workers = 8)
SMC_ECEM <- FindMarkers(AM_SMC, ident.1 = 'EC',ident.2 = 'EM', 
                          group.by = 'orig.ident',logfc.threshold = 0,min.pct = 0)
sig_dge.all <- subset(SMC_ECEM, p_val_adj<0.05&abs(avg_log2FC)>0.8) #所有差异基因
# # BP, CC和MF三种通路都一起富集
# ego_ALL <- enrichGO(gene          = row.names(sig_dge.all),
#                     #universe     = row.names(dge.celltype),
#                     OrgDb         = 'org.Hs.eg.db',
#                     keyType       = 'SYMBOL',
#                     ont           = "ALL",  #设置为ALL时BP, CC, MF都计算
#                     pAdjustMethod = "BH",
#                     pvalueCutoff  = 0.01,
#                     qvalueCutoff  = 0.05)
# ego_all <- data.frame(ego_ALL)
# 分别对BP, CC和MF进行富集
ego_CC <- enrichGO(gene          = row.names(sig_dge.all),
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)
# ego_cc <- data.frame(ego_CC)
ego_MF <- enrichGO(gene          = row.names(sig_dge.all),
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)
# ego_mf <- data.frame(ego_MF)
ego_BP <- enrichGO(gene          = row.names(sig_dge.all),
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05) 
# ego_bp <- data.frame(ego_BP)
# p_BP <- barplot(ego_BP,showCategory = 10) + ggtitle("barplot for Biological process")
# p_CC <- barplot(ego_CC,showCategory = 10) + ggtitle("barplot for Cellular component")
# p_MF <- barplot(ego_MF,showCategory = 10) + ggtitle("barplot for Molecular function")
# plotc <- p_BP/p_MF
#但是有时候我们富集到的通路太多了，不方便绘制图片，可以选择部分绘制，这时候跳过第4步，来到第5步
display_number = c(10, 0, 10)#这三个数字分别代表选取的BP、CC、MF的通路条数，这个自己设置就行了
ego_result_BP <- as.data.frame(ego_BP)[1:display_number[1], ]
ego_result_CC <- as.data.frame(ego_CC)[1:display_number[2], ]
ego_result_MF <- as.data.frame(ego_MF)[1:display_number[3], ]

##将以上我们摘取的部分通路重新组合成数据框
go_enrich_df <- data.frame(
  ID=c(ego_result_BP$ID,  ego_result_MF$ID),                         
  Description=c(ego_result_BP$Description,ego_result_MF$Description),
  GeneNumber=c(ego_result_BP$Count, ego_result_MF$Count),
  type=factor(c(rep("biological process", display_number[1]), 
                # rep("cellular component", display_number[2]),
                rep("molecular function", display_number[3])), 
              levels=c("biological process","molecular function" )))
##开始绘制GO柱状图
###横着的柱状图
go_enrich_df$type_order=factor(rev(as.integer(rownames(go_enrich_df))),labels=rev(go_enrich_df$Description))#这一步是必须的，为了让柱子按顺序显示，不至于很乱
COLS <- c("#66C3A5",  "#FD8D62")#设定颜色

plotc <- ggplot(data=go_enrich_df, aes(x=type_order,y=GeneNumber, fill=type)) + #横纵轴取值
  geom_bar(stat="identity", width=0.8) + #柱状图的宽度，可以自己设置
  scale_fill_manual(values = COLS) + ###颜色
  coord_flip() + ##这一步是让柱状图横过来，不加的话柱状图是竖着的
  xlab("GO term") +
  ylab("Gene_Number") +
  labs(title = "The Most Enriched GO Terms")+
  theme_bw()
ggsave('./Figure6/ECEM_SMC_enrichGO.png', plotc, width = 12,height = 14)

# KEGG
genelist <- bitr(row.names(sig_dge.all), fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
genelist <- pull(genelist,ENTREZID)               
ekegg <- enrichKEGG(gene = genelist, organism = 'hsa',pvalueCutoff = 0.8)
p1 <- barplot(ekegg, showCategory=15)
p2 <- dotplot(ekegg, showCategory=15)
plotc = p1/p2
ggsave("./Figure6/ECEM_SMC_enrichKEGG.png", plot = plotc, width = 12, height = 14)






#############################################################################################################3


AM_ESC <- AM[,AM$cell_type %in% "ESC"]
# CTRL VS EM
# future::plan("multicore", workers = 8)
ESC_CTRLEM <- FindMarkers(AM_ESC, ident.1 = 'EM',ident.2 = 'CTRL', 
                          group.by = 'orig.ident',logfc.threshold = 0,min.pct = 0)
sig_dge.all <- subset(ESC_CTRLEM, p_val_adj<0.05&abs(avg_log2FC)>0.8) #所有差异基因
# # BP, CC和MF三种通路都一起富集
# ego_ALL <- enrichGO(gene          = row.names(sig_dge.all),
#                     #universe     = row.names(dge.celltype),
#                     OrgDb         = 'org.Hs.eg.db',
#                     keyType       = 'SYMBOL',
#                     ont           = "ALL",  #设置为ALL时BP, CC, MF都计算
#                     pAdjustMethod = "BH",
#                     pvalueCutoff  = 0.01,
#                     qvalueCutoff  = 0.05)
# ego_all <- data.frame(ego_ALL)
# 分别对BP, CC和MF进行富集
ego_CC <- enrichGO(gene          = row.names(sig_dge.all),
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)
# ego_cc <- data.frame(ego_CC)
ego_MF <- enrichGO(gene          = row.names(sig_dge.all),
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)
# ego_mf <- data.frame(ego_MF)
ego_BP <- enrichGO(gene          = row.names(sig_dge.all),
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05) 
# ego_bp <- data.frame(ego_BP)

# p_BP <- barplot(ego_BP,showCategory = 10) + ggtitle("barplot for Biological process")
# p_CC <- barplot(ego_CC,showCategory = 10) + ggtitle("barplot for Cellular component")
# p_MF <- barplot(ego_MF,showCategory = 10) + ggtitle("barplot for Molecular function")
# plotc <- p_BP/p_CC/p_MF
#但是有时候我们富集到的通路太多了，不方便绘制图片，可以选择部分绘制，这时候跳过第4步，来到第5步
display_number = c(10, 10, 10)#这三个数字分别代表选取的BP、CC、MF的通路条数，这个自己设置就行了
ego_result_BP <- as.data.frame(ego_BP)[1:display_number[1], ]
ego_result_CC <- as.data.frame(ego_CC)[1:display_number[2], ]
ego_result_MF <- as.data.frame(ego_MF)[1:display_number[3], ]

##将以上我们摘取的部分通路重新组合成数据框
go_enrich_df <- data.frame(
  ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),                         
  Description=c(ego_result_BP$Description,ego_result_CC$Description,ego_result_MF$Description),
  GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
  type=factor(c(rep("biological process", display_number[1]), 
                rep("cellular component", display_number[2]),
                rep("molecular function", display_number[3])), 
              levels=c("biological process", "cellular component","molecular function" )))

##通路的名字太长了，我选取了通路的前五个单词作为通路的名字
# for(i in 1:nrow(go_enrich_df)){
#   description_splite=strsplit(go_enrich_df$Description[i],split = " ")
#   description_collapse=paste(description_splite[[1]][1:5],collapse = " ") #这里的5就是指5个单词的意思，可以自己更改
#   go_enrich_df$Description[i]=description_collapse
#   go_enrich_df$Description=gsub(pattern = "NA","",go_enrich_df$Description)
# }

##开始绘制GO柱状图
###横着的柱状图
go_enrich_df$type_order=factor(rev(as.integer(rownames(go_enrich_df))),labels=rev(go_enrich_df$Description))#这一步是必须的，为了让柱子按顺序显示，不至于很乱
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")#设定颜色

plotc<-ggplot(data=go_enrich_df, aes(x=type_order,y=GeneNumber, fill=type)) + #横纵轴取值
  geom_bar(stat="identity", width=0.8) + #柱状图的宽度，可以自己设置
  scale_fill_manual(values = COLS) + ###颜色
  coord_flip() + ##这一步是让柱状图横过来，不加的话柱状图是竖着的
  xlab("GO term") +
  ylab("Gene_Number") +
  labs(title = "The Most Enriched GO Terms")+
  theme_bw()


ggsave('./Figure6/CTRLEM_ESC_enrichGO.png', plotc, width = 12,height = 14)

# KEGG
genelist <- bitr(row.names(sig_dge.all), fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
genelist <- pull(genelist,ENTREZID)               
ekegg <- enrichKEGG(gene = genelist, organism = 'hsa',pvalueCutoff = 0.8)
p1 <- barplot(ekegg, showCategory=15)
p2 <- dotplot(ekegg, showCategory=15)
plotc = p1/p2
ggsave("./Figure6/CTRLEM_ESC_enrichKEGG.png", plot = plotc, width = 12, height = 14)

#########################
# EC VS EM
# future::plan("multicore", workers = 8)
ESC_ECEM <- FindMarkers(AM_ESC, ident.1 = 'EC',ident.2 = 'EM', 
                          group.by = 'orig.ident',logfc.threshold = 0,min.pct = 0)
sig_dge.all <- subset(ESC_ECEM, p_val_adj<0.05&abs(avg_log2FC)>0.8) #所有差异基因
# # BP, CC和MF三种通路都一起富集
# ego_ALL <- enrichGO(gene          = row.names(sig_dge.all),
#                     #universe     = row.names(dge.celltype),
#                     OrgDb         = 'org.Hs.eg.db',
#                     keyType       = 'SYMBOL',
#                     ont           = "ALL",  #设置为ALL时BP, CC, MF都计算
#                     pAdjustMethod = "BH",
#                     pvalueCutoff  = 0.01,
#                     qvalueCutoff  = 0.05)
# ego_all <- data.frame(ego_ALL)
# 分别对BP, CC和MF进行富集
ego_CC <- enrichGO(gene          = row.names(sig_dge.all),
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)
ego_cc <- data.frame(ego_CC)
ego_MF <- enrichGO(gene          = row.names(sig_dge.all),
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)
ego_mf <- data.frame(ego_MF)
ego_BP <- enrichGO(gene          = row.names(sig_dge.all),
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05) 
ego_bp <- data.frame(ego_BP)
display_number = c(10, 10, 10)#这三个数字分别代表选取的BP、CC、MF的通路条数，这个自己设置就行了
ego_result_BP <- as.data.frame(ego_BP)[1:display_number[1], ]
ego_result_CC <- as.data.frame(ego_CC)[1:display_number[2], ]
ego_result_MF <- as.data.frame(ego_MF)[1:display_number[3], ]

##将以上我们摘取的部分通路重新组合成数据框
go_enrich_df <- data.frame(
  ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),                         
  Description=c(ego_result_BP$Description,ego_result_CC$Description,ego_result_MF$Description),
  GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
  type=factor(c(rep("biological process", display_number[1]), 
                rep("cellular component", display_number[2]),
                rep("molecular function", display_number[3])), 
              levels=c("biological process", "cellular component","molecular function" )))

##通路的名字太长了，我选取了通路的前五个单词作为通路的名字
# for(i in 1:nrow(go_enrich_df)){
#   description_splite=strsplit(go_enrich_df$Description[i],split = " ")
#   description_collapse=paste(description_splite[[1]][1:5],collapse = " ") #这里的5就是指5个单词的意思，可以自己更改
#   go_enrich_df$Description[i]=description_collapse
#   go_enrich_df$Description=gsub(pattern = "NA","",go_enrich_df$Description)
# }

##开始绘制GO柱状图
###横着的柱状图
go_enrich_df$type_order=factor(rev(as.integer(rownames(go_enrich_df))),labels=rev(go_enrich_df$Description))#这一步是必须的，为了让柱子按顺序显示，不至于很乱
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")#设定颜色

plotc<-ggplot(data=go_enrich_df, aes(x=type_order,y=GeneNumber, fill=type)) + #横纵轴取值
  geom_bar(stat="identity", width=0.8) + #柱状图的宽度，可以自己设置
  scale_fill_manual(values = COLS) + ###颜色
  coord_flip() + ##这一步是让柱状图横过来，不加的话柱状图是竖着的
  xlab("GO term") +
  ylab("Gene_Number") +
  labs(title = "The Most Enriched GO Terms")+
  theme_bw()
ggsave('./Figure6/ECEM_ESC_enrichGO.png', plotc, width = 12,height = 14)

# KEGG
genelist <- bitr(row.names(sig_dge.all), fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
genelist <- pull(genelist,ENTREZID)               
ekegg <- enrichKEGG(gene = genelist, organism = 'hsa')
p1 <- barplot(ekegg, showCategory=20)
p2 <- dotplot(ekegg, showCategory=20)
plotc = p1/p2
ggsave("./Figure6/ECEM_ESC_enrichKEGG.png", plot = plotc, width = 12, height = 14)










#############################################################################################################3
# Unkown2
sig_dge.all <- subset(Unkown2, p_val_adj<0.05&abs(avg_log2FC)>0.3) #所有差异基因
# # BP, CC和MF三种通路都一起富集
# ego_ALL <- enrichGO(gene          = row.names(sig_dge.all),
#                     #universe     = row.names(dge.celltype),
#                     OrgDb         = 'org.Hs.eg.db',
#                     keyType       = 'SYMBOL',
#                     ont           = "ALL",  #设置为ALL时BP, CC, MF都计算
#                     pAdjustMethod = "BH",
#                     pvalueCutoff  = 0.01,
#                     qvalueCutoff  = 0.05)
# ego_all <- data.frame(ego_ALL)
# 分别对BP, CC和MF进行富集
ego_CC <- enrichGO(gene          = sig_dge.all$gene,
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)
ego_MF <- enrichGO(gene          = sig_dge.all$gene,
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)
ego_BP <- enrichGO(gene          = sig_dge.all$gene,
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05) 
display_number = c(10, 10, 5)#这三个数字分别代表选取的BP、CC、MF的通路条数，这个自己设置就行了
ego_result_BP <- as.data.frame(ego_BP)[1:display_number[1], ]
ego_result_CC <- as.data.frame(ego_CC)[1:display_number[2], ]
ego_result_MF <- as.data.frame(ego_MF)[1:display_number[3], ]

##将以上我们摘取的部分通路重新组合成数据框
go_enrich_df <- data.frame(
  ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),                         
  Description=c(ego_result_BP$Description,ego_result_CC$Description,ego_result_MF$Description),
  GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
  type=factor(c(rep("biological process", display_number[1]), 
                rep("cellular component", display_number[2]),
                rep("molecular function", display_number[3])), 
              levels=c("biological process", "cellular component","molecular function" )))

##通路的名字太长了，我选取了通路的前五个单词作为通路的名字
# for(i in 1:nrow(go_enrich_df)){
#   description_splite=strsplit(go_enrich_df$Description[i],split = " ")
#   description_collapse=paste(description_splite[[1]][1:5],collapse = " ") #这里的5就是指5个单词的意思，可以自己更改
#   go_enrich_df$Description[i]=description_collapse
#   go_enrich_df$Description=gsub(pattern = "NA","",go_enrich_df$Description)
# }

##开始绘制GO柱状图
###横着的柱状图
go_enrich_df$type_order=factor(rev(as.integer(rownames(go_enrich_df))),labels=rev(go_enrich_df$Description))#这一步是必须的，为了让柱子按顺序显示，不至于很乱
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")#设定颜色

plotc<-ggplot(data=go_enrich_df, aes(x=type_order,y=GeneNumber, fill=type)) + #横纵轴取值
  geom_bar(stat="identity", width=0.8) + #柱状图的宽度，可以自己设置
  scale_fill_manual(values = COLS) + ###颜色
  coord_flip() + ##这一步是让柱状图横过来，不加的话柱状图是竖着的
  xlab("GO term") +
  ylab("Gene_Number") +
  labs(title = "The Most Enriched GO Terms")+
  theme_bw()
ggsave('./Figure6/Unkown2_enrichGO.png', plotc, width = 12,height = 14)

# KEGG
genelist <- bitr(sig_dge.all$gene, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
genelist <- pull(genelist,ENTREZID)               
ekegg <- enrichKEGG(gene = genelist, organism = 'hsa')
p1 <- barplot(ekegg, showCategory=20)
p2 <- dotplot(ekegg, showCategory=20)
plotc = p1/p2
ggsave("./Figure6/Unkown2_enrichKEGG.png", plot = plotc, width = 12, height = 14)




