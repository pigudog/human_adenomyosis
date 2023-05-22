# ##安装copykat
# library(devtools)
# install_github("navinlabcode/copykat")
# library("copykat") #安装成功可以简单测试这个包
# copykat.test <- copykat(rawmat=exp.rawdata, sam.name="test")
##安装InferCNV
# biglm
.libPaths(c('~/R/x86_64-pc-linux-gnu-library/4.2',
            '/usr/local/lib/R/library',
            '/home/data/refdir/Rlib'))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("infercnv")
#如果没有安装Seurat的，最好装一个Seurat
# conda install -c bioconda r-seurat 
##我安装的时候一个依赖包一直报错，用conda重新创建了一个环境，才安装成功
rm(list=ls())
library(infercnv)
library(Seurat)
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
load(file = 'AM_after_cluster.Rdata')
future::plan("multicore", workers = 32)
table(AM@meta.data$cell_type)
AM@meta.data$cell_type <- AM@active.ident
sce.all.list <- SplitObject(AM , split.by = "cell_type")
sce.all.list

for (i in names(sce.all.list)) {
  epi_mat = sce.all.list[[i]]@assays$RNA@counts
  epi_phe = sce.all.list[[i]]@meta.data
  sce=CreateSeuratObject(counts = epi_mat, 
                         meta.data = epi_phe )
  sce
  table(sce@meta.data$orig.ident) 
  save(sce,file = paste0(i,'.Rdata'))
}

rm(list=ls())
options(stringsAsFactors = F)
library(Seurat)
load(file = 'first_sce.Rdata')  # 这里面仅仅是上皮细胞哦
epiMat=as.data.frame(GetAssayData( sce , slot='counts',assay='RNA'))
load('../Rdata-celltype/B.Rdata') # 这里面仅仅是B淋巴 细胞哦
B_cellMat=as.data.frame(GetAssayData( sce, slot='counts',assay='RNA'))
load('../Rdata-celltype/T.Rdata') # 这里面仅仅是T淋巴 细胞哦
T_cellsMat=as.data.frame(GetAssayData( sce, slot='counts',assay='RNA'))

B_cellMat=B_cellMat[,sample(1:ncol(B_cellMat),800)]
T_cellsMat=T_cellsMat[,sample(1:ncol(T_cellsMat),800)]
dat=cbind(epiMat,B_cellMat,T_cellsMat)

groupinfo=data.frame(v1=colnames(dat),
                     v2=c(rep('epi',ncol(epiMat)),
                          rep('spike-B_cell',300),
                          rep('ref-B_cell',500),
                          rep('spike-T_cells',300),
                          rep('ref-T_cells',500)))



# 我反而觉得对所有细胞都做更简单明了-作为初步的看的方式
#第二个文件:单细胞RNA-seq表达量的原始矩阵
library(Seurat)
library(ggplot2)
library(infercnv)
library(future)
options(future.globals.maxSize= 8912896000)
sce.all.list <- SplitObject(AM , split.by = "orig.ident")
Pre_infercnv <- function(sce.all.list){
  future::plan("multicore", workers = 48)
  # Warning message:
  # In supportsMulticoreAndRStudio(...) :
  # [ONE-TIME WARNING] Forked processing ('multicore') is not supported when running R from RStudio because it is considered unstable. For more details, how to control forked processing or not, and how to silence this warning in future R sessions, see ?parallelly::supportsMulticore
  for (j in names(sce.all.list)) {
    cnv = sce.all.list[[j]]
    dfcount = as.data.frame(cnv@assays$RNA@counts)
    groupinfo= data.frame(cellId = colnames(dfcount),cellType= cnv@meta.data$cell_type )
    # 第三文件基因注释文件
    library(AnnoProbe)
    geneInfor=annoGene(rownames(dfcount),"SYMBOL",'human')
    geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
    geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
    ## 这里可以去除性染色体
    # 也可以把染色体排序方式改变
    dfcount =dfcount [rownames(dfcount ) %in% geneInfor[,1],]
    dfcount =dfcount [match( geneInfor[,1], rownames(dfcount) ),] 
    write.table(dfcount ,file ='./inferCNV/expFile.txt',sep = '\t',quote = F)
    write.table(groupinfo,file = './inferCNV/metaFiles.txt',sep = '\t',quote = F,col.names = F,row.names = F)
    write.table(geneInfor,file ='./inferCNV/geneFile.txt',sep = '\t',quote = F,col.names = F,row.names = F)
    
    infercnv_obj = CreateInfercnvObject(raw_counts_matrix="./inferCNV/expFile.txt",
                                    annotations_file="./inferCNV/metaFiles.txt",
                                    delim="\t",
                                    gene_order_file= "./inferCNV/geneFile.txt",
                                    ref_group_names=NULL) # 如果有正常细胞的话，把正常细胞的分组填进去

    infercnv_obj = infercnv::run(infercnv_obj,
                                 cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                                 out_dir=paste0(j,'_infercnv/'), 
                                 cluster_by_groups=TRUE,  # 是否根据细胞注释文件的分组对肿瘤细胞进行
                                 denoise=TRUE, #去噪
                                 HMM=TRUE, #是否基于HMM预测CNV,选择F会加快运算速度
                                 output_format = "pdf")
    
    save(infercnv_obj,file=paste0(j,'_infercnv.Rdata'))
    print(paste0(j,'_finished'))
  }
}
Pre_infercnv(sce.all.list=sce.all.list)



dfcount = as.data.frame(AM@assays$RNA@counts)
# 第二个文件:注释文件，记录肿瘤和正常细胞
groupinfo= data.frame(cellId = colnames(dfcount),cellType= AM@meta.data$cell_type )
# 第三文件基因注释文件
library(AnnoProbe)
geneInfor=annoGene(rownames(dfcount),"SYMBOL",'human')
geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
## 这里可以去除性染色体
# 也可以把染色体排序方式改变
dfcount =dfcount [rownames(dfcount ) %in% geneInfor[,1],]
dfcount =dfcount [match( geneInfor[,1], rownames(dfcount) ),] 
# 输出
write.table(dfcount ,file ='./inferCNV/expFile.txt',sep = '\t',quote = F)
write.table(groupinfo,file = './inferCNV/metaFiles.txt',sep = '\t',quote = F,col.names = F,row.names = F)
write.table(geneInfor,file ='./inferCNV/geneFile.txt',sep = '\t',quote = F,col.names = F,row.names = F)

##
library(Seurat)
library(ggplot2)
library(infercnv)
future::plan("multicore", workers = 48)
infercnv_obj = CreateInfercnvObject(raw_counts_matrix="./inferCNV/expFile.txt",
                                    annotations_file="./inferCNV/metaFiles.txt",
                                    delim="\t",
                                    gene_order_file= "./inferCNV/geneFile.txt",
                                    ref_group_names=NULL) # 如果有正常细胞的话，把正常细胞的分组填进去

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir='infercnv_out/', 
                             cluster_by_groups=TRUE,  # 是否根据细胞注释文件的分组对肿瘤细胞进行
                             denoise=F, #去噪
                             HMM=F, #是否基于HMM预测CNV,选择F会加快运算速度
                             output_format = "pdf")
# cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
# infercnv_obj2 = infercnv::run(infercnv_obj,
#                               cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
#                               out_dir= "infercnv_output",  # dir is auto-created for storing outputs
#                               cluster_by_groups=F,   # cluster
#                               hclust_method="ward.D2", plot_steps=F)
#也就可以对热图进行改造  换颜色（用inferCNV的官方的画图函数）
save(infercnv_obj,file="infercnv.Rdata")
infercnv::plot_cnv(infercnv_obj, #上两步得到的infercnv对象
                   plot_chr_scale = T, #画染色体全长，默认只画出（分析用到的）基因
                   output_filename = "better_plot",output_format = "pdf", #保存为pdf文件
                   custom_color_pal =  color.palette(c("#8DD3C7","white","#BC80BD"), c(2, 2))) #改颜色

#第一步：提取inferCNV的结果   
obs <- read.table("infercnv_out/infercnv.observations.txt",header = T,check.names = F)     
#首先这个obs是每一个基因在每一个细胞的拷贝信息,相当于该基因的拷贝量
if(T){              #可以通过定义obs的元素数值来让差异变大   在后面画图就能够更大的差异   也可以不运行，更真实体现拷贝数（但是可能没有啥差异）
  max(obs)          #根据最大最小值来定义
  min(obs)     
  obs[obs > 0.6 & obs < 0.7] <- -2   #把0.6-0.7定义为数值2  后面依此类推
  obs[obs >= 0.7 & obs < 0.8] <- -1
  obs[obs >= 0.8 & obs < 1.0] <- 0
  obs[obs >= 1.0 & obs <= 1.1] <- 1
  obs[obs > 1.1 & obs <= 1.3] <- 2
  obs[obs > 1.3] <- 2
}

score=as.data.frame(colSums(obs))   #把obs的每一个基因拷贝量加起来，就是这个细胞的总拷贝数obs

#提取meta信息  celltype,seurat_clusters,orig.ident
meta <- subset(AM@meta.data,select = c("cell_type","seurat_clusters","orig.ident"))  #提取该细胞的其他的meta信息

#将meta信息添加给score
meta <- rownames_to_column(groupinfo)
score <- rownames_to_column(score)
score <- merge(score,meta,by.x = "rowname",by.y = "rowname")   
#这里会可能损失一些细胞   为什么呢，因为在前面infer时，有一些细胞无法推断，会被删掉，但是总体上问题不大

#我们现在可以用ggplot画图了     但是直接这样出图很丑   因为他会根据你的Y轴的大小  所以建议定义y轴范围
ggboxplot(score,"orig.ident","colSums(obs)",fill="orig.ident")   #也可以将orig.ident换成celltype
#+scale_y_continuous(limits = c(0, 1000))


#################################################################################################
usethis::create_github_token()#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
usethis::edit_r_environ()
install_github("navinlabcode/copykat")
library(devtools)
library(Seurat)
library(copykat)
library(tidyverse)
library(Seurat)
library(ggplot2)
library(infercnv)
library(future)
options(future.globals.maxSize= 8912896000)
load(file = 'AM_after_cluster.Rdata')
AM@meta.data$cell_type <- AM@active.ident
sce.all.list <- SplitObject(AM , split.by = "orig.ident")
Pre_copykat <- function(sce.all.list){
  future::plan("multicore", workers = 48)
  # Warning message:
  # In supportsMulticoreAndRStudio(...) :
  # [ONE-TIME WARNING] Forked processing ('multicore') is not supported when running R from RStudio because it is considered unstable. For more details, how to control forked processing or not, and how to silence this warning in future R sessions, see ?parallelly::supportsMulticore
  for (j in names(sce.all.list)) {
    cnv = sce.all.list[[j]]
    exp.rawdata <- as.matrix(cnv@assays$RNA@counts)
    library(copykat)
    copykat.test <- copykat(rawmat=exp.rawdata, id.type="S", ngene.chr=5, win.size=25, KS.cut=0.1, 
                            sam.name=j, distance="euclidean", 
                            norm.cell.names="",output.seg="FLASE", plot.genes="TRUE", 
                            genome="hg20",n.cores=32)
    # ngene.chr参数是过滤细胞的一个标准，它要求被用来做CNV预测的细胞，一个染色体上至少有5个基因。
    # sam.name定义样本名称 (sample name)，会给出来的文件加前缀
    save(copykat.test ,file = paste0(j,'_copykat.Rdata'))
    pred.test <- data.frame(copykat.test$prediction)
    library(ggsci)
    library(Cairo)
    color = c(pal_d3("category20")(20),
              pal_d3("category20b")(20),
              pal_d3("category20c")(20),
              pal_d3("category10")(10))
    cnv@meta.data$copykat.pred <- pred.test$copykat.pred
    p1 <- DimPlot(cnv, reduction = 'umap', 
                  group.by = "cell_type",
                  cols = color , #设置颜色
                  repel = T,  #label不重叠
                  label.size = 5, #调整label的大小
                  label = TRUE,  #是否展示label
                  pt.size = 0.5);p1
    p2 <- DimPlot(cnv, group.by = "copykat.pred") + scale_color_manual(values = c("red","blue", "grey"))
    pc <- p1 + p2
    ggsave(paste0("./copykat/",j,'_pred_mallignant.pdf'), pc, width = 12, height = 6)
    
    # pred.test <- data.frame(copykat.test$prediction)
    # CNA.test <- data.frame(copykat.test$CNAmat)
    # pred.test <- pred.test[-which(pred.test$copykat.pred=="not.defined"),]  ##remove undefined cells
    # my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu")))(n = 999)
    # 
    # chr <- as.numeric(CNA.test$chrom) %% 2+1
    # rbPal1 <- colorRampPalette(c('black','grey'))
    # CHR <- rbPal1(2)[as.numeric(chr)]
    # chr1 <- cbind(CHR,CHR)
    # 
    # rbPal5 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1])
    # com.preN <- pred.test$copykat.pred
    # pred <- rbPal5(2)[as.numeric(factor(com.preN))]
    # 
    # cells <- rbind(pred,pred)
    # col_breaks = c(seq(-1,-0.4,length=50),seq(-0.4,-0.2,length=150),seq(-0.2,0.2,length=600),seq(0.2,0.4,length=150),seq(0.4, 1,length=50))
    # 
    # p_heatmap<-heatmap.3(t(CNA.test[,4:ncol(CNA.test)]),dendrogram="r", distfun = function(x) parallelDist::parDist(x,threads =4, method = "euclidean"), hclustfun = function(x) hclust(x, method="ward.D2"),
    #                      ColSideColors=chr1,RowSideColors=cells,Colv=NA, Rowv=TRUE,
    #                      notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
    #                      keysize=1, density.info="none", trace="none",
    #                      cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
    #                      symm=F,symkey=F,symbreaks=T,cex=1, cex.main=4, margins=c(10,10))+
    #   legend("topright", paste("pred.",names(table(com.preN)),sep=""), pch=15,col=RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1], cex=0.6, bty="n")
    # ggsave(paste0("./copykat/",j,'_pred_mallignant_heatmap.pdf'),p_heatmap , width = 16, height = 14)
    print(paste0(j,'_finished'))
  }
}
Pre_copykat(sce.all.list=sce.all.list)

















load(file = 'AM_after_cluster.Rdata')
table(AM@meta.data$cell_type)
table(AM$orig.ident)
AM@meta.data$cell_type <- AM@active.ident
AM_EM <- AM[,AM$orig.ident %in% c("EM")]
exp.rawdata <- as.matrix(AM_EM@assays$RNA@counts)
write.table(exp.rawdata, file="exp.rawdata.txt", sep="\t", quote = FALSE, row.names = TRUE)

library(copykat)
copykat.test <- copykat(rawmat=exp.rawdata, id.type="S", ngene.chr=5, win.size=25, KS.cut=0.1, 
                        sam.name="AM_EM", distance="euclidean", 
                          norm.cell.names="",output.seg="FLASE", plot.genes="TRUE", 
                        genome="hg20",n.cores=32)
# ngene.chr参数是过滤细胞的一个标准，它要求被用来做CNV预测的细胞，一个染色体上至少有5个基因。
# sam.name定义样本名称 (sample name)，会给出来的文件加前缀
save(copykat.test ,file = 'copykat_AM_EM.Rdata')
load("copykat_AM_EM.Rdata")
pred.test <- data.frame(copykat.test$prediction)

library(ggsci)
library(Cairo)
color = c(pal_d3("category20")(20),
          pal_d3("category20b")(20),
          pal_d3("category20c")(20),
          pal_d3("category10")(10))
AM_EM@meta.data$copykat.pred <- pred.test$copykat.pred
p1 <- DimPlot(AM_EM, reduction = 'umap', 
              group.by = "cell_type",
              cols = color , #设置颜色
              repel = T,  #label不重叠
              label.size = 5, #调整label的大小
              label = TRUE,  #是否展示label
              pt.size = 0.5);p1
p2 <- DimPlot(AM_EM, group.by = "copykat.pred") + scale_color_manual(values = c("red","blue", "white"))
p2
pc <- p1 + p2
ggsave("./AM/AM_EM_pred_mallignant.pdf", pc, width = 12, height = 6)







pred.test <- data.frame(copykat.test$prediction)
CNA.test <- data.frame(copykat.test$CNAmat)
pred.test <- pred.test[-which(pred.test$copykat.pred=="not.defined"),]  ##remove undefined cells
my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu")))(n = 999)

chr <- as.numeric(CNA.test$chrom) %% 2+1
rbPal1 <- colorRampPalette(c('black','grey'))
CHR <- rbPal1(2)[as.numeric(chr)]
chr1 <- cbind(CHR,CHR)

rbPal5 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1])
com.preN <- pred.test$copykat.pred
pred <- rbPal5(2)[as.numeric(factor(com.preN))]

cells <- rbind(pred,pred)
col_breaks = c(seq(-1,-0.4,length=50),seq(-0.4,-0.2,length=150),seq(-0.2,0.2,length=600),seq(0.2,0.4,length=150),seq(0.4, 1,length=50))

p_heatmap<-heatmap.3(t(CNA.test[,4:ncol(CNA.test)]),dendrogram="r", distfun = function(x) parallelDist::parDist(x,threads =4, method = "euclidean"), hclustfun = function(x) hclust(x, method="ward.D2"),
          ColSideColors=chr1,RowSideColors=cells,Colv=NA, Rowv=TRUE,
          notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
          keysize=1, density.info="none", trace="none",
          cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
          symm=F,symkey=F,symbreaks=T,cex=1, cex.main=4, margins=c(10,10))+
  legend("topright", paste("pred.",names(table(com.preN)),sep=""), pch=15,col=RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1], cex=0.6, bty="n")
ggsave("./AM/AM_EM_pred_mallignant_heatmap.pdf",p_heatmap , width = 16, height = 14)

