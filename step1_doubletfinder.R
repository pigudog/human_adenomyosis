library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(stringr)
library(harmony)
library(DoubletFinder)
rm(list=ls())
options(stringsAsFactors = F)
options (warn = -1)

rm_doublet <- function(name=NULL,input=NULL,dim.usage=30,auto="false") {
  ## get the path of the matric 
  ## The matrix passed in here is in the 10X format
  inpath <- list.files(path = input,pattern = name,full.names = T)
  if (auto=="true") {
    inpath <- paste0(inpath,"/","04.Matrix/")#这是另一种矩阵路径结构
  }
  
  ## before doubeltFinder, you need cluster, just choose the default parameters
  ## This section can be written more succinctly using the pipe character %>%
  EC <- Read10X(inpath,gene.column=2)
  EC <- CreateSeuratObject(EC, min.cells = 3, min.features = 200)
  EC <- SCTransform(EC)
  EC <- RunPCA(EC,verbose=F)
  EC <- RunUMAP(EC, dims = 1:dim.usage) 
  EC <- FindNeighbors(EC, dims = 1:dim.usage) %>% FindClusters(resolution = 0.3) 
  ##DoubletFinder去双胞的标准流程，封装成一个函数
  Find_doublet <- function(data){
    ### pK Identification (no ground-truth)
    sweep.res.list <- paramSweep_v3(data, PCs = 1:dim.usage, sct = T)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
    ### Homotypic Doublet Proportion Estimate
    DoubletRate = ncol(pbmc)*8*1e-6
    homotypic.prop <- modelHomotypic(data$seurat_clusters)   # 最好提供celltype
    nExp_poi <- round(DoubletRate*ncol(data))
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    
    ### Run DoubletFinder with varying classification stringencies
    data <- doubletFinder_v3(data, PCs = 1:dim.usage, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = T)
    colnames(data@meta.data)[ncol(data@meta.data)] = "doublet_info"
    return(data)
  }
  
  ##调用上面写好的函数，返回的是一个Seurat对象，meta.data信息里会有双胞信息，需要自己手动删除
  EC<-Find_doublet(EC)
  EC<-subset(EC,subset=doublet_info=="Singlet")
  EC@meta.data$library = name #顺便打上这个样本的label
  
  ##生成的Seurat对象有个问题，会在meta.data里多了很多pANN_开头的列，需要手动删除
  c <- grep("pANN_",colnames(EC@meta.data))
  EC@meta.data <- EC@meta.data[,-c]
  ##输出此样本的细胞数
  print(paste0(name," cells: ", length(EC@meta.data$orig.ident)," is read!"," Time:",format(Sys.time(), "%Y%m%d %X"),sep = " "))
  return(EC)
}

# when the number of sample <20 (we just choose this)
folders=list.files('/home/data/t070421/adenomyosis/AM/filtered/')
write.table (folders, file ="/home/data/t070421/adenomyosis/scrna/infile.txt", sep ="", row.names =FALSE, col.names =FALSE, quote =FALSE)

## set parameters
infile = "/home/data/t070421/adenomyosis/scrna/infile.txt"
input = "/home/data/t070421/adenomyosis/AM/filtered/"
samples_name = c("CTRL", "EC" ,  "EM")
merge_seob <- function(infile,input,batch=NULL){
  filelist <- readLines(infile)
  seob_list <- lapply(filelist, rm_doublet, input=input)
  n <- length(seob_list)
  for(i in 1:n){
    #给细胞barcode加个前缀，防止合并后barcode重名
    seob_list[[i]] <- RenameCells(seob_list[[i]], add.cell.id = samples_name[i])   
    
  }
  names(seob_list) <- samples_name
  all <- merge(seob_list[[1]], seob_list[c(2:n)])
  if (!is.null(batch)) {
    all$batch <- batch
  }
  return(all)
}
scRNA <- merge_seob(infile=infile,input=input)
table(scRNA$library)

scRNA@meta.data[["orig.ident"]]<-scRNA@meta.data[["library"]]
level0 = as.factor(scRNA@meta.data[["orig.ident"]])
table(level0)
scRNA <- AddMetaData(object=scRNA,
                     metadata=level0,
                     col.name="level0")
tail(scRNA@meta.data)

# mitochondrion, ribosome、erythrocyte
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
scRNA[["percent.rb"]] <- PercentageFeatureSet(scRNA, pattern = "^RP[SL]")
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB.genes <- CaseMatch(HB.genes, rownames(scRNA))
scRNA[["percent.HB"]]<-PercentageFeatureSet(scRNA, features=HB.genes)

### QC violin picture
# possible theme
theme.set2 = theme(axis.title.x=element_blank())
# the elements of picture
plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb", "percent.HB")
group = "orig.ident"
# vlnplot before qc
plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(scRNA, group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=2)    
dir.create("AM")
ggsave("AM/vlnplot_before_qc.pdf", plot = violin, width = 12, height = 8)


minGene=500
maxGene=10000
maxUMI=15000
pctMT=20
pctHB=1

### vlnplot after QC
scRNA <- subset(scRNA, subset = nCount_RNA < maxUMI & nFeature_RNA > minGene & nFeature_RNA < maxGene & percent.mt < pctMT & percent.HB <pctHB)
plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(scRNA, group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=2)     
ggsave("AM/vlnplot_after_qc.pdf", plot = violin, width = 14, height = 8) 

# select useful metadata
cellinfo <- subset(scRNA@meta.data, select = c("orig.ident", "percent.mt", "percent.rb", "percent.HB","level0"))
AM <- CreateSeuratObject(scRNA@assays$RNA@counts, meta.data = cellinfo)


save(AM,file = 'AM.Rdata')
