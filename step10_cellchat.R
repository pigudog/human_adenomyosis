options(timeout=10000)
# install.packages("pak", repos = sprintf(
#   "https://r-lib.github.io/p/pak/stable/%s/%s/%s",
#   .Platform$pkgType,
#   R.Version()$os,
#   R.Version()$arch
# ))
install.packages("glue")
install.packages("pak", INSTALL_opts = '--no-lock')
# stream <- "stable"
# install.packages("pak", repos = sprintf(
#   "https://r-lib.github.io/p/pak/%s/%s/%s/%s",
#   stream,
#   .Platform$pkgType,
#   R.Version()$os,
#   R.Version()$arch
# ))
pkgs <- c("BiocNeighbors", "ComplexHeatmap", "circlize", "NMF")
install.packages("pak")
pak::pkg_install(pkgs, upgrade = FALSE, dependencies = TRUE)
pak::pkg_install("sqjin/CellChat")
## R环境安装NMF和ComplexHeatmap包
devtools::install_github("renozao/NMF@devel")
devtools::install_github("jokergoo/ComplexHeatmap")
devtools::install_github("sqjin/CellChat")
reticulate::py_install(packages ='umap-learn')
library(CellChat)
library(patchwork)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
options(stringsAsFactors = FALSE)
rm(list=ls())
load(file = 'AM_after_cluster.Rdata')
table(AM@meta.data$cell_type)
AM@meta.data$cell_type <- AM@active.ident
AM <- NormalizeData(object = AM)
cellchat <- createCellChat(object = AM, group.by = "ident", assay = "RNA")
levels(cellchat@idents) # show factor levels of the cell labels

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB
# set the used database in the object
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat) 
# subset the expression data of signaling genes for saving computation cost
# future::plan("multiprocess", workers = 4) # do parallel
#Error: Strategy 'multiprocess' is defunct in future (>= 1.32.0) [2023-03-06]. 
# Instead, explicitly specify either 'multisession' (recommended) or 'multicore'.
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
pdf("./AM/cellchat1.pdf")
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
pdf("./AM/cellchat2.pdf")
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

cellchat@netP$pathways
# WNT
pathways.show <- c("WNT") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
pdf("./AM/cellchat_WNT.pdf")
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
dev.off()

# Circle plot
pdf("./AM/cellchat_WNT_cirvle.pdf")
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
dev.off()

# Chord diagram
pdf("./AM/cellchat_WNT_chord.pdf")
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.
dev.off()

# Heatmap
pdf("./AM/cellchat_WNT_heatmap.pdf")
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object
dev.off()

# PTN
pathways.show <- c("PTN") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
pdf("./AM/cellchat_PTN.pdf")
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
dev.off()

# Heatmap
pdf("./AM/cellchat_PTN_heatmap.pdf")
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object
dev.off()

# VISFATIN
pathways.show <- c("VISFATIN") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
pdf("./AM/cellchat_VISFATIN.pdf")
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
dev.off()

# Heatmap
pdf("./AM/cellchat_VISFATIN_heatmap.pdf")
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object
dev.off()

# egf
pathways.show <- c("EGF") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
pdf("./AM/cellchat_EGF.pdf")
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
dev.off()

# Heatmap
pdf("./AM/cellchat_EGF_heatmap.pdf")
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object
dev.off()

# MIF
pathways.show <- c("MIF") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
pdf("./AM/cellchat_MIF.pdf")
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
dev.off()

# Heatmap
pdf("./AM/cellchat_MIF_heatmap.pdf")
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object
dev.off()


# FGF
pathways.show <- c("FGF") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
pdf("./AM/cellchat_FGF.pdf")
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
dev.off()

# Heatmap
pdf("./AM/cellchat_FGF_heatmap.pdf")
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object
dev.off()

# TGFb
pathways.show <- c("TGFb") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
pdf("./AM/cellchat_TGFb.pdf")
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
dev.off()

# Heatmap
pdf("./AM/cellchat_TGFb_heatmap.pdf")
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object
dev.off()

# Heatmap
pdf("./AM/cellchat_FGF_heatmap.pdf")
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object
dev.off()

# PDGF
pathways.show <- c("PDGF") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
pdf("./AM/cellchat_PDGF.pdf")
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
dev.off()

# Heatmap
pdf("./AM/cellchat_PDGF_heatmap.pdf")
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object
dev.off()


# CCL
pathways.show <- c("CCL") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
pdf("./AM/cellchat_CCL.pdf")
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
dev.off()

# Circle plot
pdf("./AM/cellchat_CCL_cirvle.pdf")
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
dev.off()

# Chord diagram
pdf("./AM/cellchat_CCL_chord.pdf")
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.
dev.off()

# Heatmap
pdf("./AM/cellchat_CCL_heatmap.pdf")
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object
dev.off()

# CXCL
pathways.show <- c("CXCL") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
pdf("./AM/cellchat_CXCL.pdf")
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
dev.off()

# Circle plot
pdf("./AM/cellchat_CXCL_cirvle.pdf")
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
dev.off()

# Chord diagram
pdf("./AM/cellchat_CXCL_chord.pdf")
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.
dev.off()

# Heatmap
pdf("./AM/cellchat_CXCL_heatmap.pdf")
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object
dev.off()

netAnalysis_contribution(cellchat, signaling = pathways.show)

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
# the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
dev.off()

pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)


########################################################################
# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = seq(1,4)
setwd("./cellchatAM2/")
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}
setwd("../")



saveRDS(cellchat, file = "cellchat.rds")



############################################################################333333
# compare
remotes::install_github("HenrikBengtsson/future", ref="develop")
library(future)
library(CellChat)
library(patchwork)
library(Seurat)
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
options(stringsAsFactors = FALSE)
rm(list=ls())
# 这儿给了四个函数
# 其中ps是画图用的
# Pre_cellchat是将数据转化成cellchat
# draw——cellchat是分别对样本去绘图
# compare 是两两比较，输入样本的名称即可
# lift, template_name是你选择要参考的那个样本的名称
ps <- function(filename, plot = FALSE, w = 9, h = 6) {
  if (is.object(plot)) {
    print(plot)
  }
  plot <- recordPlot()
  pdf(file = paste0(filename,".pdf"), onefile = T, width = w, height = h)
  replayPlot(plot)
  dev.off()
}
Pre_cellchat <- function(sce.all.list){
  CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
  CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
  plan("multicore", workers = 32)
  # Warning message:
  # In supportsMulticoreAndRStudio(...) :
  # [ONE-TIME WARNING] Forked processing ('multicore') is not supported when running R from RStudio because it is considered unstable. For more details, how to control forked processing or not, and how to silence this warning in future R sessions, see ?parallelly::supportsMulticore
  for (i in names(sce.all.list)) {
    cellchat = sce.all.list[[i]]
    cellchat <- createCellChat(object = cellchat, group.by = "ident", assay = "RNA")
    cellchat@DB <- CellChatDB.use
    cellchat <- subsetData(cellchat) 
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- projectData(cellchat, PPI.human)
    cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
    # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat) 
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
    save(cellchat,file = paste0(i,'cellchat.Rdata'))
  }
}
dir_create_list <- function(object.list) {
  for(i in  names(object.list)){
    if (!file.exists(paste0(i,'_cellchat'))){
      dir.create(paste0(i,'_cellchat'))
    }
  }
}
draw_cellchat <- function(object.list){
  dir_create_list(object.list=object.list)
  for(j in  names(object.list)){
    cellchat <- object.list[[j]]
    groupSize <- as.numeric(table(cellchat@idents))
    par(mfrow = c(1,2), xpd=TRUE)
    netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
    netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
    ps(paste0(j,'_cellchat/cellchat1.pdf'))
    
    mat <- cellchat@net$weight
    par(mfrow = c(3,4), xpd=TRUE)
    for (i in 1:nrow(mat)) {
      mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
      mat2[i, ] <- mat[i, ]
      netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
    }
    pdf(paste0(j,'_cellchat/cellchat2.pdf'))
    
    pathways.show.all <- cellchat@netP$pathways
    vertex.receiver = seq(1,4)
    setwd(paste0(j,'_cellchat/'))
    for (i in 1:length(pathways.show.all)) {
      # Visualize communication network associated with both signaling pathway and individual L-R pairs
      netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
      # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
      gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
      ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
    }
    setwd("../")
  }
}
compare_cellchat <- function(object.list,seq=c("CTRL","EC")){
  
  object.list1 <- list(object.list[[seq[1]]], object.list[[seq[2]]])
  names(object.list1) <- seq
  cellchat <- mergeCellChat(object.list1, add.names = names(object.list1))
  gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
  gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
  gg1 + gg2
  ps(paste0("./comparison/Strength_of_interaction_",seq[1],seq[2]))
 
  
  par(mfrow = c(1,2), xpd=TRUE)
  netVisual_diffInteraction(cellchat, weight.scale = T)
  netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
  ps(paste0("./comparison/circle_of_interaction_",seq[1],seq[2]))
  
  
  gg1 <- netVisual_heatmap(cellchat)
  #> Do heatmap based on a merged object
  gg2 <- netVisual_heatmap(cellchat, measure = "weight")
  #> Do heatmap based on a merged object
  gg1 + gg2
  ps(paste0("./comparison/heatmap_of_interaction_",seq[1],seq[2]))
  
  
  gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
  gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
  gg1 + gg2
  ps(paste0("./comparison/rankNet_of_interaction_",seq[1],seq[2]))
  
}
# Define the cell labels to lift up
liftcellchat <- function(object.list,template_name){
  group.new = levels(object.list[[template_name]]@idents)
  for(j in  names(object.list)){
    object.list[[j]] <- liftCellChat(object.list[[j]], group.new)
  }
  
}

data.dir <- './comparison'
dir.create(data.dir)
load(file = 'AM_after_cluster.Rdata')
AM@meta.data$cell_type <- AM@active.ident
# 按照需要对比的组进行区分
sce.all.list <- SplitObject(AM, split.by = "orig.ident")
Pre_cellchat(sce.all.list= sce.all.list)

load("CTRLcellchat.Rdata")
cellchat_CTRL <- cellchat
load("EMcellchat.Rdata")
cellchat_EM <- cellchat
load("ECcellchat.Rdata")
cellchat_EC <- cellchat
object.list <- list(CTRL = cellchat_CTRL, EM = cellchat_EM, EC =cellchat_EC)

draw_cellchat(object.list = object.list)
liftcellchat(object.list=object.list,template_name="EM")



compare_cellchat(object.list,seq=c("CTRL","EM"))
compare_cellchat(object.list,seq=c("CTRL","EC"))
compare_cellchat(object.list,seq=c("EM","EC"))






