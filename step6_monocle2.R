library(monocle)
detach("package:monocle3")
library(monocle3)
# 
pkgs <- c("url::https://github.com/cole-trapnell-lab/monocle-release/files/10134172/monocle_2.26.0.tar.gz")
pak::pkg_install(pkgs)
rm(list=ls())
load(file = 'AM_after_cluster.Rdata')
AM@meta.data$cell_type <- AM@active.ident




# AM1, unkown1 to (Mac and Mast)
AM1 <- AM[,AM$cell_type %in% c("Mast","Mac","Unkown1")]
pbmc <- AM1 #导入注释好的seurat对象（已注释）
# AM2, 第二个monocle：fibro-like，unkonwn2，SMA
AM2 <- AM[,AM$cell_type %in% c("Unkown2", "Fibro-like","SMC")]
pbmc <-AM2
# AM3, 第二个monocle：Epi，unkonwn1，En
AM3 <- AM[,AM$cell_type %in% c("Unkown1", "Endo","Epi")]
pbmc <-AM3

##提取表型信息--细胞信息(建议载入细胞的聚类或者细胞类型鉴定信息、实验条件等信息)
expr_matrix <- as(as.matrix(pbmc@assays$RNA@counts), 'sparseMatrix')
##提取表型信息到p_data(phenotype_data)里面 
p_data <- pbmc@meta.data 
# p_data$celltype <- pbmc@active.ident  ##整合每个细胞的细胞鉴定信息到p_data里面。如果已经添加则不必重复添加
##提取基因信息 如生物类型、gc含量等
f_data <- data.frame(gene_short_name = row.names(expr_matrix),row.names = row.names(expr_matrix))
##expr_matrix的行数与f_data的行数相同(gene number), expr_matrix的列数与p_data的行数相同(cell number)

#构建CDS对象
pd <- new('AnnotatedDataFrame', data = p_data) 
fd <- new('AnnotatedDataFrame', data = f_data)

#将p_data和f_data从data.frame转换AnnotatedDataFrame对象。
# 如果我们想要从Seurat对象或SCESet中导入所有的插槽，
# 我们可以设置参数'import_all'为TRUE。#(默认为FALSE或只保留最小数据集)
cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())


# size facotr帮助我们标准化细胞之间的mRNA的差异。
# 离散度值可以帮助我们进行后续的差异分析。
# （类似于seurat的数据归一化处理）
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
# 下面是monocle对新构建的CellDataSet 对象的标准操作， 注意estimateDispersions这步的时间和电脑配置密切相关，
# 甚至如果电脑内存不够，还会报错！
# 与seurat把标准化后的表达矩阵保存在对象中不同，
# monocle只保存一些中间结果在对象中，需要用时再用这些中间结果转化。
# 经过上面三个函数的计算，mycds对象中多了SizeFactors、Dipersions、num_cells_expressed和num_genes_expressed等信息。



# 因为Seurat已经完成细胞过滤，此步可省略
# 但由于Seurat是通过基因表达量来对细胞进行过滤。
# 因此，我们可以通过下面的代码用表达某基因的细胞的数目对基因进行过滤，
# 从而得到后续操作需要的基因。
cds <- detectGenes(cds, min_expr = 0.1) #这一操作会在fData(cds)中添加一列num_cells_expressed
print(head(fData(cds)))#此时有13714个基因
expressed_genes <- row.names(subset(fData(cds),
                                    num_cells_expressed >= 10)) #过滤掉在小于10个细胞中表达的基因，还剩11095个基因。

#The ordering workflow
# Step 1: choosing genes that define progress
# Step 2: reducing the dimensionality of the data
# Step 3: ordering the cells in pseudotime

## Step 1: 选择定义过程的基因
#Monocle官网教程提供了4个选择方法：
# 选择发育差异表达基因
# 选择clusters差异表达基因
# 选择离散程度高的基因
# 自定义发育marker基因
# ##使用seurat选择的高变基因⚠️
# express_genes <- VariableFeatures(pbmc)
# cds <- setOrderingFilter(cds, express_genes)
# plot_ordering_genes(cds)
# ##使用clusters差异表达基因
# deg.cluster <- FindAllMarkers(pbmc)
# express_genes <- subset(deg.cluster,p_val_adj<0.05)$gene
# cds <- setOrderingFilter(cds, express_genes)
# plot_ordering_genes(cds)
# ##使用monocle选择的高变基因⚠️
# disp_table <- dispersionTable(cds)
# disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
# cds <- setOrderingFilter(cds, disp.genes)
# plot_ordering_genes(cds)


#理想情况下，我们希望尽可能少地使用正在研究的系统生物学的先验知识。
# 我们希望从数据中发现重要的排序基因，而不是依赖于文献和教科书，
# 因为这可能会在排序中引入偏见。我们将从一种更简单的方法开始，
# 但是我们通常推荐一种更复杂的方法，称为“dpFeature”。



#这一步输入的expressed_genes来自于步骤4。
#⚠️⚠️后续分析使用的是该方法
#也可输入seurat筛选出的高变基因：
library(Seurat)
expressed_genes <- VariableFeatures(pbmc) 
diff <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr="~cell_type",cores=1) 
#~后面是表示对谁做差异分析的变量，理论上可以为p_data的任意列名
head(diff)

##差异表达基因作为轨迹构建的基因,
# 差异基因的选择标准是qval<0.01,decreasing=F表示按数值增加排序
deg <- subset(diff, qval < 0.01) #选出2829个基因
deg <- deg[order(deg$qval,decreasing=F),]
head(deg)

##差异基因的结果文件保存
write.table(deg,file="train.monocle.DEG.csv",col.names=T,row.names=F,sep="\t",quote=F)
write.table(deg,file="./AM/AM1train.monocle.DEG.csv",col.names=T,row.names=F,sep="\t",quote=F)
write.table(deg,file="./AM/AM2train.monocle.DEG.csv",col.names=T,row.names=F,sep="\t",quote=F)
write.table(deg,file="./AM/AM3train.monocle.DEG.csv",col.names=T,row.names=F,sep="\t",quote=F)
## 轨迹构建基因可视化
# 选择的用于排序的基因数目一般在2000左右比较合适
# gene数太多的话也可以选择top基因
# ordergene <- row.names(deg)[order(deg$qval)][1:400]
ordergene <- rownames(deg) 
cds <- setOrderingFilter(cds, ordergene)  
#这一步是很重要的，在我们得到想要的基因列表后，
# 我们需要使用setOrderingFilter将它嵌入cds对象，
# 后续的一系列操作都要依赖于这个list。
# setOrderingFilter之后，
# 这些基因被储存在cds@featureData@data[["use_for_ordering"]]，
# 可以通过table(cds@featureData@data[["use_for_ordering"]])查看
pdf("./AM/AM1train.ordergenes.pdf")
plot_ordering_genes(cds)
dev.off()
pdf("./AM/AM2train.ordergenes.pdf")
plot_ordering_genes(cds)
dev.off()
pdf("./AM/AM3train.ordergenes.pdf")
plot_ordering_genes(cds)
dev.off()
#出的图黑色的点表示用来构建轨迹的差异基因，
# 灰色表示背景基因。
# 红色的线是根据第2步计算的基因表达大小和离散度分布的趋势(可以看到，找到的基因属于离散度比较高的基因)


# Step 2: 降维
# 一旦细胞有序排列，我们就可以在降维空间中可视化轨迹。
# 所以首先选择用于细胞排序的基因，
# 然后使用反向图嵌入(DDRTree)算法对数据进行降维
cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')


# Step 3: 拟时间轴轨迹构建和在拟时间内排列细胞
# 将表达数据投射到更低的维度空间，
# 通过机器学习描述细胞如何从一种状态过渡到另一种状态的轨迹。\
# 假设轨迹具有树状结构，一端是“根”，另一端是“叶”。
# 尽可能地将最佳树与数据匹配起来。
# 这项任务被称为“歧管学习”，在生物过程的开始阶段，细胞从根部开始，
# 沿着主干前进，直到到达第一个分支（如果有的话）。
# 然后，细胞必须选择一条路径，沿着树走得越来越远，直到到达一片叶子。
# 一个细胞的伪时间值是它回到根的距离。

cds <- orderCells(cds)
save(cds,file = "cdsAM1")
save(cds,file = "cdsAM2")
save(cds,file = "cdsAM3")
# load("cds")
# ⚠️使用root_state参数可以设置拟时间轴的根，如下面的拟时间着色图中可以看出，
# 左边是根。根据state图可以看出，根是State1，
# 若要想把另一端设为根，可以按如下操作
# #cds <- orderCells(cds, root_state = 5) #把State5设成拟时间轴的起始点
# install.packages("igraph")
# library(igraph)
# source(file = "order_cells.R")
# if(T){
#   root_state = NULL
#   num_paths = NULL
#   reverse = NULL
#   root_cell <- select_root_cell(cds, root_state, reverse)
#   cds@auxOrderingData <- new.env(hash = TRUE)
# 
#   if (cds@dim_reduce_type == "DDRTree") {
#     if (is.null(num_paths) == FALSE) {
#       message("Warning: num_paths only valid for method 'ICA' in reduceDimension()")
#     }
#     cc_ordering <- extract_ddrtree_ordering(cds, root_cell)
#     pData(cds)$Pseudotime <- cc_ordering[row.names(pData(cds)), ]$pseudo_time
#     K_old <- reducedDimK(cds)
#     old_dp <- cellPairwiseDistances(cds)
#     old_mst <- minSpanningTree(cds)
#     old_A <- reducedDimA(cds)
#     old_W <- reducedDimW(cds)
#     cds <- project2MST(cds, project_point_to_line_segment)
#     minSpanningTree(cds) <- cds@auxOrderingData[[cds@dim_reduce_type]]$pr_graph_cell_proj_tree
#     root_cell_idx <- which(V(old_mst)$name == root_cell, arr.ind = T)
#     cells_mapped_to_graph_root <- which(cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex == root_cell_idx)
#     if (length(cells_mapped_to_graph_root) == 0) {
#       cells_mapped_to_graph_root <- root_cell_idx
#     }
#     cells_mapped_to_graph_root <- V(minSpanningTree(cds))[cells_mapped_to_graph_root]$name
#     tip_leaves <- names(which(degree(minSpanningTree(cds)) == 1))
#     root_cell <- cells_mapped_to_graph_root[cells_mapped_to_graph_root %in% tip_leaves][1]
#     if (is.na(root_cell)) {
#       root_cell <- select_root_cell(cds, root_state, reverse)
#     }
#     cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell <- root_cell
#     cc_ordering_new_pseudotime <- extract_ddrtree_ordering(cds, root_cell)
#     pData(cds)$Pseudotime <- cc_ordering_new_pseudotime[row.names(pData(cds)), ]$pseudo_time
#     if (is.null(root_state) == TRUE) {
#       closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
#       pData(cds)$State <- cc_ordering[closest_vertex[, 1], ]$cell_state
#     }
#   }
# }

# 可视化
## 以pseudotime值上色 
# (Pseudotime是monocle2基于细胞基因表达信息计算的概率，表示时间的先后。)
pdf("./AM/AM1train.monocle.pseudotime.pdf",width = 7,height = 7)
plot_cell_trajectory(cds,color_by="Pseudotime", size=1,show_backbone=TRUE) 
dev.off()
pdf("./AM/AM2train.monocle.pseudotime.pdf",width = 7,height = 7)
plot_cell_trajectory(cds,color_by="Pseudotime", size=1,show_backbone=TRUE) 
dev.off()
pdf("./AM/AM3train.monocle.pseudotime.pdf",width = 7,height = 7)
plot_cell_trajectory(cds,color_by="Pseudotime", size=1,show_backbone=TRUE) 
dev.off()

##  以细胞类型上色
pdf("./AM/AM1train.monocle.celltype.pdf",width = 7,height = 7)
plot_cell_trajectory(cds,color_by="cell_type", size=1,show_backbone=TRUE)
dev.off()
pdf("./AM/AM2train.monocle.celltype.pdf",width = 7,height = 7)
plot_cell_trajectory(cds,color_by="cell_type", size=1,show_backbone=TRUE)
dev.off()
pdf("./AM/AM3train.monocle.celltype.pdf",width = 7,height = 7)
plot_cell_trajectory(cds,color_by="cell_type", size=1,show_backbone=TRUE)
dev.off()

## 以细胞状态上色
pdf("./AM/AM1train.monocle.state.pdf",width = 7,height = 7)
plot_cell_trajectory(cds, color_by = "State",size=1,show_backbone=TRUE)
dev.off()
pdf("./AM/AM2train.monocle.state.pdf",width = 7,height = 7)
plot_cell_trajectory(cds, color_by = "State",size=1,show_backbone=TRUE)
dev.off()
pdf("./AM/AM3train.monocle.state.pdf",width = 7,height = 7)
plot_cell_trajectory(cds, color_by = "State",size=1,show_backbone=TRUE)
dev.off()
# state的多少是monocle算出来的，不能调整，
# 与输入的用于轨迹学习的基因有关。
# 分叉和顶点之间或者顶点和顶点之间为一个state，
# 与发育轨迹时间先后没有关系，与细胞类型也不完全相关

## 以细胞状态上色（拆分）“分面”轨迹图，以便更容易地查看每个状态的位置。
pdf("./AM/AM1train.monocle.state.faceted.pdf",width = 10,height = 7)
plot_cell_trajectory(cds, color_by = "State") + facet_wrap("~State", nrow = 1)
dev.off()
pdf("./AM/AM2train.monocle.state.faceted.pdf",width = 10,height = 7)
plot_cell_trajectory(cds, color_by = "State") + facet_wrap("~State", nrow = 1)
dev.off()
pdf("./AM/AM3train.monocle.state.faceted.pdf",width = 10,height = 7)
plot_cell_trajectory(cds, color_by = "State") + facet_wrap("~State", nrow = 1)
dev.off()

## scale_color_manual()自己设置颜色
library(ggsci)
p1=plot_cell_trajectory(cds, color_by = "cell_type")  + scale_color_npg() 
p2=plot_cell_trajectory(cds, color_by = "State")  + scale_color_nejm()
colour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080")
p3=plot_cell_trajectory(cds, color_by = "State")  + scale_color_manual(values = colour)
p1|p2|p3

## 觉得这种轨迹不太直观的也可以画成树形图
library(patchwork)
p1 <- plot_cell_trajectory(cds, x = 1, y = 2, color_by = "cell_type") + 
  theme(legend.position='none',panel.border = element_blank())  #去掉第一个的legend
  # scale_color_manual(values = colour) 
p2 <- plot_complex_cell_trajectory(cds, x = 1, y = 2,
                                   color_by = "cell_type")+
  # scale_color_manual(values = colour) +
  theme(legend.title = element_blank()) 
p1|p2

# 还可以画沿时间轴的细胞密度图
library(ggpubr)
load("cdsAM2")
df <- pData(cds) 
## pData(cds)取出的是cds对象中cds@phenoData@data的内容
View(df)
p <- ggplot(df, aes(Pseudotime, colour = cell_type, fill=cell_type)) +
  geom_density(bw=0.5,size=1,alpha = 0.5)+theme_classic2()
ggsave("./AM/AM1density_pseudotimeplot.pdf", plot = p, width = 16, height = 8)
p <- ggplot(df, aes(Pseudotime, colour = cell_type, fill=cell_type)) +
  geom_density(bw=0.5,size=1,alpha = 0.5)+theme_classic2()
ggsave("./AM/AM2density_pseudotimeplot.pdf", plot = p, width = 16, height = 8)
p <- ggplot(df, aes(Pseudotime, colour = cell_type, fill=cell_type)) +
  geom_density(bw=0.5,size=1,alpha = 0.5)+theme_classic2()
ggsave("./AM/AM3density_pseudotimeplot.pdf", plot = p, width = 16, height = 8)

# 手动设置颜色注意
ClusterName_color_panel <- c(
  "Naive CD4 T" = "#DC143C", "Memory CD4 T" = "#0000FF", "CD14+ Mono" = "#20B2AA",
  "B" = "#FFA500", "CD8 T" = "#9370DB", "FCGR3A+ Mono" = "#98FB98",
  "NK" = "#F08080", "DC" = "#0000FF", "Platelet" = "#20B2AA"
)
ggplot(df, aes(Pseudotime, colour = cell_type, fill=cell_type)) +
  geom_density(bw=0.5,size=1,alpha = 0.5)+
    theme_classic2()+ 
    scale_fill_manual(name = "", values = ClusterName_color_panel)+
    scale_color_manual(name = "", values = ClusterName_color_panel)

# 提取感兴趣的细胞（进行后续分析）
#比如对State7的细胞感兴趣
pdata <- Biobase::pData(cds)
s.cells <- subset(pdata, State=="7") %>% rownames()
save(s.cells, file = "Monocle_state7.rda")

# save cds
write.csv(pData(cds), "pseudotime.csv")
save(cds, file = "cds.rda")

# 指定基因的可视化
##选择前4个top基因并将其对象取出
keygenes <- head(ordergene,4)
cds_subset <- cds[keygenes,]
##可视化：以state/celltype/pseudotime进行
p1 <- plot_genes_in_pseudotime(cds_subset, color_by = "State")
p2 <- plot_genes_in_pseudotime(cds_subset, color_by = "cell_type")
p3 <- plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime")
plotc <- p1|p2|p3
ggsave("./AM/AM1Genes_pseudotimeplot.pdf", plot = plotc, width = 16, height = 8)


#指定基因
s.genes <- c("SELL","CCR7","IL7R", "CD84","CCL5","S100A4")
p1 <- plot_genes_jitter(cds[s.genes,], grouping = "State", color_by = "State")
p2 <- plot_genes_violin(cds[s.genes,], grouping = "State", color_by = "State")
p3 <- plot_genes_in_pseudotime(cds[s.genes,], color_by = "State")
plotc <- p1|p2|p3
ggsave("Genes_Jitterplot.pdf", plot = plotc, width = 16, height = 8)

# 由于state1-7并不代表分化时间先后，因此这张图还要结合其他的图来看
# 拟时序展示单个基因表达量
colnames(pData(cds))
pData(cds)$CCL5 = log2( exprs(cds)['CCL5',]+1)
p1=plot_cell_trajectory(cds, color_by = "CCL5")  + scale_color_gsea()
pData(cds)$S100A4 = log2(exprs(cds)['S100A4',]+1)
p2=plot_cell_trajectory(cds, color_by = "S100A4")    + scale_color_gsea()
library(patchwork)
p1+p2

# 寻找拟时相关的基因（拟时差异基因）
# 使用回归算法
# 注意：不要使用多核运算，经常会出现警告
# Monocle的主要工作是通过生物过程（如细胞分化）将细胞按顺序排列，
# 而不知道要提前查看哪些基因。
# 一旦这样做了，你就可以分析细胞，找到随着细胞进展而变化的基因
# 官方给出的差异分析有三大方法：
# 1、Basic Differential Analysis
# 2、Finding Genes that Distinguish Cell Type or State
# 3、Finding Genes that Change as a Function of Pseudotime
# 我们重点关注第三个：根据伪时间功能寻找差异基因
# sm.ns函数指出Monocle应该通过表达式值拟合自然样条曲线，
# 以帮助它将表达式的变化描述为进程的函数。
# 寻找拟时差异基因（qvalue体现基因与拟时的密切程度）绘制热图
#这里是把排序基因（ordergene）提取出来做回归分析，来找它们是否跟拟时间有显著的关系
#如果不设置，就会用所有基因来做它们与拟时间的相关性
Time_diff <- differentialGeneTest(cds[ordergene,], cores = 1, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
Time_diff <- Time_diff[,c(5,2,3,4,1,6,7)] #把gene放前面，也可以不改
write.csv(Time_diff, "Time_diff_all.csv", row.names = F)
Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()
p=plot_pseudotime_heatmap(cds[Time_genes,], num_clusters=4, show_rownames=T, return_heatmap=T)
ggsave("Time_heatmapAll.pdf", p, width = 5, height = 10)

