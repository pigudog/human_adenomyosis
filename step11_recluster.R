reticulate::py_install(packages ='umap-learn')
reticulate::py_install(packages ='scikit-learn')
py_module_available(module = 'umap')
py_module_available(module = 'sklearn')
Sys.which("python")
library(reticulate)
library(patchwork)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)

options(stringsAsFactors = FALSE)
rm(list=ls())
load(file = 'AM_after_cluster.Rdata')
load(file="AM2_HVG.Rdata")
AM@meta.data$cell_type <- AM@active.ident

ps <- function(filename, plot = FALSE, w = 9, h = 6) {
  if (is.object(plot)) {
    print(plot)
  }
  plot <- recordPlot()
  pdf(file = paste0(filename,".pdf"), onefile = T, width = w, height = h)
  replayPlot(plot)
  dev.off()
}
n_jobs <- 48 # 使用的线程数

# AM2, 第二个monocle：fibro-like，unkonwn2，SMA
AM2 <- AM[,AM$cell_type %in% c("VSMC", "Pericyte","SMC","ESC","Unkown2")]
AM2 <- SCTransform(AM2, vst.flavor = "v2")
AM2 <- RunPCA(AM2, npcs=50, verbose=FALSE)
ElbowPlot(AM2)
pc.num=1:20
AM2 <- RunUMAP(AM2, dims=pc.num)
# 0.5
AM2 <- FindNeighbors(AM2, dims = pc.num) 
AM2<-FindClusters(AM2,resolution=0.8)
DimPlot(AM2)
p <- DimPlot(AM2, group.by = "orig.ident")
ggsave("AM2/UMAP_recluster_before_singler.pdf", p, width = 8, height = 6)
p <- DimPlot(AM2, group.by = "orig.ident", split.by = "orig.ident", ncol = 4)
ggsave("AM2/UMAP_Samples_Split_before_singler.pdf", p, width = 18, height = 12)
p <- DimPlot(AM2,  label = T, group.by = "seurat_clusters")

AM2.markers <- FindAllMarkers(object = AM2,assay = "SCT", only.pos = TRUE, min.pct = 0.25,thresh.use = 0.25)
# wriet to placental_decidual/HVGmarkers_res_0.5
write.csv(AM2.markers,file=paste0("AM2",'/HVGmarkers_res_0.8.csv'))

DimPlot(AM2,label = T)
#  Need to look at the picture, determine the cell subsets:
celltype=data.frame(ClusterID=0:15,
                    celltype='unkown')
celltype[celltype$ClusterID %in% c(3,5,11,14),2]='SMC'

# celltype[celltype$ClusterID %in% c(2,11,12),2]='Fibro-like'
# celltype[celltype$ClusterID %in% c(0,1,8,10,15),2]='Actin-Fibro'
# celltype[celltype$ClusterID %in% c(3),2]='Infra-Fibro'
celltype[celltype$ClusterID %in% c(4,9,10),2]='ESC'
celltype[celltype$ClusterID %in% c(0),2]='Unkown2'
celltype[celltype$ClusterID %in% c(2,6,7,8,13),2]='Pericyte'
celltype[celltype$ClusterID %in% c(1,12,15),2]='VSMC'
table(celltype)
new.cluster.ids <- celltype$celltype
names(new.cluster.ids) <- levels(AM2)
AM2 <- RenameIdents(AM2, new.cluster.ids)

library(dplyr) 
top10 <- AM2.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
p <- DoHeatmap(AM2,top10$gene,size=3)
ggsave(p,filename=paste0("AM2",'/top10_heatmap.pdf'),width = 24,height = 18)
library(cowplot)
library(tidydr)
library(stringr)
library(viridis)
library(scCustomize)
genes_to_check = c(
  # "ITGB1", "MME","ANPEP",#MSC
                   "PDGFRB","GJA4","RGS5", #Pericyte
                   "DES", "ACTA2",# SMC
                    
                   "PDGFRA","COL1A2" # ESC
                   
                   # "COL1A1", "COL3A1","MDK" # Fibro

)
pal <- viridis_magma_dark_high
plots <- list()
for (i in 1:length(genes_to_check)){
  plots[[i]]=FeaturePlot_scCustom(seurat_object = AM2,
                                  colors_use = pal,
                                  features = genes_to_check[i])+NoAxes()
}
p<-wrap_plots(plots, ncol = 3)
ggsave(p,filename=paste0("AM2",'/second _fearureplot_sepcify.pdf'),width = 16,height = 12)
pal2 <- pal[50:100]
p <- DotPlot_scCustom(AM2, features = genes_to_check,colors_use = pal2,flip_axes = T,x_lab_rotate = TRUE,
                      remove_axis_titles = FALSE) + coord_flip()
ggsave(p,filename=paste0("AM2",'/second_dotplot_sepcify.pdf'),width = 6,height = 5)
library(ggsci)
library(Cairo)
AM2@meta.data$cell_type <- AM2@active.ident
color = c(pal_d3("category20")(20),
          pal_d3("category20b")(20),
          pal_d3("category20c")(20),
          pal_d3("category10")(10))
#设置group顺序
# AM@meta.data$cell_type <- factor(AM@active.ident )
df=AM2@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  cbind(cell_type=AM2@meta.data$cell_type)

CairoPDF("./AM2/umap_ellipse.pdf", width=12, height=9) 
#用stat_ellipse加置信区间
p_UMAP1=ggplot(df, aes(UMAP_1, UMAP_2, color=cell_type))+
  geom_point(size = 1) + #细胞亚群点的大小
  scale_color_manual(values=color)+
  theme(panel.border = element_blank(), #边框
        axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),
        # panel.grid = element_blank(),
        # panel.grid.major = element_blank(), #网格线
        # panel.grid.minor = element_blank(), #网格线
        panel.background = element_rect(fill = 'white'), #背景色
        plot.background=element_rect(fill="white"),
        legend.title = element_blank(), #去掉legend.title 
        legend.key=element_rect(fill='white'), 
        legend.text = element_text(size=14), #设置legend标签的大小
        legend.key.size=unit(0.6,'cm') )+
  stat_ellipse(aes(x = UMAP_1, y = UMAP_2, fill = cell_type),
               geom = "polygon",#多边形
               linetype=2, # 椭圆，虚线
               alpha = 0.05,  #椭圆不透明度
               show.legend = FALSE, #去掉椭圆对应的图例
               level = 0.93)+ #置信区间
  guides(fill= guide_legend(override.aes = list(size = 4)))+ #图例点的大小
  scale_color_manual(values=color); p_UMAP1 
dev.off()
#theme_dr加左下角的坐标
#呈现效果如图，如果还想更加精确的加上细胞群外圈的框框的话，PPT和AI或许是更好的选择。
#theme_dr
CairoPDF("./AM2/UMAP_theme_dr.pdf", width=12, height=9) 
p_UMAP2<-p_UMAP1+theme_dr(xlength = 0.22, ylength = 0.22,
                          arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")
)+theme(panel.grid = element_blank());p_UMAP2
dev.off()

#先提取出tSNE_1的位置信息
celltype_position <- df %>%
  group_by(cell_type) %>%
  dplyr::summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2))

library(ggrepel)
CairoPDF("./AM2/UMAP_label.pdf", width=12, height=9) 
p_UMAP3<-p_UMAP2+#geom_point(aes(color=factor(celltype)), size=3)+
  geom_label_repel(data = celltype_position,aes(label=cell_type), 
                   fontface="bold",point.padding=unit(0.1, "lines"))+theme(legend.position = "none");p_UMAP3

dev.off()
#分样本
sample_table <- as.data.frame(table(AM2@meta.data$orig.ident,AM2@meta.data$cell_type))
names(sample_table) <- c("Samples","celltype","CellNumber")

plot_sample<-ggplot(sample_table,aes(x=Samples,weight=CellNumber,fill=celltype))+
  geom_bar(position="fill")+
  scale_fill_manual(values=color) + 
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        axis.line.x = element_line(colour = "black") ,
        axis.line.y = element_line(colour = "black") ,
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16)
  )+labs(y="Percentage")+RotatedAxis()
plot_sample
p1 <- DimPlot(AM2, reduction = 'umap', 
              group.by = "seurat_clusters",
              cols = color , #设置颜色
              repel = T,  #label不重叠
              label.size = 5, #调整label的大小
              label = TRUE,  #是否展示label
              pt.size = 0.5);p1
p2 <- DimPlot(AM2, reduction = 'umap', 
              group.by = "cell_type",
              cols = color , 
              label = F, pt.size = 0.5);p2
p3 <- DimPlot(AM2, reduction = 'umap', group.by = "orig.ident",cols = color , 
              label = F, pt.size = 0.5);p3
P_final <- (p_UMAP3 +plot_sample)/(p1+ p3)
ggsave(P_final,filename="./AM2/P_final.png",width=10, height=10)



save(AM2,AM2.markers,file = 'AM2_HVG.Rdata')
