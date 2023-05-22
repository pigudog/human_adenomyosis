rm(list=ls())
library(Seurat) 
library(tidyverse)
library(patchwork) #拼图
library(scCustomize)
library(ggplot2)
library(cowplot)
library(dplyr)
library(ggunchull)
library(tidydr)
library(ggsci)
library(Cairo)
load(file = 'AM_after_cluster.Rdata')
ps <- function(filename, plot = FALSE, w = 9, h = 6) {
  if (is.object(plot)) {
    print(plot)
  }
  plot <- recordPlot()
  pdf(file = paste0(filename,".pdf"), onefile = T, width = w, height = h)
  replayPlot(plot)
  dev.off()
}

# color
color = c(pal_d3("category20")(20),
          pal_d3("category20b")(20),
          pal_d3("category20c")(20),
          pal_d3("category10")(10))

## picture B
AM@meta.data$cell_type <- AM@active.ident
df=AM@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  cbind(cell_type=AM@meta.data$cell_type)
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
  scale_color_manual(values=color) 
#theme_dr加左下角的坐标
#呈现效果如图，如果还想更加精确的加上细胞群外圈的框框的话，PPT和AI或许是更好的选择。
#theme_dr
p_UMAP2<-p_UMAP1+theme_dr(xlength = 0.22, ylength = 0.22,
                          arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")
)+theme(panel.grid = element_blank())


#先提取出tSNE_1的位置信息
celltype_position <- df %>%
  group_by(cell_type) %>%
  dplyr::summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2))

library(ggrepel)
CairoPDF("./Figure1/picture_B.pdf", width=8, height=6) 
p_UMAP3<-p_UMAP2+geom_point(aes(color=factor(cell_type)), size=3)+
  geom_label_repel(data = celltype_position,aes(label=cell_type), 
                   fontface="bold",point.padding=unit(0.2, "lines"))+theme(legend.position = "none");p_UMAP3

dev.off()


## pictureC
#分样本
sample_table <- as.data.frame(table(AM@meta.data$orig.ident,AM@active.ident))
names(sample_table) <- c("Samples","celltype","CellNumber")

plot_sample<-ggplot(sample_table,aes(x=Samples,weight=CellNumber,fill=celltype))+
  geom_bar(position="fill")+
  scale_fill_manual(values=color) + 
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        axis.line.x = element_line(colour = "black") ,
        axis.line.y = element_line(colour = "black") ,
        axis.title = element_text(lineheight=2, face="bold", hjust=2, size =15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        plot.title = element_text(size = 20)
  )+labs(y="Percentage")+RotatedAxis()
plot_sample
ps("./Figure1/picture_C",w = 6, h =8)
dev.off()

## pictureD
p1 <- DimPlot(AM, reduction = 'umap', 
              group.by = "seurat_clusters",
              cols = color , #设置颜色
              repel = T,  #label不重叠
              label.size = 5, #调整label的大小
              label = TRUE,  #是否展示label
              pt.size = 0.5);p1
p2 <- DimPlot(AM, reduction = 'umap', 
              group.by = "cell_type",
              cols = color , 
              label = F, pt.size = 0.5);p2
p3 <- DimPlot(AM, reduction = 'umap', group.by = "orig.ident",cols = color , 
              label = F, pt.size = 0.5);p3

P_d <- (p1+ p3)
ps("./Figure1/picture_D",w = 10, h =5)
dev.off()

# pictureE
library(cowplot)
genes_to_check = c("PDGFRB","GJA4","RGS5", #Pericyte
                   "PDGFRA","COL1A2", # ESC
                   "PECAM1", "CDH5",  # Endo
                   "KRT8","KRT18", # Epi 
                   "DES", "ACTA2",# SMC
                   
                   # "COL1A1", "COL3A1","MDK", # Fibro
                   "GZMB", # NK
                   "CD3G","CD3D",  # T cells
                   "TPSB2","TPSAB1", # Mast
                   "CD14","MS4A7" # Mc
)
pal <- viridis_magma_dark_high
plots <- list()
for (i in 1:length(genes_to_check)){
  plots[[i]]=FeaturePlot_scCustom(seurat_object = AM,
                                  colors_use = pal,
                                  features = genes_to_check[i])+NoAxes()
}
p_featureplot<-wrap_plots(plots, ncol = 3)
ggsave(p_featureplot,filename="./Figure1/picture_F.png",width = 16,height = 12)

pal2 <- pal[0:100]
p_dotplot <- DotPlot_scCustom(AM, features = genes_to_check,colors_use = pal2,flip_axes = T,x_lab_rotate = TRUE,
                              remove_axis_titles = FALSE) + coord_flip()+ 
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        axis.line.x = element_line(colour = "black") ,
        axis.line.y = element_line(colour = "black") ,
        axis.title = element_text(lineheight=2, face="bold", hjust=2, size =15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        plot.title = element_text(size = 20)
  );p_dotplot
ps("./Figure1/picture_E",w = 6, h =6)

# 设置可能用到的主题
library(ggpubr)
library(ggsignif)
my_comparisons <- list(c("CTRL","EC"),c("EC","EM"),c("CTRL","EM"))
theme.set2 = theme(axis.title.x=element_blank())
# 设置绘图元素
# PDGFRA和PDGFC
plot.featrures = c("MMP1","CSF3","CXCL8","IGF1")
group = "orig.ident"
# 质控前小提琴图
plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(AM, group.by=group, pt.size = 0,
                       features = plot.featrures[i])+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    stat_compare_means(comparisons = my_comparisons) +
    ylim(-1, 8)+
    geom_boxplot(position=position_dodge(0.6),
                 size=0.4,
                 width=0.1,
                 # fill="gold",
                 # color="blue",
                 # outlier.color = "blue",
                 # outlier.fill = "red",
                 outlier.shape = 19,
                 outlier.size = 1.5,
                 outlier.stroke = 0.5,
                 outlier.alpha = 45,
                 notch = F,
                 notchwidth = 0.5)+
    theme.set2 + NoLegend()}
p <- wrap_plots(plots = plots, nrow=1)     
ggsave(p,filename=paste0("Figure1",'/umap_MMP1.pdf'),width = 16,height = 5)
