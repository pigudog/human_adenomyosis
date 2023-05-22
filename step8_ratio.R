library(Seurat) 
library(tidyverse)
library(patchwork) #拼图
library(scCustomize)
rm(list=ls())
# devtools::install_github("sajuukLyu/ggunchull", type = "source")
# install.packages("tidydr")
library(ggplot2)
library(cowplot)
library(dplyr)
library(ggunchull)
library(tidydr)
library(ggsci)
library(Cairo)
#读取数据
load(file = 'AM_after_cluster.Rdata')


color = c(pal_d3("category20")(20),
          pal_d3("category20b")(20),
          pal_d3("category20c")(20),
          pal_d3("category10")(10))
color2=c("#D9AB42","#A35E47","#0F4C3A","#563F2E","#78C2C4","#C73E3A","#B47157","#2B5F75")
AM@meta.data$cell_type <- AM@active.ident
DimPlot(AM, reduction = "umap", group.by = "cell_type",
        cols = color,
        pt.size = 1.5,
        label = T,label.box = T # 拥有方框的标签
) 

df=AM@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  cbind(cell_type=AM@meta.data$cell_type)

CairoPDF("./AM/umap_ellipse.pdf", width=12, height=9) 
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
CairoPDF("./AM/UMAP_theme_dr.pdf", width=12, height=9) 
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
CairoPDF("./AM/UMAP_label.pdf", width=12, height=9) 
p_UMAP3<-p_UMAP2+#geom_point(aes(color=factor(celltype)), size=3)+
  geom_label_repel(data = celltype_position,aes(label=cell_type), 
                   fontface="bold",point.padding=unit(0.1, "lines"))+theme(legend.position = "none");p_UMAP3

dev.off()







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
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16)
  )+labs(y="Percentage")+RotatedAxis()
plot_sample

# #分组
# group_table <- as.data.frame(table(AM@meta.data$cell_type,AM@active.ident))
# names(group_table) <- c("group","celltype","CellNumber")
# plot_group<-ggplot(AM@meta.data,aes(x=group,fill=celltype))+
#   geom_bar(position="fill")+
#   scale_fill_manual(values=color) + 
#   theme(panel.grid = element_blank(),
#         panel.background = element_rect(fill = "transparent",colour = NA),
#         axis.line.x = element_line(colour = "black") ,
#         axis.line.y = element_line(colour = "black") ,
#         plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16)
#   )+labs(y="Percentage")
# plot_group

library(paletteer)
pal <- paletteer_d("ggsci::nrc_npg")[c(1,3,4,9,5)]
plot1 <- DimPlot(AM, label = T, pt.size = 1,cols = pal)+
  NoLegend()+labs(x = "UMAP1", y = "UMAP2",title = "Celltype") +
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

library(RColorBrewer)
colour <- c(brewer.pal(9, "Set1"), 
            "#FF34B3","#BC8F8F","#20B2AA","#00F5FF","#FFA500","#ADFF2F",
            "#FF6A6A","#7FFFD4", "#AB82FF","#90EE90","#00CD00","#008B8B",
            "#6495ED","#FFC1C1","#CD5C5C","#8B008B","#FF3030", "#7CFC00")  
            
plot7 <- DimPlot(AM, label = T, pt.size = 1,cols = cell_type_cols)+
  NoLegend()+labs(x = "UMAP1", y = "UMAP2",title = "Celltype") +
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
plot_grid(plot6,plot7)


#设置group顺序
AM@meta.data$cell_type <- factor(AM@active.ident )

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
#简单的拼图
p_UMAP3 /(p1+p2 + p3)


plot_sample + 
  # plot_group + 
  plot_layout(widths = c(1, 1))

(p_UMAP3 +plot_sample)/(p1+p2 + p3)+ 
  plot_annotation(title = "FIG1", tag_levels = "A")



library(clustree)
library(cowplot)
genes_to_check = c("PDGFRB","GJA4","RGS5", #Pericyte
                   "PDGFRA","COL1A2", # ESC
                   "PECAM1", "CDH5",  # Endo
                   "DES", "ACTA2",# SMC
                   "KRT8","KRT18", # Epi 
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
pal2 <- pal[50:100]
p_dotplot <- DotPlot_scCustom(AM, features = genes_to_check,colors_use = pal2,flip_axes = T,x_lab_rotate = TRUE,
                      remove_axis_titles = FALSE) + coord_flip()


P_final <- (p_UMAP3 +plot_sample)/(p1+ p3)
ggsave(P_final,filename="./AM/P_final.png",width=12, height=12)
