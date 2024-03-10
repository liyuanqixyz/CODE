##Fig1
load("../data/integrated.final.Rdata")
Idents(integrated) <- integrated$Global_Cluster

library(ggplot2)
library(ggalluvial)
library(RColorBrewer)
library(ggpubr)
mycolor <- brewer.pal(12, 'Set3')[c(1,3:8,10,11)]
xx <- c("#8DD3C7","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#BC80BD","#CCEBC5")
sce <- integrated
DimPlot(sce,label = T,cols =mycolor,pt.size = 0.5,label.size = 8)+theme(axis.line = element_line(size = 1),
                                                                        legend.text = element_text(size = 16))

DimPlot(sce,label = F,cols =mycolor,split.by = "Tissue",label.size = 16)+theme(axis.line = element_line(size = 1),
                                                                               legend.text = element_text(size = 16))
DimPlot(sce,label = F,cols =c("#6C9D88","#6980AC"),group.by = "Tissue")+theme(axis.line = element_line(size = 0.5),
                                                                              legend.text = element_text(size = 12))

DimPlot(sce,label = F,cols =mycolor,group.by = "Group",label.size = 16,pt.size = 0.005)

table(sce$Group)
prop.table(table(Idents(sce)))
table(Idents(sce), sce$Tissue)
Cellratio <- prop.table(table(Idents(sce), sce$Group), margin = 2)
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio)[1]<- "celltype"
colnames(Cellratio)[2]<- "Samples"
Cellratio$group <- ifelse(substr(Cellratio$Samples,4,4)=="N","Normal","Tumor")

library(ggplot2)
Cellratio$Samples <- factor(Cellratio$Samples,levels = c("P01T","P02T","P03T","P05T","P01N","P02N","P03N","P04N"))

p1 <- ggplot(Cellratio) + 
  geom_bar(aes(x =Samples, y= Freq, fill = celltype),stat = "identity",width = 0.9,size = 0.2)+
  theme_classic() +
  labs(x=NULL,y = 'Ratio of cell types')+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size = 15),
                                               axis.line = element_line(size = 0.8),
                                               legend.text = element_text(size = 12),
                                               axis.text.y = element_text(size = 12),
                                               axis.title.y = element_text(size = 15))+
  scale_fill_manual(values = mycolor)

ggsave(p1,width = 8.02,height = 6.68,file="../PDF/Fig1/barplot1.pdf")


p<-ggplot(Cellratio,aes(x=Samples,y=Freq,fill=celltype,stratum=celltype,alluvium=celltype))+
  geom_col(width = 0.4,color=NA)+
  geom_flow(width=0.4,alpha=0.5,knot.pos=0)+theme_bw()+theme(axis.title.x = element_blank(),axis.text.x = element_text(size = 10))+
  ylab(label = "Percentages")+scale_fill_manual(values = mycolor)

p+mytheme


ggplot(Cellratio, aes(x = celltype, y = Freq))+
  geom_boxplot(aes(fill=group))+ 
  stat_compare_means(aes(group = group),
                     label = "p.signif", 
                     method = "wilcox.test",
                     hide.ns = T,
                     size = 10)+theme_bw()+
  labs(x='Celltype',y = 'Percentages')+theme(axis.title.x = element_blank(),
                                             axis.text.x = element_text(hjust = 1,vjust = 1,angle = 45,size = 17),
                                             axis.line = element_line(size = 0.8),
                                             legend.text = element_text(size = 17),
                                             axis.text.y = element_text(size = 17),
                                             axis.title.y = element_text(size = 17))+
  scale_fill_manual(values = c("#6C9D88","#6980AC"))

###Dotplot###
genes_to_check = c('KLRB1','NKG7',#NK
                   'PECAM1', 'RAMP2', ## endo
                   'MKI67','TOP2A',#cycling
                   'CD68','C1QB', ##macro
                   'CD3D', 'CD3E',#T
                   "RGS5","COL1A1",#hsc
                   "ALB","SERPINA1",#hcc
                   'CD79A', 'MS4A1' ,  #B
                   'IGHG1', 'MZB1') #plasma cells

genes_to_check=unique(genes_to_check) 
p <- DotPlot(sce, features = unique(genes_to_check),
             assay='RNA' ,group.by = 'Global_Cluster' ,cols = "RdYlBu")+
  labs(y = "Cell type", x = "Features")+ coord_flip()+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size = 15),
                                                            axis.title.x = element_blank(),
                                                            axis.title.y = element_blank(),
                                                            axis.text.y = element_text(size = 15))

#######
load("../data/integrated.final.Rdata")

Idents(integrated) <- integrated$Tissue
integrated <- subset(integrated,idents = "Tumor")

Idents(integrated) <- integrated$Sub_Cluster
integrated <- subset(integrated,idents = c("Endothelial cells","Cycling-MKI67","Fibroblasts","Cyclings","Hepatocytes"),invert=T)
markers <- read.csv("../result/integrated.final.csv",row.names = 1)
markers <- markers[markers$p_val_adj<0.05,]
markers <- markers[markers$cluster%in%c("NK cells","DC-CD1C","DC-CLEC9A","Macrophages-C1QC",
                                        "Macrophages-SPP1","Monocytes-S100A8","Mast cells-TPSAB1",
                                        "CD8-GNLY","CD8-GZMK","CD4-ANXA1","CD8-XCL1","CD4-IL7R","CD8-SLC4A10",
                                        "CD4-FOXP3","B cells","Plasma cells"),]

meta <- integrated@meta.data
meta$Sub_Cluster <- as.character(meta$Sub_Cluster)
genes <- unique(markers$gene)

a <- integrated@assays$RNA@data[genes,]
dat <- as.data.frame(a)
#dat<- as.data.frame(integrated@assays$RNA@data)
colnames(dat) <- as.character(meta$Sub_Cluster)
dat2 <- rownames_to_column(dat,var = "gene")
write.table(dat2,file = "../result/scRNA-Seq_reference_tumor_imm.txt",sep = "\t",quote =F,row.names = F)


########
sce<- subset(integrated,idents = c("Cyclings"),invert = T)
imm <- subset(sce,idents = c("NK cells","Myeloids","T cells","B cells","Plasma cells"))
nonimm <- subset(sce,idents = c("Endothelial cells","Fibroblasts","Hepatocytes"))

imm_meta <- imm@meta.data
imm_meta$immunity <- "Immune cells"
nonimm_meta <- nonimm@meta.data
nonimm_meta$immunity <- "Non-immune cells"
meta <- rbind(imm_meta,nonimm_meta)
sce <- AddMetaData(sce,metadata = meta)

table(sce$immunity)
Idents(sce) <- sce$immunity
Cellratio <- prop.table(table(Idents(sce), sce$Group), margin = 2)
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio)[1]<- "celltype"
colnames(Cellratio)[2]<- "Samples"
Cellratio$group <- ifelse(substr(Cellratio$Samples,4,4)=="N","Normal","Tumor")
Cellratio$Samples <- factor(Cellratio$Samples,levels = c("P01T","P02T","P03T","P05T","P01N","P02N","P03N","P04N"))


p1 <- ggplot(Cellratio) + 
  geom_bar(aes(x =Samples, y= Freq, fill = celltype),stat = "identity",width = 0.9,size = 0.2)+
  theme_classic() +
  labs(x=NULL,y = 'Ratio of cell types')+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size = 15),
                                               axis.line = element_line(size = 0.8),
                                               legend.text = element_text(size = 12),
                                               axis.text.y = element_text(size = 12),
                                               axis.title.y = element_text(size = 15))+
  scale_fill_manual(values = mycolor)
ggsave(p1,width = 8.02,height = 6.68,file="../PDF/Fig1//imm_nonimm.pdf")
########
Cellratio <- prop.table(table(Idents(imm), imm$Group), margin = 2)
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio)[1]<- "celltype"
colnames(Cellratio)[2]<- "Samples"
Cellratio$group <- ifelse(substr(Cellratio$Samples,4,4)=="N","Normal","Tumor")
Cellratio$Samples <- factor(Cellratio$Samples,levels = c("P01T","P02T","P03T","P05T","P01N","P02N","P03N","P04N"))

p2 <- ggplot(Cellratio) + 
  geom_bar(aes(x =Samples, y= Freq, fill = celltype),stat = "identity",width = 0.9,size = 0.2)+
  theme_classic() +
  labs(x=NULL,y = 'Ratio of cell types')+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size = 15),
                                               axis.line = element_line(size = 0.8),
                                               legend.text = element_text(size = 12),
                                               axis.text.y = element_text(size = 12),
                                               axis.title.y = element_text(size = 15))+
  scale_fill_manual(values = c("#8DD3C7","#80B1D3","#FDB462","#BC80BD","#CCEBC5"))
ggsave(p2,width = 8.02,height = 6.68,file="../PDF/Fig1/imm_sub.pdf")


Cellratio <- prop.table(table(Idents(nonimm), nonimm$Group), margin = 2)
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio)[1]<- "celltype"
colnames(Cellratio)[2]<- "Samples"
Cellratio$group <- ifelse(substr(Cellratio$Samples,4,4)=="N","Normal","Tumor")
Cellratio$Samples <- factor(Cellratio$Samples,levels = c("P01T","P02T","P03T","P05T","P01N","P02N","P03N","P04N"))

p3 <- ggplot(Cellratio) + 
  geom_bar(aes(x =Samples, y= Freq, fill = celltype),stat = "identity",width = 0.9,size = 0.2)+
  theme_classic() +
  labs(x=NULL,y = 'Ratio of cell types')+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size = 15),
                                               axis.line = element_line(size = 0.8),
                                               legend.text = element_text(size = 12),
                                               axis.text.y = element_text(size = 12),
                                               axis.title.y = element_text(size = 15))+
  scale_fill_manual(values = c("#BEBADA","#B3DE69","#FCCDE5"))
ggsave(p3,width = 8.02,height = 6.68,file="../PDF/Fig1/nonimm_sub.pdf")

#######Fig2######
Bcell<- subset(integrated,idents = c("B cells","Plasma cells"))
scRNA <- CreateSeuratObject(Bcell@assays$RNA@counts)

Sample <- substr(colnames(scRNA),1,4)
patient<- substr(Sample,1,3)
tissue <- ifelse(substr(Sample,4,4)=="N","Normal","Tumor")
group <- data.frame(Samples=Sample,Tissue=tissue,Patient=patient)
rownames(group) <- colnames(scRNA)
scRNA <- AddMetaData(scRNA,metadata =group)
table(scRNA@meta.data$Samples)
table(scRNA@meta.data$Tissue)

sce <- scRNA
sce <- NormalizeData(sce, 
                     normalization.method = "LogNormalize",
                     scale.factor = 1e4) 
sce <- FindVariableFeatures(sce,nfeatures = 3000)
sce <- ScaleData(sce)
sce <- RunPCA(sce, features = VariableFeatures(object = sce))
ElbowPlot(sce,ndims = 50)

##########subset B cells#######
library(harmony)
sce<- RunHarmony(sce, "Samples",max.iter.harmony = 20)
names(sce@reductions)
ElbowPlot(sce,ndims = 50,reduction = "harmony")

ndims <- 1:30
sce <- RunUMAP(sce,  dims = ndims, 
               reduction = "harmony")
sce <- FindNeighbors(sce, reduction = "harmony",
                     dims = ndims) 
sce <- FindClusters(sce,resolution = 0.3)
DimPlot(sce,label=T,pt.size = 0.5)

marker <- FindAllMarkers(Bcell,only.pos = T)
top10 = marker%>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)


DimPlot(Bcell,label = T,pt.size = 0.5)
DimPlot(Bcell,pt.size = 0.5,group.by = "Samples")
DimPlot(Bcell,pt.size = 0.5,group.by = "Tissue")

VlnPlot(Bcell,features = c("CD79A","MZB1","MS4A1","BCL6","IGHA1","IGHG1"))
FeaturePlot(Bcell,features = c("CD79A","MZB1","MS4A1","IGHG1"))
FeaturePlot(Bcell,features = c("IGKC","IGLC3","IGHG1","IGHA1"))

F1 <- FeaturePlot(Bcell,features = c("IL10"),cols = c("lightgrey","#E41A1C"))
ggsave(F1,filename = "../PDF/Fig2/featureplot_IL10.pdf",width = 4.46,height =3.28 )

F2 <- FeaturePlot(Bcell,features = c("IL12A"),cols = c("lightgrey","#E41A1C"))
ggsave(F2,filename = "../PDF/Fig2/featureplot_IL12A.pdf",width = 4.46,height =3.28 )

F3 <- FeaturePlot(Bcell,features = c("EBI3"),cols = c("lightgrey","#E41A1C"))
ggsave(F3,filename = "../PDF/Fig2/featureplot_EBI3.pdf",width = 4.46,height =3.28 )


Maincelltype=data.frame(ClusterID=0:4,
                        Maincelltype=0:4)   

Maincelltype[Maincelltype$ClusterID %in% c(0,1,4),2]='B cells'
Maincelltype[Maincelltype$ClusterID %in% c(2,3),2]='Plasma cells'

for(i in 1:nrow(Maincelltype)){
  Bcell@meta.data[which(Bcell@meta.data$RNA_snn_res.0.3 == Maincelltype$ClusterID[i]),'Maincelltype'] <- Maincelltype$Maincelltype[i]}

table(Bcell$Maincelltype)
Idents(Bcell) <- Bcell$Maincelltype
mycolor <- brewer.pal(12, 'Set3')[c(5,10)]

DimPlot(Bcell,label = T,cols = mycolor,pt.size = 0.5,label.size = 6)+theme(axis.line = element_line(size = 1),legend.text = element_text(size = 16))
DimPlot(Bcell,label = F,cols = c("#6C9D88","#6980AC"),pt.size = 0.5,label.size = 6,group.by = "Tissue")+theme(axis.line = element_line(size = 1),legend.text = element_text(size = 16))


prop.table(table(Idents(Bcell)))
table(Idents(Bcell), Bcell$Tissue)#

Cellratio <- prop.table(table(Idents(Bcell), Bcell$Tissue), margin = 2)#
Cellratio <- as.data.frame(Cellratio)

Tumor <- Cellratio[Cellratio$Var2=="Tumor",]

Tumor$percent <- paste0(format(Tumor$Freq*100,digits = 2),"%")

pie(x = Tumor$Freq,labels = Tumor$percent,main = "Tumor",col = mycolor)

Normal<- Cellratio[Cellratio$Var2=="Normal",]

Normal$percent <- paste0(format(Normal$Freq*100,digits = 2),"%")
pie(x = Normal$Freq,labels = Normal$percent,main = "Normal",col =mycolor)


#############
DimPlot(Bcell,split.by = "Tissue",cols =mycolor,pt.size = 0.5)
Idents(Bcell) <- Bcell$Tissue
table(Idents(Bcell), Bcell$Maincelltype)
Cellratio <- prop.table(table(Idents(Bcell), Bcell$Maincelltype), margin = 2)
Cellratio <- as.data.frame(Cellratio)

Plasma <- Cellratio[Cellratio$Var2=="Plasma cells",]
Plasma$percent <- paste0(format(Plasma$Freq*100,digits = 2),"%")
pie(x = Plasma$Freq,labels = Plasma$percent,main = "Plasma cells",col = c("#885F9C","#6980AC"))


BB <- Cellratio[Cellratio$Var2=="B cells",]
BB$percent <- paste0(format(BB$Freq*100,digits = 2),"%")
pie(x = BB$Freq,labels = BB$percent,main = "B cells",col = c("#885F9C","#6980AC"))

p1=FeaturePlot(Bcell,features = c("MS4A1"),cols =c("lightgrey","#E41A1C"),pt.size = 0.5)+theme(axis.line = element_line(size = 1))
p2=FeaturePlot(Bcell,features = c("MZB1"),cols =c("lightgrey","#E41A1C"),pt.size = 0.5)+theme(axis.line = element_line(size = 1))

p3 <- p1|p2
ggsave(p3,width =6.24,height = 2.42,filename = "../result/test.pdf")


#################
Idents(Bcell)=Bcell$RNA_snn_res.0.3
markers <- FindAllMarkers(Bcell,only.pos = T)
markers <- markers[markers$p_val_adj<0.05,]
top10 <- markers%>%group_by(cluster)%>%top_n(n = 10,wt = avg_log2FC)


celltype=data.frame(ClusterID=0:4,
                    celltype=0:4)   

celltype[celltype$ClusterID %in% c(0),2]='B cells-TNFRSF13B'
celltype[celltype$ClusterID %in% c(1,5),2]='B cells-IGHD'
celltype[celltype$ClusterID %in% c(2),2]='Plasma cells-JCHAIN'
celltype[celltype$ClusterID %in% c(3),2]='Plasma cells-IGHG1'
celltype[celltype$ClusterID %in% c(4),2]='B cells-RGS2'

for(i in 1:nrow(celltype)){
  Bcell@meta.data[which(Bcell@meta.data$RNA_snn_res.0.3 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}

table(Bcell@meta.data$celltype)

Idents(Bcell) <- Bcell$celltype
yanse <- brewer.pal(12, 'Set3')[c(1,3,4,5,7)]

p1 <- DimPlot(Bcell,label = F,cols = yanse,pt.size = 0.5,label.size = 6)+theme(axis.line = element_line(size = 1),legend.text = element_text(size = 16))
ggsave(p1,filename = "../PDF/newFig3/dimcelltype.pdf",width = 7.42,height = 4.34)

library(MySeuratWrappers)
p <- VlnPlot(Bcell, features = genes,  
             stacked=T,pt.size=0,  
             direction = "horizontal", # 
             x.lab = '', y.lab = '',cols = yanse)+# 
  theme(axis.text.x = element_blank(),   
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12),
        title = element_text(size = 12))#

ggsave(p,filename = "../PDF/newFig3/vlnplot.pdf",width = 8.24,height = 4.85)


table(Bcell$Samples)
prop.table(table(Idents(Bcell)))
table(Idents(Bcell), Bcell$Tissue)
Cellratio <- prop.table(table(Idents(Bcell), Bcell$Tissue), margin = 2)
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio)[1]<- "celltype"
colnames(Cellratio)[2]<- "Tissue"


p<-ggplot(Cellratio,aes(x=Tissue,y=Freq,fill=celltype,stratum=celltype,alluvium=celltype))+
  geom_col(width = 0.4,color=NA)+
  geom_flow(width=0.4,alpha=0.5,knot.pos=0)+theme_bw()+theme(axis.title.x = element_blank(),axis.text.x = element_text(size = 12))+
  ylab(label = "Percentages")+scale_fill_manual(values = yanse)+theme(axis.line = element_line(size = 0.6),
                                                                      axis.text.y = element_text(size = 12),
                                                                      axis.title.y = element_text(size = 12),
                                                                      legend.text = element_text(size = 12))
ggsave(p,filename = "../PDF/0805/New_fig3/chongji.pdf",width = 6.24,height = 4.85)

#####Fig3###
#####
library(Seurat)
library(monocle)
library(RColorBrewer)
load("../data/Bcell.Rdata")
Idents(Bcell) <- Bcell$celltype
data <- as(as.matrix(Bcell@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = Bcell@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())


mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=4, relative_expr = TRUE)

###seurat
var.genes <- VariableFeatures(Bcell)
mycds <- setOrderingFilter(mycds, var.genes)
p3 <- plot_ordering_genes(mycds)

####
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
####
mycds <- orderCells(mycds)

#State###
yanse <- brewer.pal(12, 'Set3')[c(1,3,4,5,7)]
tissue_color <- c("#6C9D88","#6980AC")

plot1 <- plot_cell_trajectory(mycds, color_by = "State")+scale_color_manual(values=yanse)+
  theme(title = element_text(size = 15),
        legend.text = element_text(size = 15),
        axis.line.x = element_line(size = 1),
        axis.line.y =element_line(size = 1),
        axis.text = element_text(size = 14))+
  guides(color=guide_legend(override.aes = list(size=3)))

mytheme <- theme(title = element_text(size = 15),
                 legend.text = element_text(size = 15),
                 axis.line.x = element_line(size = 1),
                 axis.line.y =element_line(size = 1),
                 axis.text = element_text(size = 14))


plot3 <- plot_cell_trajectory(mycds, color_by = "Pseudotime")+mytheme

plot4 <- plot_cell_trajectory(mycds, color_by = "celltype")+scale_color_manual(values=yanse)+
  mytheme+guides(color=guide_legend(override.aes = list(size=3)))

plot5 <- plot_cell_trajectory(mycds, color_by = "Tissue")+scale_color_manual(values=tissue_color)+mytheme+guides(color=guide_legend(override.aes = list(size=3)))

ggsave(plot1,filename = "../PDF/monocle/mono_state.pdf",width = 6.24,height = 4.42)
ggsave(plot3,filename = "../PDF/monocle/mono_Pseudotime.pdf",width = 6.24,height = 4.42)
ggsave(plot5,filename = "../PDF/monocle/mono_tissue.pdf",width = 6.24,height = 4.42)
ggsave(plot4,filename = "../PDF/monocle/mono_celltype.pdf",width =12.24,height = 6.62)



plot_cell_trajectory(mycds, color_by = "Tissue") + facet_wrap(~Tissue, nrow = 2)

Naive <- c("IGHD","TCL1A","BANK1")
memory <- c("CD38","MZB1","MS4A1")
Pla <- c("SDC1","IGHG1","JCHAIN")

p1 <- plot_genes_in_pseudotime(mycds[Naive,], color_by = "State")+scale_color_manual(values=yanse)+mytheme+NoLegend()+xlab(label = NULL)+theme(text = element_text(size = 15))
p2 <- plot_genes_in_pseudotime(mycds[memory,], color_by = "State")+mytheme+scale_color_manual(values=yanse)+NoLegend()+ylab(label = NULL)+theme(text = element_text(size = 15))
p3 <- plot_genes_in_pseudotime(mycds[Pla,], color_by = "State")+mytheme+scale_color_manual(values=yanse)+xlab(label = NULL)+ylab(label = NULL)+theme(text = element_text(size = 15))
p4 <- p1|p2|p3
p5 <- p4+guides(color=guide_legend(override.aes = list(size=3)))

ggsave(p5,filename = "../PDF/monocle/gene_state.pdf",width =7.24,height = 5.62)

####
disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= 0.5&dispersion_empirical >= 1*dispersion_fit)
disp.genes <- as.character(disp.genes$gene_id)
diff_test <- differentialGeneTest(mycds[disp.genes,], cores = 4, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test, qval < 1e-04))

p2 = plot_pseudotime_heatmap(mycds[sig_gene_names,], num_clusters=4,
                             show_rownames=F, return_heatmap=T)



p6= plot_pseudotime_heatmap(mycds[sig_gene_names,], num_clusters=4,
                            show_rownames=T, return_heatmap=T)

##BEAM###
disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= 0.5&dispersion_empirical >= 1*dispersion_fit)
disp.genes <- as.character(disp.genes$gene_id)
mycds_sub <- mycds[disp.genes,]
plot_cell_trajectory(mycds_sub, color_by = "State")

beam_res <- BEAM(mycds_sub, branch_point = 1, cores = 8)
beam_res <- beam_res[order(beam_res$qval),]
beam_res <- beam_res[,c("gene_short_name", "pval", "qval")]
mycds_sub_beam <- mycds_sub[row.names(subset(beam_res, qval < 1e-4)),]

plot_genes_branched_heatmap(mycds_sub_beam,  branch_point = 1, num_clusters = 5, show_rownames = T)
plot_genes_branched_heatmap(mycds_sub_beam,  branch_point = 1, num_clusters = 5, show_rownames = F)



write.csv(beam_res,"../data/beam_res_point1.csv")
save(mycds_sub_beam,file = "../data/mycds_sub_beam.Rdata")

####
BEAM_res=BEAM(mycds_sub,branch_point = 1,cores = 1)
BEAM_res=BEAM_res[,c("gene_short_name","pval","qval")]
saveRDS(BEAM_res, file = "BEAM_res.rds")

#####
tmp1=plot_genes_branched_heatmap(mycds_sub[row.names(subset(BEAM_res,qval<1e-4)),],
                                 branch_point = 1,
                                 num_clusters = 5, #
                                 cores = 1,
                                 branch_labels = c("Cell fate 1", "Cell fate 2"),
                                 #hmcols = NULL, #
                                 hmcols = colorRampPalette(rev(brewer.pal(9, "PRGn")))(62),
                                 branch_colors = c("#979797", "#F05662", "#7990C8"), #pre-branch, Cell fate 1, Cell fate 2
                                 use_gene_short_name = T,
                                 show_rownames = F,
                                 return_heatmap = T #
)

pdf("../PDF/Fig4/branched_heatmap.pdf",width = 8,height = 6)
pdf("../PDF/Fig4/branched_heatmap1.pdf",width = 6,height = 8)
tmp1$ph_res
dev.off()

#######
gene_group=tmp1$annotation_row
gene_group$gene=rownames(gene_group)

library(clusterProfiler)
library(org.Hs.eg.db)
allcluster_go=data.frame()
for (i in unique(gene_group$Cluster)) {
  small_gene_group=filter(gene_group,gene_group$Cluster==i)
  df_name=bitr(small_gene_group$gene, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
  go <- enrichGO(gene         = unique(df_name$ENTREZID),
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.2,
                 readable      = TRUE)
  go_res=go@result
  if (dim(go_res)[1] != 0) {
    go_res$cluster=i
    allcluster_go=rbind(allcluster_go,go_res)
  }
}
head(allcluster_go[,c("ID","Description","qvalue","cluster")])
all_go <- allcluster_go[allcluster_go$p.adjust<0.05,]
save(all_go,file = "../result/all_go.Rdata")


C1 <- all_go%>%group_by(cluster)%>%top_n(n = 10,wt = -p.adjust)

########
p1 <- ggplot(data = C1[1:10,],
             aes(x = Count,
                 y = reorder(Description,(-log10(pvalue))),
                 fill = -log10(pvalue))) +
  scale_fill_distiller(palette = "Blues",direction = 1) +##"RdPu"
  geom_bar(stat = "identity",width = 0.8) +
  theme_bw() +
  labs(x = "Count",
       y = "",
       title = "Gene cluster1")

p2 <- ggplot(data = C1[C1$cluster=="2",],
             aes(x = Count,
                 y = reorder(Description,(-log10(pvalue))),
                 fill = -log10(pvalue))) +
  scale_fill_distiller(palette = "Blues",direction = 1) +##"RdPu"
  geom_bar(stat = "identity",width = 0.8) +
  theme_bw() +
  labs(x = "Count",
       y = "",
       title = "Gene cluster2")


p3 <- ggplot(data = C1[C1$cluster=="3",],
             aes(x = Count,
                 y = reorder(Description,(-log10(pvalue))),
                 fill = -log10(pvalue))) +
  scale_fill_distiller(palette = "Blues",direction = 1) +##"RdPu"
  geom_bar(stat = "identity",width = 0.8) +
  theme_bw() +
  labs(x = "Count",
       y = "",
       title = "Gene cluster3")




p4 <- ggplot(data = C1[C1$cluster=="4",],
             aes(x = Count,
                 y = reorder(Description,(-log10(pvalue))),
                 fill = -log10(pvalue))) +
  scale_fill_distiller(palette = "Blues",direction = 1) +##"RdPu"
  geom_bar(stat = "identity",width = 0.8) +
  theme_bw() +
  labs(x = "Count",
       y = "",
       title = "Gene cluster4")


p5 <- ggplot(data = C1[C1$cluster=="5",],
             aes(x = Count,
                 y = reorder(Description,(-log10(pvalue))),
                 fill = -log10(pvalue))) +
  scale_fill_distiller(palette = "Blues",direction = 1) +##"RdPu"
  geom_bar(stat = "identity",width = 0.8) +
  theme_bw() +
  labs(x = "Count",
       y = "",
       title = "Gene cluster5")

######Fig4##############################
load("../data/integrated.final.Rdata")
load("../data/Bcell.new.Rdata")

library(CellChat)
library(ggplot2)
library(ggalluvial)
library(Seurat)
stringAsFactors=F

Idents(integrated) <- integrated$Global_Cluster
table(integrated$Global_Cluster)

No_B <- subset(integrated,invert = TRUE,idents=c("B cells","Plasma cells"))
table(No_B$Global_Cluster)
No_B_meta <- No_B@meta.data[,c("Global_Cluster","Sub_Cluster")]
No_B_meta$New_celltype <- No_B_meta$Sub_Cluster

B_meta <- Bcell@meta.data[c("Sub_Cluster","celltype")]
colnames(B_meta) <- c("Global_Cluster","Sub_Cluster")
B_meta$New_celltype <- B_meta$Sub_Cluster

All_meta <- rbind(No_B_meta,B_meta)
integrated <- AddMetaData(integrated,metadata = All_meta)
table(integrated$New_celltype)
integrated$Sub_Cluster <- factor(integrated$Sub_Cluster,levels = c("B cells-IGHD","B cells-TNFRSF13B","B cells-RGS2",
                                                                   "Plasma cells-JCHAIN","Plasma cells-IGHG1",
                                                                   "DC-CLEC9A","DC-CD1C","Macrophages-C1QC","Macrophages-SPP1",
                                                                   "Monocytes-S100A8","Mast cells-TPSAB1","Cycling-MKI67",
                                                                   "CD4-IL7R","CD4-ANXA1","CD4-FOXP3","CD8-GZMK",
                                                                   "CD8-GNLY","CD8-XCL1","CD8-SLC4A10",
                                                                   "NK cells","Hepatocytes","Fibroblasts","Endothelial cells","Cyclings"))



Idents(integrated) <- integrated$Tissue
sce<- subset(integrated,idents = "Tumor")

cellchat <- createCellChat(sce@assays$RNA@data)
meta <- data.frame(cellType = sce$Sub_Cluster, row.names =  Cells(sce))
meta$cellType <- as.character(meta$cellType)
meta$cellType <-factor(meta$cellType,levels=c("B cells-IGHD","B cells-TNFRSF13B","B cells-RGS2",
                                              "Plasma cells-JCHAIN","Plasma cells-IGHG1",
                                              "DC-CLEC9A","DC-CD1C","Macrophages-C1QC","Macrophages-SPP1",
                                              "Monocytes-S100A8","Mast cells-TPSAB1","Cycling-MKI67",
                                              "CD4-IL7R","CD4-ANXA1","CD4-FOXP3","CD8-GZMK",
                                              "CD8-GNLY","CD8-XCL1","CD8-SLC4A10",
                                              "NK cells","Hepatocytes","Fibroblasts","Endothelial cells","Cyclings"))

cellchat <- addMeta(cellchat, meta = meta, meta.name = "cellType")
cellchat <- setIdent(cellchat, ident.use = "cellType") # set "labels" as default cell identity
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

##

CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted signaling for cell-cell communication analysis
cellchat@DB <- CellChatDB.use # set the used database in the object

#####
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
cellchat <- identifyOverExpressedGenes(cellchat) 
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)


######
#USER can set type = "truncatedMean" and trim = 0.1. 
#The function computeAveExpr can help to check the average expression of signaling genes of interest, 
#e.g, computeAveExpr(cellchat, features = c("CXCL12","CXCR4"), type =  "truncatedMean", trim = 0.1).
cellchat <- computeCommunProb(cellchat)#0.25 ###

#cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat,thresh = 0.05)
cellchat <- aggregateNet(cellchat)

cellchat_tumor <- cellchat
save(cellchat_tumor,file = "../0304_cellchat/cellchat_tumor_factor.Rdata")

#######
sce<- subset(integrated,idents = "Normal")
cellchat <- createCellChat(sce@assays$RNA@data)
meta <- data.frame(cellType = sce$Sub_Cluster, row.names =  Cells(sce))
meta$cellType <- as.character(meta$cellType)
meta$cellType <-factor(meta$cellType,levels=c("B cells-IGHD","B cells-TNFRSF13B","B cells-RGS2",
                                              "Plasma cells-JCHAIN","Plasma cells-IGHG1",
                                              "DC-CLEC9A","DC-CD1C","Macrophages-C1QC","Macrophages-SPP1",
                                              "Monocytes-S100A8","Mast cells-TPSAB1","Cycling-MKI67",
                                              "CD4-IL7R","CD4-ANXA1","CD4-FOXP3","CD8-GZMK",
                                              "CD8-GNLY","CD8-XCL1","CD8-SLC4A10",
                                              "NK cells","Hepatocytes","Fibroblasts","Endothelial cells","Cyclings"))

cellchat <- addMeta(cellchat, meta = meta, meta.name = "cellType")
cellchat <- setIdent(cellchat, ident.use = "cellType") # set "labels" as default cell identity
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

##

CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted signaling for cell-cell communication analysis
cellchat@DB <- CellChatDB.use # set the used database in the object

#####
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
cellchat <- identifyOverExpressedGenes(cellchat) 
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)


######
#USER can set type = "truncatedMean" and trim = 0.1. 
#The function computeAveExpr can help to check the average expression of signaling genes of interest, 
#e.g, computeAveExpr(cellchat, features = c("CXCL12","CXCR4"), type =  "truncatedMean", trim = 0.1).
cellchat <- computeCommunProb(cellchat)#0.25 ###

#cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat,thresh = 0.05)
cellchat <- aggregateNet(cellchat)

cellchat_normal <- cellchat
save(cellchat_normal,file = "../0304_cellchat/cellchat_normal_factor.Rdata")

###
load("../0304_cellchat/cellchat_tumor_factor.Rdata")
load("../0304_cellchat/cellchat_normal_factor.Rdata")

object.list <- list(Normal = cellchat_normal,Tumor= cellchat_tumor)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

#####
gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)

####hepatocytes vs TIL-Bs###
Bh<- netVisual_bubble(cellchat, sources.use = 1:5, targets.use = c(21),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Tumor", angle.x = 45, remove.isolate = T)
hB<- netVisual_bubble(cellchat, sources.use = 21, targets.use = c(1:5),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Tumor", angle.x = 45, remove.isolate = T)

#######T cells vs TIL-Bs####
BT<- netVisual_bubble(cellchat, sources.use = c(1:5), targets.use = c(13:19),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Tumor", angle.x = 45, remove.isolate = T)
TB<- netVisual_bubble(cellchat, sources.use = c(13:19), targets.use = c(1:5),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Tumor", angle.x = 45, remove.isolate = T)

####myeloid vs TIL-Bs#####

BM<- netVisual_bubble(cellchat, sources.use = c(1:5), targets.use = c(6:12),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Tumor", angle.x = 45, remove.isolate = T)
MB<- netVisual_bubble(cellchat, sources.use = c(6:12), targets.use = c(1:5),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Tumor", angle.x = 45, remove.isolate = T)



####FigS1####
library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggalluvial)
library(RColorBrewer)
load("../data/integrated.final.Rdata")
set.seed(123)

#######T cells
Tcell <- subset(integrated,idents = c("T cells")) %>%
  FindNeighbors(dims = 1:30, verbose = FALSE, reduction = "harmony") %>%
  FindClusters(resolution = 0.8,verbose = FALSE) %>%
  RunUMAP(dims = 1:30, verbose = FALSE, reduction = "harmony")

Idents(Tcell) <- Tcell$Sub_Cluster

yanse <- brewer.pal(12, 'Set3')[c(1,3,4,5,6,7,10)]
p1 <- DimPlot(Tcell,cols = yanse,label.size = 4,pt.size = 0.005)+theme(axis.line = element_line(size = 0.5),legend.text = element_text(size = 12))

ggsave(p1,filename = "../PDF/Tcell/Dimplot_new.pdf",width = 4.84,height = 3.85)

###
table(Tcell$Group)
prop.table(table(Idents(Tcell)))
table(Idents(Tcell), Tcell$Tissue)
Cellratio <- prop.table(table(Idents(Tcell), Tcell$Tissue), margin = 2)
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio)[1]<- "celltype"
colnames(Cellratio)[2]<- "Tissue"


p2<-ggplot(Cellratio,aes(x=Tissue,y=Freq,fill=celltype,stratum=celltype,alluvium=celltype))+
  geom_col(width = 0.4,color=NA)+
  geom_flow(width=0.4,alpha=0.5,knot.pos=0)+theme_bw()+theme(axis.title.x = element_blank(),axis.text.x = element_text(size = 15))+
  ylab(label = "Percentages")+scale_fill_manual(values = yanse)+theme(axis.line = element_line(size = 0.6),
                                                                      axis.text.y = element_text(size = 15),
                                                                      axis.title.y = element_text(size = 15),
                                                                      legend.text = element_text(size = 15))

ggsave(p2,filename = "../PDF/Tcell/ggalluvial.pdf",width = 6.24,height = 3.85)


table(Idents(Tcell), Tcell$Group)
Cellratio <- prop.table(table(Idents(Tcell), Tcell$Group), margin = 2)
Cellratio <- as.data.frame(Cellratio)
Cellratio <- Cellratio[Cellratio$Var1=="CD4-FOXP3",]

Cellratio$Tissue <- rep(c("Normal","Tumor"),4)

mycolor <- c("#6C9D88","#6980AC")

p3 <- ggviolin(Cellratio,x = "Tissue",y = "Freq",fill = "Tissue",palette = mycolor,add = "boxplot",add.params =list(fill = "white"))+
  stat_compare_means(method = "wilcox.test",paired = F, 
                     comparisons=list(c("Tumor", "Normal")),label = "p.signif",vjust = 0.4,label.y = 0.3,size=10)+
  NoLegend()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 20),
        axis.title.y =element_text(size = 20))+
  ylab(label = "Percentages")

ggsave(p3,filename = "../PDF/Tcell/spp1.pdf",width = 4.54,height = 5.25)


#####
stimulatory <- c("ICOS","TNFRSF4","TNFRSF9","TNFRSF8",'TNFRSF18',
                 "CD40LG","CD27","CD28","CD226",
                 "TNFSF14")

inhibitory<- c("PDCD1",'HAVCR2','LAIR1','CTLA4','TIGIT','CD160','BTLA')

######
response <- c("LAG3",'CD244','CCL5','CXCL9','CXCL13','MSR1','CSF3R','CYBB','KLRK1')

p1 <- DotPlot(Tcell,features = inhibitory,cols = "RdYlBu")+
  theme(axis.text.x = element_text(hjust = 1,vjust = 1,angle = 45),axis.title.x = element_blank(),axis.title.y = element_blank())

p2 <- DotPlot(Tcell,features = stimulatory,cols = "RdYlBu")+
  theme(axis.text.x = element_text(hjust = 1,vjust = 1,angle = 45),axis.title.x = element_blank(),axis.title.y = element_blank())

ggsave(p1,file="../PDF/Tcell/inhibit.pdf",width = 5.84,height = 4.84)
ggsave(p2,file="../PDF/Tcell/stimulatory.pdf",width = 6.34,height = 4.84)

##
Genes<- c("CD4","CD8A",#conventional
          "GZMK","LEF1","LTB","SELL","IL7R","TCF7","CCR7",#naive
          "IFNG","GZMB","GZMH","PRF1","NKG7","GNLY",#Effector
          "ENPD1","HAVCR2","CTLA4","TNFRSF9","TOX2","TIGIT","LAG3",#Exhausted 
          "CD69",'NR4A1','MYADM','GATA3','TBX21','ANXA1',#Memory
          'SLC4A10','NCR3','RORC','RORA',#mait
          "FOXP3",'IL2RA')#Treg 
x <- as.data.frame(AverageExpression(Tcell,features = Genes))
library(pheatmap)


Celltype<- c('CD4-IL7R','CD4-ANXA1','CD4-FOXP3','CD8-GZMK','CD8-GNLY','CD8-XCL1','CD8-SLC4A10')
Celltype <- as.data.frame(Celltype)
rownames(Celltype) <- colnames(x)

mat_colors <- list(Celltype =yanse)
names(mat_colors$Celltype) <- unique(Celltype[,1])

p <- pheatmap(x,cluster_rows = F,cluster_cols = F,scale = "row",
              annotation_colors = mat_colors,show_colnames = F,
              annotation_col = Celltype,
              color = colorRampPalette(c("navy", "white", "firebrick3"))(100))

ggsave(p,file="../PDF/Tcell/heatmap.pdf",width = 5.74,height = 6.84)




#######

table(sce$orig.ident)#
prop.table(table(Idents(sce)))
table(Idents(sce), sce$orig.ident)#

Cellratio <- prop.table(table(Idents(sce), sce$orig.ident), margin = 2)#
Cellratio <- as.data.frame(Cellratio)

colnames(Cellratio)[1]<- "celltype"
colnames(Cellratio)[2]<- "Sample"

Cellratio$group <- ifelse(substr(Cellratio$Sample,6,6)=="N","Normal","Tumor")


library(ggplot2)
Cellratio$Samples <- factor(Cellratio$Sample,levels = c("HCC01T","HCC02T","HCC03T","HCC04T","HCC05T","HCC06T","HCC07T","HCC08T","HCC09T","HCC10T","HCC03N","HCC04N","HCC05N","HCC06N","HCC07N","HCC08N","HCC09N","HCC10N"))

p1 <- ggplot(Cellratio) + 
  geom_bar(aes(x =Samples, y= Freq, fill = celltype),stat = "identity",width = 0.9,size = 0.2)+
  theme_classic() +
  labs(x=NULL,y = 'Ratio of celltype')+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size = 20),
                                             axis.text = element_text(size = 20, color = 'black'),
                                             axis.line = element_line(colour = 'black', size = 2),
                                             axis.title = element_text(size = 20, color = 'black'),
                                             text = element_text(size = 20))
####FigS2######
library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggalluvial)
library(RColorBrewer)
load("../data/integrated.final.Rdata")

Myeloids<- subset(integrated,idents = c("Myeloids")) %>%
  FindNeighbors(dims = 1:30, verbose = FALSE, reduction = "harmony") %>%
  FindClusters(resolution = 0.5,verbose = FALSE) %>%
  RunUMAP(dims = 1:30, verbose = FALSE, reduction = "harmony")

Idents(Myeloids) <- Myeloids$Sub_Cluster

yanse <- brewer.pal(12, 'Set3')[c(1,3,4,5,6,7,10)]
p1 <- DimPlot(Myeloids,cols = yanse,pt.size = 0.05,label.size = 6)+theme(axis.line = element_line(size = 1),legend.text = element_text(size = 16))

ggsave(p1,filename = "../PDF/Myeloid/Dimplot_new.pdf",width = 5.24,height = 2.85)

####
table(Myeloids$Group)
prop.table(table(Idents(Myeloids)))
table(Idents(Myeloids), Myeloids$Tissue)
Cellratio <- prop.table(table(Idents(Myeloids), Myeloids$Tissue), margin = 2)
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio)[1]<- "celltype"
colnames(Cellratio)[2]<- "Tissue"


p2<-ggplot(Cellratio,aes(x=Tissue,y=Freq,fill=celltype,stratum=celltype,alluvium=celltype))+
  geom_col(width = 0.4,color=NA)+
  geom_flow(width=0.4,alpha=0.5,knot.pos=0)+theme_bw()+theme(axis.title.x = element_blank(),axis.text.x = element_text(size = 15))+
  ylab(label = "Percentages")+scale_fill_manual(values = yanse)+theme(axis.line = element_line(size = 0.6),
                                                                      axis.text.y = element_text(size = 15),
                                                                      axis.title.y = element_text(size = 15),
                                                                      legend.text = element_text(size = 15))

ggsave(p2,filename = "../PDF/Myeloid/ggalluvial.pdf",width = 6.24,height = 3.85)


table(Idents(Myeloids), Myeloids$Group)
Cellratio <- prop.table(table(Idents(Myeloids), Myeloids$Group), margin = 2)
Cellratio <- as.data.frame(Cellratio)
Cellratio <- Cellratio[Cellratio$Var1=="Macrophages-SPP1",]

Cellratio$Tissue <- rep(c("Normal","Tumor"),4)

mycolor <- c("#6C9D88","#6980AC")

p3 <- ggviolin(Cellratio,x = "Tissue",y = "Freq",fill = "Tissue",palette = mycolor,add = "boxplot",add.params =list(fill = "white"))+
  stat_compare_means(method = "t.test",paired = F, 
                     comparisons=list(c("Tumor", "Normal")),label = "p.signif",vjust = 0.4,label.y = 0.4,size=10)+
  NoLegend()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 20),
        axis.title.y =element_text(size = 20))+
  ylab(label = "Percentages")

ggsave(p3,filename = "../PDF/Myeloid/spp1.pdf",width = 4.54,height = 5.25)


genes <- c("LYZ","S100A8","S100A9","MKI67","C1QC",'CXCL10','SPP1','VEGFA','MRC1','CLEC9A','XCR1',
           'CD1C','CD1D','ITGAX','CSF3R','IL1B','CPA3','KIT')

library(MySeuratWrappers)
p <- VlnPlot(Myeloids, features = genes,  
             stacked=T,pt.size=0,  
             direction = "horizontal", #  
             x.lab = '', y.lab = '',cols = yanse)+#  
  theme(axis.text.x = element_blank(),   
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12),
        title = element_text(size = 15))# 

ggsave(p,filename = "../PDF/Myeloid/vlnplot.pdf",width = 8.24,height = 4.85)

#######TCGA-LIHC###
###TCGA#####
library(tidyverse)
library(openxlsx)
library(estimate)
#1）
data <- read.table("../TCGA-LIHC.htseq_fpkm.tsv.gz",sep = "\t" , header = T,
                   #row.names = "Ensembl_ID",
                   stringsAsFactors = FALSE ,
                   check.names = FALSE)
#2）
probeMap <- read.table("../gencode.v22.annotation.gene.probeMap",sep = "\t" , header = T,
                       stringsAsFactors = FALSE ,
                       check.names = FALSE)
# 3）
exp <- data %>%
  inner_join(probeMap, by = c("Ensembl_ID" = "id")) %>%
  select(gene , starts_with("TCGA") ) 
#4）
expr <- exp %>%
  mutate(rowMean =rowMeans(.[grep("TCGA", names(.))])) %>%
  arrange(desc(rowMean)) %>%
  distinct(gene,.keep_all = T) %>%
  select(-rowMean) %>%
  column_to_rownames(var = "gene")
expr[1:4,1:4]

phenotype <- read.delim("../TCGA-LIHC.GDC_phenotype.tsv/TCGA-LIHC.GDC_phenotype.tsv")
xx <- phenotype[,c("submitter_id.samples",'radiation_therapy','history_prior_surgery_type_other','history_of_neoadjuvant_treatment','sample_type.samples','tumor_stage.diagnoses')]
xx$submitter_id.samples <- gsub("-",".",xx$submitter_id.samples,fixed = T)
xx <- xx[substr(xx$submitter_id.samples,14,15)<10,]
xx <- xx[xx$history_of_neoadjuvant_treatment=="No",]
xx <- xx[xx$sample_type.samples=="Primary Tumor",]
xx <- xx[xx$radiation_therapy!="YES",]


b <- read.table("../UCSC-LIHC.tpm.txt",header = T)
expr <- b
rownames(expr) <- expr$gene
expr[1] <-NULL 
expr[1:5,1:5]
expr <- expr[,substr(colnames(expr),14,15)<10]

TCGA_clin <- read.delim("../TCGA-LIHC.survival.tsv")
TCGA_clin$sample <- gsub("-",".",TCGA_clin$sample,fixed = T)
expr <- expr[,colnames(expr)%in%TCGA_clin$sample]
expr <- expr[,colnames(expr)%in%xx$submitter_id.samples]

#########fpkm##########
fpkm <- 2^expr - 1
#
fpkm <- fpkm %>% 
  rownames_to_column("gene")
#
write.table(fpkm,"../UCSC-LIHC.fpkm.txt",sep = "\t" ,quote = F, row.names=F)
write.table(a,"../a.txt",sep = "\t" ,quote = F, row.names=F)

###
expr <- read.table("../UCSC-LIHC.fpkm.txt",header = T)

fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

tpms <- apply(expr[,2:425],2,fpkmToTpm)
tpms <- as.data.frame(tpms)
tpms <- cbind(expr$gene,tpms)
colnames(tpms)[1] <- "gene"
write.table(tpms,"../UCSC-LIHC.tpm.txt",sep = "\t" ,quote = F, row.names=F)

######Fig5A-C#######
expr <- write.table("../UCSC-LIHC.tpm.txt")
library(CIBERSORT)
expr <- as.matrix(expr)
res1 <- cibersort(sig_matrix = LM22, mixture_file = expr,perm = 500)
res1 <- as.data.frame(res1)
res1$sample <- rownames(res1)

dat1 <- merge(TCGA_clin,res1,by = "sample")
dat1$naive <- dat1$`B cells naive`
dat1$memory <- dat1$`B cells memory`
dat1$Plasma <- dat1$`Plasma cells`
dat1$OS.time <- dat1$OS.time/30

library(survival)
library(survminer)
res.cut1 <- surv_cutpoint(dat1, time = "OS.time", event = "OS",
                          variables = c("naive","memory","Plasma"))                                  
summary(res.cut1)
res.cat1 <- surv_categorize(res.cut1)
head(res.cat1)
table(res.cat1$Plasma)

fit1 <- survfit(Surv(OS.time, OS) ~Plasma, data = res.cat1)
ggsurvplot(fit1, main= "Plasma cells",data = res.cat1,palette= c("red","blue"),risk.table =T, 
           conf.int = F,pval = T,legend=c(0.8,0.8))+xlab(label = "Months")

fit2 <- survfit(Surv(OS.time, OS) ~naive, data = res.cat1)
ggsurvplot(fit2, main= "naive cells",data = res.cat1,palette= c("red","blue"),risk.table =T, 
           conf.int = F,pval = T,legend=c(0.8,0.8))+xlab(label = "Months")


fit3 <- survfit(Surv(OS.time, OS) ~memory, data = res.cat1)
ggsurvplot(fit3, main= "memory cells",data = res.cat1,palette= c("red","blue"),risk.table =T, 
           conf.int = F,pval = T,legend=c(0.8,0.8))+xlab(label = "Months")

###Fig6 A-C####
library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggalluvial)
library(RColorBrewer)
load("../data/integrated.final.Rdata")
markers <- FindAllMarkers(integrated,min.pct = 0.2,only.pos = T)
markers1 <-markers[markers$p_val_adj<0.05,]
plasma_markers <- markers1[markers1$cluster=="Plasma cells",]
save(markers1,plasma_markers,file = "../0810/deg_all.Rdata")

#####UNICOX####
library(KEGGREST)
library(GSVA)
library(survival)
library(survminer)
library(limma)
expr <- readRDS("../data/expr_356samples.Rdata")
expr[1:5,1:5]
clinical <- read.delim("../data/UCSC-LIHC.survival.tsv")
clinical$OS.time <- clinical$OS.time/30
colnames(clinical)[1] <- "id"

expr <- expr[,substr(colnames(expr),14,15)<10]
expr <- expr[,!duplicated(substr(colnames(expr),1,12))]
load("../0810/deg_all.Rdata")
expr<- expr[rownames(expr)%in%plasma_markers$gene,]

expr <- as.data.frame(t(expr))

rownames(expr)<- gsub(".","-",rownames(expr),fixed = T)
expr$id <- rownames(expr)

td <- merge(clinical,expr,by = "id")
rownames(td) <- td$id
td[1] <- NULL

pFilter=0.05 #
outResult=data.frame() #
sigGenes=c("OS","OS.time") #

for(i in colnames(td[,4:ncol(td)])){ #
  tdcox <- coxph(Surv(OS.time, OS) ~ td[,i], data = td)#
  tdcoxSummary = summary(tdcox) #
  pvalue=tdcoxSummary$coefficients[,"Pr(>|z|)"] #
  if(pvalue<pFilter){ # 
    sigGenes=c(sigGenes,i)
    outResult=rbind(outResult,#
                    cbind(id=i,#
                          HR=tdcoxSummary$conf.int[,"exp(coef)"],#
                          L95CI=tdcoxSummary$conf.int[,"lower .95"],#
                          H95CI=tdcoxSummary$conf.int[,"upper .95"],#
                          pvalue=tdcoxSummary$coefficients[,"Pr(>|z|)"])#
    )
  }
}

write.table(outResult,file="../0810/UniCoxSurvival.txt",sep="\t",row.names=F,quote=F)

UniCoxSurSigGeneExp=td[,sigGenes] #
UniCoxSurSigGeneExp=cbind(id=row.names(UniCoxSurSigGeneExp),UniCoxSurSigGeneExp)#以id也就是样品名命名行名
write.table(UniCoxSurSigGeneExp,file="../0810/UniCoxSurSigGeneExp.txt",sep="\t",row.names=F,quote=F)

#####LASSO####
library(glmnet)
library(caret)
library(survival)
library(survminer)

outResult <- read.table("../0810/UniCoxSurvival.txt",header = T)
UniCoxSurSigGeneExp <- read.table("../0810/UniCoxSurSigGeneExp.txt",header = T)

x<- UniCoxSurSigGeneExp[,4:64]
x <- x[,outResult[outResult$HR>1,]$id]
rownames(x) <- UniCoxSurSigGeneExp$id
x <- as.matrix(x)

####
y <- as.matrix(UniCoxSurSigGeneExp[,2])
rownames(y) <- UniCoxSurSigGeneExp$id
identical(rownames(x),rownames(y))
####
library(glmnet)
model_lasso <- glmnet(x,y,alpha = 1,family = "binomial")
print(model_lasso)
plot(model_lasso)
plot(model_lasso, xvar = "lambda", label = TRUE)

set.seed(123)
cv_fit <- cv.glmnet(x=x, y=y,alpha = 1)
plot(cv_fit)

#### 2.2 
model_lasso_min <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.min)
model_lasso_1se <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.1se)


head(model_lasso_min$beta)
choose_gene_min=rownames(model_lasso_min$beta)[as.numeric(model_lasso_min$beta)!=0]
choose_gene_1se=rownames(model_lasso_1se$beta)[as.numeric(model_lasso_1se$beta)!=0]
length(choose_gene_min)
length(choose_gene_1se)

lasso.prob <- predict(cv_fit,  #predict
                      newx=x ,  #
                      s=c(cv_fit$lambda.min,cv_fit$lambda.1se) )
re=cbind(y ,lasso.prob)  #cbind
head(re)

##

re=as.data.frame(re)
colnames(re)=c('event','prob_min','prob_1se')
#####
re$event=as.factor(re$event)  #
library(ggpubr) 
p1 = ggboxplot(re, x = "event", y = "prob_min",
               color = "event", palette = "jco",
               add = "jitter")+ stat_compare_means()
p2 = ggboxplot(re, x = "event", y = "prob_1se",
               color = "event", palette = "jco",
               add = "jitter")+ stat_compare_means()
library(patchwork)
p1+p2

##
#### 3.3####

library(ROCR)
library(caret)

pred_min <- prediction(re[,2], re[,1])
auc_min = performance(pred_min,"auc")@y.values[[1]]
perf_min <- performance(pred_min,"tpr","fpr")
plot(perf_min,colorize=FALSE, col="blue") 
lines(c(0,1),c(0,1),col = "gray", lty = 4 )
text(0.8,0.2, labels = paste0("AUC = ",round(auc_min,2)))

#1se
pred_1se <- prediction(re[,3], re[,1])
auc_1se = performance(pred_1se,"auc")@y.values[[1]]
perf_1se <- performance(pred_1se,"tpr","fpr")
plot(perf_1se,colorize=FALSE, col="red") 
lines(c(0,1),c(0,1),col = "gray", lty = 4 )
text(0.8,0.2, labels = paste0("AUC = ",round(auc_1se,2)))

#
plot(perf_min,colorize=FALSE, col="blue") 
plot(perf_1se,colorize=FALSE, col="red",add = T) 
lines(c(0,1),c(0,1),col = "gray", lty = 4 )
text(0.8,0.3, labels = paste0("AUC = ",round(auc_min,3)),col = "blue")
text(0.8,0.2, labels = paste0("AUC = ",round(auc_1se,3)),col = "red")
#

library(ROCR)
library(caret)
# 
#min
pred_min <- prediction(re[,2], re[,1])
auc_min = performance(pred_min,"auc")@y.values[[1]]
perf_min <- performance(pred_min,"tpr","fpr")

#1se
pred_1se <- prediction(re[,3], re[,1])
auc_1se = performance(pred_1se,"auc")@y.values[[1]]
perf_1se <- performance(pred_1se,"tpr","fpr")

tpr_min = performance(pred_min,"tpr")@y.values[[1]]
tpr_1se = performance(pred_1se,"tpr")@y.values[[1]]
dat = data.frame(tpr_min = perf_min@y.values[[1]],
                 fpr_min = perf_min@x.values[[1]],
                 tpr_1se = perf_1se@y.values[[1]],
                 fpr_1se = perf_1se@x.values[[1]])

ggplot()+ 
  geom_line(data = dat,aes(x = fpr_min, y = tpr_min),color = "blue") + 
  geom_line(data = dat,aes(x = fpr_1se, y = tpr_1se),color = "red")+
  geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey")+
  theme_bw()+
  annotate("text",x = .75, y = .25,
           label = paste("AUC of min = ",round(auc_min,2)),color = "blue")+
  annotate("text",x = .75, y = .15,label = paste("AUC of 1se = ",round(auc_1se,2)),color = "red")+
  scale_x_continuous(name  = "fpr")+
  scale_y_continuous(name = "tpr")

########
sig_gene_mult_cox <- choose_gene_min
DEG_met_expr.lasso_cox <- UniCoxSurSigGeneExp[,c("OS","OS.time",sig_gene_mult_cox)]
rownames(DEG_met_expr.lasso_cox) <- UniCoxSurSigGeneExp$id

multiCox <- coxph(Surv(OS.time, OS) ~ ., data =  DEG_met_expr.lasso_cox)
summary(multiCox)

#predict
riskScore=apply(DEG_met_expr.lasso_cox[,names(multiCox$coefficients)], 1,function(k){
  sum(multiCox$coefficients * k)
}) 

riskScore<-as.data.frame(riskScore)
riskScore$sample <- rownames(riskScore)
head(riskScore,2)

riskScore_cli <- riskScore
riskScore_cli$OS <- DEG_met_expr.lasso_cox$OS
riskScore_cli$OS.time <- DEG_met_expr.lasso_cox$OS.time
riskScore_cli$riskScore2 <- ifelse(riskScore_cli$riskScore > median(riskScore_cli$riskScore),
                                   "High","Low")

#KM
fit <- survfit(Surv(OS.time, as.numeric(OS)) ~ riskScore2, data=riskScore_cli)

lasso_KM <- ggsurvplot(fit, data = riskScore_cli,
                       pval = T,
                       risk.table = T,
                       #surv.median.line = "hv", #
                       palette=c("red", "blue"),  #
                       #legend.labs=c("High risk","Low risk"), #
                       legend.title="RiskScore",
                       title="Overall survival", #
                       ylab="Cumulative survival (percentage)",xlab = " Time (Months)", #更改横纵坐标
                       censor.shape = 124,censor.size = 2,conf.int = FALSE)
lasso_KM

####ROC###
library(timeROC)
with(riskScore_cli,
     ROC_riskscore <<- timeROC(T = OS.time,
                               delta = OS,
                               marker = riskScore,
                               cause = 1,
                               weighting = "marginal",
                               times = c(12,24,36),
                               ROC = TRUE,
                               iid = TRUE))

plot(ROC_riskscore, time = 12, col = "red", add = F,title = "")
plot(ROC_riskscore, time = 24, col = "blue", add = T)
plot(ROC_riskscore, time = 36, col = "purple", add = T)
legend("bottomright",c("1-Year","3-Year","5-Year"),col=c("red","blue","purple"),lty=1,lwd=2)
text(0.5,0.2,paste("1-Year AUC = ",round(ROC_riskscore$AUC[1],3)))
text(0.5,0.15,paste("2-Year AUC = ",round(ROC_riskscore$AUC[2],3)))
text(0.5,0.1,paste("3-Year AUC = ",round(ROC_riskscore$AUC[3],3)))

####Fig7#####
expr <- readRDS("../data/expr_356samples.Rdata")
expr[1:5,1:5]
genes <- c("PDCD1",'CD274','PDCD1LG2','HAVCR2','LGALS9','LAIR1','PTPN6','PTPN11',
           'CTLA4','TIGIT','BTLA','CD160')

dat <- expr[rownames(expr)%in%genes,]
dat <- as.data.frame(t(dat))

dat$sample <- gsub(x=rownames(dat),pattern = '.',replacement ="-",fixed = T)
dat <-merge(riskScore_cli,dat,by = "sample")
colnames(dat)[5] <- "Group"
group <- dat$Group
rownames(dat) <- dat$sample
dat <- dat[,-c(1:5)]

xx <- melt(dat)
xx$group <- rep(group,12)

ggplot(xx, aes(x = variable, y = log2(value+1)))+
  geom_boxplot(aes(fill=group))+ 
  stat_compare_means(aes(group = group),
                     label = "p.signif", 
                     method = "wilcox.test",
                     hide.ns =F)+theme_bw()+theme(axis.text.x = element_text(hjust = 1,vjust = 1,angle = 45))+scale_fill_manual(values = c("#FB8072","#80B1D3"))
coord_flip()

######TMB###
setwd("../TMB/gdc_download_20230823_071234.693255/")
##
files <- list.files(pattern = '*.gz',recursive = TRUE)
all_mut <- data.frame()
for (file in files) {
  mut <- read.delim(file,skip = 7, header = T, fill = TRUE,sep = "\t")
  all_mut <- rbind(all_mut,mut)
}


###
a <- all_mut
write.csv(a,file = "../a.csv")

##TMB##
get_TMB <- function(file) {
  # 
  use_cols <- c(
    "Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode", 
    "HGVSc", "t_depth", "t_alt_count"
  )
  # 
  mut_type <- c(
    "5'UTR", "3'UTR", "3'Flank", "5'Flank", "Intron", "IGR",
    "Silent", "RNA", "Splice_Region"
  )
  # 
  df <- read_csv(file, col_select = use_cols)
  data <- df %>% subset(!Variant_Classification %in% mut_type) %>%
    # 
    mutate(vaf = t_alt_count / t_depth) %>%
    group_by(Tumor_Sample_Barcode) %>%
    summarise(mut_num = n(), TMB = mut_num / 30, MaxVAF = max(vaf))
  return(data)
}

lihc.csv <- get_TMB('../a.csv')
write.csv(lihc.csv,file = "../lihc.csv")

####################
a <- read.csv("../0810/TMB/a.csv")
load("../0810/model.Rdata")

a$Tumor_Sample_Barcode <- substr(a$Tumor_Sample_Barcode,1,16)
a$Matched_Norm_Sample_Barcode <- substr(a$Matched_Norm_Sample_Barcode,1,16)

High <- riskScore_cli[riskScore_cli$riskScore2=="High",]
High_dat <- a[a$Tumor_Sample_Barcode%in%High$sample,]

Low <- riskScore_cli[riskScore_cli$riskScore2=="Low",]
Low_dat <- a[a$Tumor_Sample_Barcode%in%Low$sample,]
save(High,High_dat,Low,Low_dat,file = "../0810/TMB/Group_oncoplots.Rdata")

all <- a[a$Tumor_Sample_Barcode%in%riskScore_cli$sample,]
########
library(maftools)
Hi_mut <- read.maf(High_dat)
Low_mut <- read.maf(Low_dat)
all_mut <- read.maf(all)

plotmafSummary(maf =Low_mut, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

###oncoplots##
vc_cols <- RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) <- c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)
oncoplot(maf = all_mut, colors = vc_cols, top = 20)
oncoplot(maf = Low_mut, colors = vc_cols,genes =all_mut@gene.summary$Hugo_Symbol[1:20],keepGeneOrder = T)
oncoplot(maf = Hi_mut, colors = vc_cols,genes =all_mut@gene.summary$Hugo_Symbol[1:20],keepGeneOrder = T)

#######
lihc.csv <- read.csv("../0810/TMB/lihc.csv")
lihc.csv$Tumor_Sample_Barcode <- substr(lihc.csv$Tumor_Sample_Barcode,1,16)
table(duplicated(lihc.csv$Tumor_Sample_Barcode))
colnames(lihc.csv)[2] <- "sample"

load("../0810/model.Rdata")
expr1 <- merge(riskScore_cli,lihc.csv,by = "sample")
expr1 <- expr1[,c(5,8)]
colnames(expr1)[1] <- "Group"

ks.test(expr1$TMB,"pnorm")#####
library(ggplot2)
library(ggpubr)
library(ggsci)

ggplot(expr1,aes(x = Group,y = log2(TMB+1),fill=Group))+geom_boxplot(color=c('red','blue'))+
  stat_compare_means(aes(group = Group),
                     label = "p.signif", 
                     method = "wilcox.test",
                     hide.ns = F)+theme_bw()
##TIDE###
expr <-  read.table("../CIBERSORT/UCSC-LIHC.tpm.txt",header = T)
rownames(expr) <- expr$gene
expr[1] <- NULL
expr <- expr[,substr(colnames(expr),14,15)<10]
expr =expr[,!duplicated(substr(colnames(expr),1,12))]

mydat<-t(apply(expr,1,function(x)x-(mean(x))))
head(mydat[,1:5])
write.table(mydat,file='../result/UCSC_TIDE.txt',sep='\t')##

#####
patient_ICBresponse <- read.csv("../result/UCSC_LIHC_TIDEres.csv",header = TRUE)
ICBresponse <- patient_ICBresponse
##
table(ICBresponse$Responder=="False")

ICBresponse$Responder <- ifelse(str_detect(ICBresponse$Responder,"False"),"NR","R")
ICBresponse1 <- dplyr::pull(ICBresponse, 3) ##
names(ICBresponse1) <- ICBresponse$Patient
ICBresponse1[1:10]
save(ICBresponse,file = "ICBresponse.Rdata")


load("ICBresponse.Rdata")
##
## barplot
dat_plot <- data.frame(id = ICBresponse$Patient,
                       t = ICBresponse$TIDE)

# 
dat_plot$threshold = factor(ifelse(dat_plot$t > 0,'NR','R'),levels=c('R','NR'))
# 
dat_plot <- dat_plot %>% arrange(t)
#
dat_plot$id <- factor(dat_plot$id,levels = dat_plot$id)
#
library(ggplot2)
library(ggthemes)
library(ggprism)
p5 <- ggplot(data = dat_plot,aes(x = id,y = t,fill = threshold)) +
  geom_col()+
  # scale_x_continuous(breaks = seq(0,50,100))+ # X 
  scale_fill_manual(values = c('R'= '#B44C4C','NR'='#3478AC')) +
  geom_hline(yintercept = c(-0.5,0.5),color = 'white',size = 0.5,lty='dashed') +
  xlab('') + 
  ylab('TIDE score') + 
  guides(fill=guide_legend(key.width = 3, key.height = 5, nrow = 2, ncol = 1, byrow = TRUE))+ #显示图例
  theme_bw()+
  theme_prism(border = T) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = c(.95, .95),#
    legend.justification = c("right", "top")#
  )
p5

##


##
rescat <- readRDS("../result/UCSC_rescat.Rdata")
rescat$id <- rownames(rescat)
rescat$id <- gsub(pattern = "-",replacement = ".",rescat$id,fixed = T)

dat_plot <- merge(dat_plot,rescat,by = "id")

table(dat_plot$Plasma[dat_plot$threshold=="R"])
table(dat_plot$Plasma[dat_plot$threshold=="NR"])

group <- c("Plasma-high","Plasma-high","Plasma-low","Plasma-low")
Response <- c("NR","R","NR","R")
Percent <- c(0.4679,0.5321,0.3400,0.6600)
num <- c(124,141,34,66)
data <- data.frame(group,Response,Percent,num)


##
NR <- c(124,34)
R <- c(141,66)
dat <- data.frame(NR,R)
rownames(dat) <- c("Plasma-high","Plasma-low")
chisq.test(dat)
# X-squared = 24.263, df = 1, p-value = 8.403e-07

p6 = ggplot(data, aes( x = group, weight = Percent, fill = Response))+
  geom_bar( position = "stack")+xlab("X-squared = 24.263, df = 1, p-value = 0.03739")

p7 <- p6+theme(panel.grid = element_blank())+theme_bw()+
  scale_fill_manual( values = c('#B44C4C','#3478AC'))
p7

##
library(patchwork)
p <- (p5 +p7+p3)+plot_annotation(tag_levels = 'A')
p
ggsave(p,filename = "../rawdata/immune-response.pdf",width = 9,height = 5)

#####Fig8###
options(stringsAsFactors = F)
library(oncoPredict)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)
th=theme(axis.text.x = element_text(angle = 45,vjust = 0.5))
dir='../815/DataFiles/DataFiles/Training Data/'

GDSC2_Expr = readRDS(file=file.path(dir,'GDSC2_Expr (RMA Normalized and Log Transformed).rds'))
GDSC2_Res = readRDS(file = file.path(dir,"GDSC2_Res.rds"))
GDSC2_Res <- exp(GDSC2_Res)#

expr1 <- readRDS("../data/expr_356samples.Rdata")
expr1 <- log2(expr1+1)
testExpr <- as.matrix(expr1)

###
setwd("../815/")
calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = testExpr,
              batchCorrect = 'eb',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData' )

#######
library(data.table)
library(tidyverse)
testPtype <- fread("../815/TCGA_calcPhenotype_Output/DrugPredictions.csv", data.table = F)
testPtype[1:4, 1:4]
load("../0810/model.Rdata")
x <- column_to_rownames(testPtype,var = "V1")
x <- x[,c('Cisplatin_1005','Oxaliplatin_1089','Oxaliplatin_1806','5-Fluorouracil_1073',
          'Paclitaxel_1080','Irinotecan_1088','Vincristine_1818','Sorafenib_1085')]

riskScore_cli$sample <- gsub(pattern = '-',replacement = ".",x =riskScore_cli$sample)
x$sample <- rownames(x)
xx <- merge(x,riskScore_cli,by = "sample")
yy <- xx
group <- xx$riskScore2
colnames(xx)
xx <- xx[,1:9]
xx <- column_to_rownames(xx,var = "sample")

xx <- melt(xx)
xx$group <- rep(group,8)


ggplot(xx, aes(x = variable, y =log2(value+1)))+
  geom_boxplot(aes(fill=group))+ 
  stat_compare_means(aes(group = group),
                     label = "p.signif", 
                     method = "wilcox.test",
                     hide.ns =F)+theme_bw()+theme(axis.text.x = element_text(hjust = 1,vjust = 1,angle = 45))+scale_fill_manual(values = c("#FB8072","#80B1D3"))

###Fig9###
library(survival)
library(survminer)
library(tidyverse)
load("../0810/model.Rdata")
phe <- read.csv("../820_enrichment/clinical_infor.csv",row.names = 1)
phe <- phe[-which(phe$age_at_index=="'--"),]
phe <- phe[phe$ajcc_pathologic_stage!="not report",]
colnames(phe)
phe <- phe[,c("ajcc_pathologic_stage","case_submitter_id","age_at_index","gender")]
colnames(phe) <- c('stage','sample','age','gender')

riskScore_cli$sample <- substr(riskScore_cli$sample,1,12)
expr <- merge(phe,riskScore_cli,"sample")
expr <- column_to_rownames(expr,var = "sample")
expr <- expr[,1:6]
a <- expr

a <- a[c(5,6,1,2,3,4)]
a$stage <- gsub(pattern = "I",replacement = "0",x = a$stage)
a$stage <- gsub(pattern = "00",replacement = "0",x = a$stage)
a$stage <- gsub(pattern = "00",replacement = "1",x = a$stage)
a$stage <- gsub(pattern = "0V",replacement = "1",x = a$stage)
table(a$stage)

a$stage <- as.factor(a$stage)
a$gender <- as.factor(a$gender)
a$age <- as.numeric(a$age)
a$riskScore <- as.numeric(a$riskScore)
#####
pFilter=10 #
outResult=data.frame() #
sigGenes=c("OS","OS.time") #
for(i in colnames(a[,3:ncol(a)])){ 
  tdcox <- coxph(Surv(OS.time, OS) ~ a[,i], data = a)#
  tdcoxSummary = summary(tdcox) #summary
  pvalue=tdcoxSummary$coefficients[,"Pr(>|z|)"] #
  if(pvalue<pFilter){ 
    sigGenes=c(sigGenes,i)
    outResult=rbind(outResult,
                    cbind(id=i,
                          HR=tdcoxSummary$conf.int[,"exp(coef)"],#
                          L95CI=tdcoxSummary$conf.int[,"lower .95"],#
                          H95CI=tdcoxSummary$conf.int[,"upper .95"],#
                          pvalue=tdcoxSummary$coefficients[,"Pr(>|z|)"])#
    )
  }
}

#####
write.table(outResult,file="../0810/forest/unicox.txt",sep="\t",row.names=F,quote=F)
tdmcs <- read.table("../0810/forest/unicox.txt",header=T,sep="\t",row.names=1)
gene <- rownames(tdmcs)
hr <- sprintf("%.3f",tdmcs$"HR")
hrLow  <- sprintf("%.3f",tdmcs$"L95CI")
hrHigh <- sprintf("%.3f",tdmcs$"H95CI")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pValue <- ifelse(tdmcs$pvalue<0.001, "<0.001", sprintf("%.3f", tdmcs$pvalue))

pdf(file="../0810/forest/UniCoxSurForestPlot.pdf", width = 6,height = 3)
n <- nrow(tdmcs)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(2.5,2))

xlim = c(0,2.5)
par(mar=c(4,2.5,2,1))
plot(0,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.2-0.5*0.2,n:1,pValue,adj=1,cex=text.cex);text(1.2-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(2.5,n:1,Hazard.ratio,adj=1,cex=text.cex);text(2.5,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)

par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(0,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'green')
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
axis(1)

dev.off()

############
tdmultiCox=coxph(Surv(OS.time, OS) ~ ., data =a)
tdmultiCoxSum=summary(tdmultiCox)
outResult=data.frame()
outResult=cbind(
  HR=tdmultiCoxSum$conf.int[,"exp(coef)"],
  L95CI=tdmultiCoxSum$conf.int[,"lower .95"],
  H95CIH=tdmultiCoxSum$conf.int[,"upper .95"],
  pvalue=tdmultiCoxSum$coefficients[,"Pr(>|z|)"])
outResult=cbind(id=row.names(outResult),outResult)
write.table(outResult,file="../0810/forest/multiCoxClinical.txt",sep="\t",row.names=F,quote=F)

####
tdmcs <- read.table("../0810/forest/multiCoxClinical.txt",header=T,sep="\t",row.names=1)
gene <- rownames(tdmcs)
hr <- sprintf("%.3f",tdmcs$"HR")
hrLow  <- sprintf("%.3f",tdmcs$"L95CI")
hrHigh <- sprintf("%.3f",tdmcs$"H95CI")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pValue <- ifelse(tdmcs$pvalue<0.001, "<0.001", sprintf("%.3f", tdmcs$pvalue))

pdf(file="../0810/forest/mutiCoxSurForestPlot.pdf", width = 6,height = 3)
n <- nrow(tdmcs)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(2.5,2))

xlim = c(0,2.5)
par(mar=c(4,2.5,2,1))
plot(0,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.2-0.5*0.2,n:1,pValue,adj=1,cex=text.cex);text(1.2-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(2.5,n:1,Hazard.ratio,adj=1,cex=text.cex);text(2.5,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)

par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(0,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'green')
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
axis(1)
dev.off()


######
library(survival)
library(survminer)
library(rms)

dd=datadist(a)
options(datadist="dd")

f2<-psm(Surv(OS.time,OS)~stage+age+gender+riskScore,data = a,dist='lognormal')
med<-Quantile(f2)#
surv<-Survival(f2)#

##
nom<-nomogram(f2,fun=function(x)med(lp=x),
              funlabel="Median Survival Time")
plot(nom)

##
nom<-nomogram(f2,lp = F,fun=list(function(x)surv(12,x),
                                 function(x)surv(24,x),
                                 function(x)surv(36,x)),
              funlabel=c("1-year Survival Probability",
                         "2-year Survival Probability",
                         '3-year Survival Probability'))
plot(nom,xfrac=.6)



####

rcorrcens(Surv(OS.time,OS)~predict(f2),data=a)

f2 <- psm(Surv(OS.time,OS) ~ stage+age+gender+riskScore,data = a,x=T, y=T, dist='lognormal') 
## 

cal1 <- calibrate(f2, 
                  cmethod='KM', 
                  method="boot", 
                  u=24, # 
                  m=66, #
                  B=1000) 
cal2 <- calibrate(f2, 
                  cmethod='KM', 
                  method="boot", 
                  u=48, 
                  m=66, 
                  B=1000)

cal3 <- calibrate(f2, 
                  cmethod='KM', 
                  method="boot", 
                  u=36, 
                  m=66, 
                  B=1000)



plot(cal1,lwd=1,lty=1,
     conf.int=T,
     errbar.col="black",
     col="red", 
     xlim=c(0,1),ylim=c(0,1),
     xlab="Nomogram-Predicted Probability DFS",
     ylab="Actual DFS (proportion)",
     subtitles = F,
     lines(cal1[,c('mean.predicted',"KM")],
           type = 'b', lwd = 3, col = c("red"), pch = 16))

plot(cal2,lwd=1,lty=1,
     conf.int=T,# 
     errbar.col="black",
     col="blue", 
     xlim=c(0,1),ylim=c(0,1),
     xlab="Nomogram-Predicted Probability DFS",
     ylab="Actual DFS (proportion)",
     subtitles = F,
     lines(cal1[,c('mean.predicted',"KM")],
           type = 'b', lwd = 3, col = c("blue"), pch = 16))

plot(cal3,lwd=1,lty=1,
     conf.int=T,
     errbar.col="black",
     col="purple", 
     xlim=c(0,1),ylim=c(0,1),
     xlab="Nomogram-Predicted Probability DFS",
     ylab="Actual DFS (proportion)",
     subtitles = F)

###figS3#####
########
load("../data/455sample.Rdata")
dat1 <- dat[,group_list=="Tumor"]
load("../0810/model.Rdata")

dat1 <- dat1[,colnames(dat1)%in%riskScore_cli1$sample]

riskScore_cli1<- riskScore_cli1 %>% 
  mutate(sample=fct_relevel(sample,colnames(dat1))) %>% 
  arrange(sample)
######
library(limma)
library(ggplot2)
group = riskScore_cli1$riskScore2
design <- model.matrix(~0+factor(group))
colnames(design)=levels(factor(group))   
rownames(design)=colnames(dat1)
contrast.matrix <- makeContrasts("High-Low",levels = design) 

fit <- lmFit(dat1,design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2) 
tempOutput = topTable(fit2, coef=1, n=Inf)
DEG = na.omit(tempOutput)
DEG_LCI <- DEG
save(DEG_LCI,file = "../820_enrichment/DEG_LCI.Rdata")


###KEGG###
library(clusterProfiler)
library(dplyr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
library(RColorBrewer)

#####KEGG###
load("../820_enrichment/DEG_UCSC.Rdata")

up<-  mapIds(org.Hs.eg.db,rownames(DEG_UCSC[DEG_UCSC$logFC>1,]),'ENTREZID', 'SYMBOL')    
up_KEGG<- enrichKEGG(gene = up,
                     organism = "hsa",
                     pvalueCutoff =0.05,
                     qvalueCutoff = 0.2)
write.table(up_KEGG,file="../820_enrichment/UCSC_KEGG.txt",sep="\t",quote=F,row.names = F)

up_KEGG2= setReadable(up_KEGG, 
                      OrgDb = "org.Hs.eg.db", 
                      keyType = "ENTREZID")


######
library(tidyverse)
ek.rt = read.table("../820_enrichment/UCSC_KEGG.txt",header=TRUE,sep="\t",quote = "")  
ek.rt <- separate(data=ek.rt, col=GeneRatio, into = c("GR1", "GR2"), sep = "/") 
ek.rt <- separate(data=ek.rt, col=BgRatio, into = c("BR1", "BR2"), sep = "/") 
ek.rt <- mutate(ek.rt, enrichment_factor = (as.numeric(GR1)/as.numeric(GR2))/(as.numeric(BR1)/as.numeric(BR2))) 

ek.rt10 <- ek.rt %>% filter(row_number() >= 1,row_number() <= 10)


p <- ggplot(ek.rt,aes(enrichment_factor, fct_reorder(factor(Description), enrichment_factor))) + 
  geom_point(aes(size=Count,color=-1*log10(pvalue))) +
  scale_color_gradient(low="green",high ="red") + 
  labs(color=expression(-log[10](p_value)),size="Count",
       x="Enrichment Factor",y="KEGG term",title="KEGG enrichment") + 
  theme_bw()
p

####GO####
up_GO= enrichGO(gene =up,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "ALL",
                pvalueCutoff = 0.05)

up_GO2= setReadable(up_GO, 
                    OrgDb = "org.Hs.eg.db",
                    keyType = "ENTREZID")
GO_res <- up_GO2@result

dotplot(up_GO, split="ONTOLOGY",font.size =10,showCategory =5)+ 
  facet_grid(ONTOLOGY~.,scale="free")+scale_color_gradient(low="green",high ="red")


######LCI##################
##LCI#######
rm(list = ls())
options(stringsAsFactors = F)

getOption('timeout')
options(timeout=10000)

library(AnnoProbe)
library(GEOquery) 
gset <- geoChina("GSE14520")
gset
gset[[1]]
a=gset[[1]] 
dat=exprs(a) 
dim(dat)
dat[1:4,1:4] 
boxplot(dat[,1:4],las=2)  

#dat=log2(dat)
boxplot(dat[,1:4],las=2)  
library(limma)
dat=normalizeBetweenArrays(dat)
boxplot(dat[,1:4],las=2)  
pd=pData(a) 



library(stringr) 
# pd=pd[1:43,]
# dat=dat[,1:43]
phe=pd 
table(phe$characteristics_ch1)
group_list=ifelse(grepl('Non-Tumor',phe$characteristics_ch1),'Normal','Tumor')
table(group_list)

dat[1:4,1:4]  
dim(dat)
a  
ids <- select(hthgu133a.db, keys=row.names(dat), columns = c("SYMBOL","GENETYPE"))
head(ids)
ids=ids[ids$SYMBOL != '',]
dat=dat[rownames(dat) %in% ids$PROBEID,]
ids=ids[match(rownames(dat),ids$PROBEID),]
head(ids) 
ids <- ids[,1:2]
colnames(ids)=c('probe_id','symbol')  
ids$probe_id=as.character(ids$probe_id)
rownames(dat)=ids$probe_id
dat[1:4,1:4] 

ids=ids[ids$probe_id %in%  rownames(dat),]
dat[1:4,1:4]   
dat=dat[ids$probe_id,] 

ids$median=apply(dat,1,median)
ids=ids[order(ids$symbol,ids$median,decreasing = T),]
ids=ids[!duplicated(ids$symbol),]
dat=dat[ids$probe_id,] 
rownames(dat)=ids$symbol
dat[1:4,1:4]

dat['ALB',]
dat['CXCR5',]
save(dat,group_list,phe,file = 'Download data/Fudan/LCI/455sample.Rdata')

#####Fig5D-F###
library(CIBERSORT)
load("../data/455sample.Rdata")
dat <- as.matrix(dat)
b <- cibersort(sig_matrix = LM22,mixture_file = dat)


result2 <- as.data.frame(b)
result2 <- result2[,1:22]

xx <- data.frame(ID=colnames(dat),tissue=group_list)
xx <- xx[xx$tissue=="Tumor",]

result2$ID <- rownames(result2)
result2 <- merge(xx,result2,"ID")
LCI_clin <- read.csv("../data/LCI_all_tumor_clinical.csv")
colnames(LCI_clin)[3] <- "ID"
LCI_clin[1:2] <- NULL

sur2 <- merge(LCI_clin,result2,"ID")
sur2$naive <- sur2$`B cells naive`
sur2$memory <- sur2$`B cells memory`
sur2$Plasmacells <- sur2$`Plasma cells`


res.cut2 <- surv_cutpoint(sur2, time = "Survival.months", event = "Survival.status",
                          variables = c("naive","memory","Plasmacells"))                                  
summary(res.cut2)
res.cat2 <- surv_categorize(res.cut2)
head(res.cat2)

fit2 <- survfit(Surv(Survival.months, Survival.status) ~Plasmacells, data = res.cat2)
ggsurvplot(fit2, data = res.cat2,palette= c("red","blue"),risk.table =T, conf.int = F,pval = T,
           legend.labs=NULL)+xlab(label = "Months")

fit3 <- survfit(Surv(Survival.months, Survival.status) ~naive, data = res.cat2)
ggsurvplot(fit3, data = res.cat2,palette= c("red","blue"),risk.table =T, conf.int = F,pval = T,
           legend.labs=NULL)+xlab(label = "Months")

fit4 <- survfit(Surv(Survival.months, Survival.status) ~memory, data = res.cat2)
ggsurvplot(fit4, data = res.cat2,palette= c("red","blue"),risk.table =T, conf.int = F,pval = T,
           legend.labs=NULL)+xlab(label = "Months")


ggsurvplot(fit2, data = res.cat2,palette= c("red","blue"),risk.table = F, conf.int = F,pval = T,
           legend.labs=NULL)+xlab(label = "Months")  #3.36,4.42

###Fig6D-E##########
library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggalluvial)
library(RColorBrewer)
load("../data/integrated.final.Rdata")
markers <- FindAllMarkers(integrated,min.pct = 0.2,only.pos = T)
markers1 <-markers[markers$p_val_adj<0.05,]
plasma_markers <- markers1[markers1$cluster=="Plasma cells",]
save(markers1,plasma_markers,file = "../0810/deg_all.Rdata")

#####UNICOX####
library(KEGGREST)
library(GSVA)
library(survival)
library(survminer)
library(limma)
expr <- readRDS("../data/expr_356samples.Rdata")
expr[1:5,1:5]
clinical <- read.delim("../data/UCSC-LIHC.survival.tsv")
clinical$OS.time <- clinical$OS.time/30
colnames(clinical)[1] <- "id"

expr <- expr[,substr(colnames(expr),14,15)<10]
expr <- expr[,!duplicated(substr(colnames(expr),1,12))]
load("../0810/deg_all.Rdata")
expr<- expr[rownames(expr)%in%plasma_markers$gene,]

expr <- as.data.frame(t(expr))

rownames(expr)<- gsub(".","-",rownames(expr),fixed = T)
expr$id <- rownames(expr)

td <- merge(clinical,expr,by = "id")
rownames(td) <- td$id
td[1] <- NULL

pFilter=0.05 
outResult=data.frame() 
sigGenes=c("OS","OS.time") 

for(i in colnames(td[,4:ncol(td)])){ 
  tdcox <- coxph(Surv(OS.time, OS) ~ td[,i], data = td)
  tdcoxSummary = summary(tdcox) 
  pvalue=tdcoxSummary$coefficients[,"Pr(>|z|)"] 
  if(pvalue<pFilter){ 
    sigGenes=c(sigGenes,i)
    outResult=rbind(outResult,
                    cbind(id=i,
                          HR=tdcoxSummary$conf.int[,"exp(coef)"],
                          L95CI=tdcoxSummary$conf.int[,"lower .95"],
                          H95CI=tdcoxSummary$conf.int[,"upper .95"],
                          pvalue=tdcoxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
}

write.table(outResult,file="../0810/UniCoxSurvival.txt",sep="\t",row.names=F,quote=F)

UniCoxSurSigGeneExp=td[,sigGenes] 
UniCoxSurSigGeneExp=cbind(id=row.names(UniCoxSurSigGeneExp),UniCoxSurSigGeneExp)
write.table(UniCoxSurSigGeneExp,file="../0810/UniCoxSurSigGeneExp.txt",sep="\t",row.names=F,quote=F)

#####LASSO####
library(glmnet)
library(caret)
library(survival)
library(survminer)

outResult <- read.table("../0810/UniCoxSurvival.txt",header = T)
UniCoxSurSigGeneExp <- read.table("../0810/UniCoxSurSigGeneExp.txt",header = T)

####59risk
x<- UniCoxSurSigGeneExp[,4:64]
x <- x[,outResult[outResult$HR>1,]$id]
rownames(x) <- UniCoxSurSigGeneExp$id
x <- as.matrix(x)

####
y <- as.matrix(UniCoxSurSigGeneExp[,2])
rownames(y) <- UniCoxSurSigGeneExp$id
identical(rownames(x),rownames(y))
###
library(glmnet)
model_lasso <- glmnet(x,y,alpha = 1,family = "binomial")
print(model_lasso)
plot(model_lasso)
plot(model_lasso, xvar = "lambda", label = TRUE)

set.seed(123)
cv_fit <- cv.glmnet(x=x, y=y,alpha = 1)
plot(cv_fit)

####
model_lasso_min <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.min)
model_lasso_1se <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.1se)


head(model_lasso_min$beta)
choose_gene_min=rownames(model_lasso_min$beta)[as.numeric(model_lasso_min$beta)!=0]
choose_gene_1se=rownames(model_lasso_1se$beta)[as.numeric(model_lasso_1se$beta)!=0]
length(choose_gene_min)
length(choose_gene_1se)



lasso.prob <- predict(cv_fit, 
                      newx=x , 
                      s=c(cv_fit$lambda.min,cv_fit$lambda.1se) )
re=cbind(y ,lasso.prob) 
head(re)

##
re=as.data.frame(re)
colnames(re)=c('event','prob_min','prob_1se')
#
re$event=as.factor(re$event) 
library(ggpubr) 
p1 = ggboxplot(re, x = "event", y = "prob_min",
               color = "event", palette = "jco",
               add = "jitter")+ stat_compare_means()
p2 = ggboxplot(re, x = "event", y = "prob_1se",
               color = "event", palette = "jco",
               add = "jitter")+ stat_compare_means()
library(patchwork)
p1+p2

##

library(ROCR)
library(caret)

pred_min <- prediction(re[,2], re[,1])
auc_min = performance(pred_min,"auc")@y.values[[1]]
perf_min <- performance(pred_min,"tpr","fpr")
plot(perf_min,colorize=FALSE, col="blue") 
lines(c(0,1),c(0,1),col = "gray", lty = 4 )
text(0.8,0.2, labels = paste0("AUC = ",round(auc_min,2)))

#1se
pred_1se <- prediction(re[,3], re[,1])
auc_1se = performance(pred_1se,"auc")@y.values[[1]]
perf_1se <- performance(pred_1se,"tpr","fpr")
plot(perf_1se,colorize=FALSE, col="red") 
lines(c(0,1),c(0,1),col = "gray", lty = 4 )
text(0.8,0.2, labels = paste0("AUC = ",round(auc_1se,2)))


plot(perf_min,colorize=FALSE, col="blue") 
plot(perf_1se,colorize=FALSE, col="red",add = T) 
lines(c(0,1),c(0,1),col = "gray", lty = 4 )
text(0.8,0.3, labels = paste0("AUC = ",round(auc_min,3)),col = "blue")
text(0.8,0.2, labels = paste0("AUC = ",round(auc_1se,3)),col = "red")


library(ROCR)
library(caret)

pred_min <- prediction(re[,2], re[,1])
auc_min = performance(pred_min,"auc")@y.values[[1]]
perf_min <- performance(pred_min,"tpr","fpr")

#1se
pred_1se <- prediction(re[,3], re[,1])
auc_1se = performance(pred_1se,"auc")@y.values[[1]]
perf_1se <- performance(pred_1se,"tpr","fpr")

tpr_min = performance(pred_min,"tpr")@y.values[[1]]
tpr_1se = performance(pred_1se,"tpr")@y.values[[1]]
dat = data.frame(tpr_min = perf_min@y.values[[1]],
                 fpr_min = perf_min@x.values[[1]],
                 tpr_1se = perf_1se@y.values[[1]],
                 fpr_1se = perf_1se@x.values[[1]])

ggplot()+ 
  geom_line(data = dat,aes(x = fpr_min, y = tpr_min),color = "blue") + 
  geom_line(data = dat,aes(x = fpr_1se, y = tpr_1se),color = "red")+
  geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey")+
  theme_bw()+
  annotate("text",x = .75, y = .25,
           label = paste("AUC of min = ",round(auc_min,2)),color = "blue")+
  annotate("text",x = .75, y = .15,label = paste("AUC of 1se = ",round(auc_1se,2)),color = "red")+
  scale_x_continuous(name  = "fpr")+
  scale_y_continuous(name = "tpr")


sig_gene_mult_cox <- choose_gene_min
DEG_met_expr.lasso_cox <- UniCoxSurSigGeneExp[,c("OS","OS.time",sig_gene_mult_cox)]
rownames(DEG_met_expr.lasso_cox) <- UniCoxSurSigGeneExp$id

multiCox <- coxph(Surv(OS.time, OS) ~ ., data =  DEG_met_expr.lasso_cox)
summary(multiCox)
###
load("../data/455sample.Rdata")
dat1 <- dat[,group_list=="Tumor"]
LCI_clin <- read.csv("../data/LCI_all_tumor_clinical.csv")
colnames(LCI_clin)[3] <- "ID"
LCI_clin[1:2] <- NULL
LCI_clin <- LCI_clin[,c("ID","Survival.status","Survival.months")]

dat1 <- dat1[rownames(dat1)%in%names(multiCox$coefficients),]
dat1 <- as.data.frame(t(dat1))
dat1$ID <- rownames(dat1)

dat2 <- merge(dat1,LCI_clin,by = "ID")

riskScore1=apply(dat2[,names(multiCox$coefficients)], 1,function(k){
  sum(multiCox$coefficients * k)
}) 

riskScore1<-as.data.frame(riskScore1)
riskScore1$sample <- dat2$ID
head(riskScore1,2)

riskScore_cli1 <- riskScore1
riskScore_cli1$OS <- dat2$Survival.status
riskScore_cli1$OS.time <- dat2$Survival.months
riskScore_cli1$riskScore2 <- ifelse(riskScore_cli1$riskScore1 > median(riskScore_cli1$riskScore1),
                                    "High","Low")

res.cut <- surv_cutpoint(riskScore_cli1, time = "OS.time", event = "OS",variables = c("riskScore1"))  
summary(res.cut)
res.cat <- surv_categorize(res.cut)
head(res.cat)
table(res.cat$riskScore1)
fit <- survfit(Surv(OS.time, OS) ~riskScore1, data = res.cat)

ggsurvplot(fit, data = res.cat,
           pval = T,
           risk.table = T,
           #surv.median.line = "hv", 
           palette=c("red", "blue"), 
           #legend.labs=c("High risk","Low risk"), 
           legend.title="RiskScore",
           title="Overall survival", #
           ylab="Cumulative survival (percentage)",xlab = " Time (Months)", 
           censor.shape = 124,censor.size = 2,conf.int = FALSE)
res.cat$sample <- riskScore_cli1$sample
save(riskScore_cli,riskScore_cli1,DEG_met_expr.lasso_cox,dat2,file = "../0810/model.Rdata")
save(res.cat,file = "../0810/Cutpoint/cutpoint.Rdata")

library(timeROC)
with(riskScore_cli1,
     ROC_riskscore <<- timeROC(T = OS.time,
                               delta = OS,
                               marker = riskScore,
                               cause = 1,
                               weighting = "marginal",
                               times = c(12,24,36),
                               ROC = TRUE,
                               iid = TRUE))

plot(ROC_riskscore, time = 12, col = "red", add = F,title = "")
plot(ROC_riskscore, time = 24, col = "blue", add = T)
plot(ROC_riskscore, time = 36, col = "purple", add = T)
legend("bottomright",c("1-Year","3-Year","5-Year"),col=c("red","blue","purple"),lty=1,lwd=2)
text(0.5,0.2,paste("1-Year AUC = ",round(ROC_riskscore$AUC[1],3)))
text(0.5,0.15,paste("2-Year AUC = ",round(ROC_riskscore$AUC[2],3)))
text(0.5,0.1,paste("3-Year AUC = ",round(ROC_riskscore$AUC[3],3)))

######
###GSE149614####
library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggsci)
library(patchwork)

load("../01_data/PT.Rdata")
scRNA <- CreateSeuratObject(PT)

###
Sample <- substr(colnames(PT),1,6)
patient<- substr(Sample,1,5)
tissue <- ifelse(substr(Sample,6,6)=="N","Normal","Tumor")
table(patient)

group <- data.frame(Samples=Sample,Tissue=tissue,Patient=patient)
rownames(group) <- colnames(scRNA)

scRNA <- AddMetaData(scRNA,metadata =group)


table(scRNA@meta.data$Samples)
table(scRNA@meta.data$Tissue)

saveRDS(scRNA,"../01_data/scRNA_meta.Rdata")

#####
scRNA[["percent.mt"]]= PercentageFeatureSet(scRNA,pattern = "^MT-")

feats <- c("nFeature_RNA", "nCount_RNA","percent.mt")
p1=VlnPlot(scRNA, group.by = "Samples", features = feats, pt.size =0.01)+NoLegend()
ggsave(p1,filename = "../03_result/nFeature_RNA_QC_before.pdf",width = 16.86,height = 12.61)


####
minGene=200
maxGene=8000
pctMT=10
minUMI=200

scRNA <- subset(scRNA, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene & percent.mt < pctMT&nCount_RNA>minUMI)
p1=VlnPlot(scRNA, group.by = "Samples", features = feats, pt.size =0.01)+NoLegend()

ggsave(p1,filename = "../03_result/nFeature_RNA_QC_after.pdf",width = 16.86,height = 12.61)



###
sce <- scRNA
sce <- NormalizeData(sce, 
                     normalization.method = "LogNormalize",
                     scale.factor = 1e4) 
sce <- FindVariableFeatures(sce,nfeatures = 3000)
sce <- ScaleData(sce)
sce <- RunPCA(sce, features = VariableFeatures(object = sce))
ElbowPlot(sce,ndims = 50)


#####
library(harmony)
seuratObj <- RunHarmony(sce, "orig.ident",max.iter.harmony = 20)
names(seuratObj@reductions)
ElbowPlot(seuratObj,ndims = 50,reduction = "harmony")

ggsave(filename = "../03_result/harmony.pdf",width = 5.86,height = 4.51)

ndims <- 1:40

seuratObj <- RunUMAP(seuratObj,  dims = ndims, 
                     reduction = "harmony")

seuratObj <- RunTSNE(seuratObj,  dims = ndims, 
                     reduction = "harmony")

DimPlot(seuratObj,reduction = "umap",label=F,pt.size = 0.5) 

sce=seuratObj
sce <- FindNeighbors(sce, reduction = "harmony",
                     dims = ndims) 
sce <- FindClusters(sce,resolution = 0.8)
DimPlot(sce,reduction = 'tsne',label = T,pt.size = 0.5)

###feature plot###
endothelial <- FeaturePlot(sce,features = "PECAM1",reduction = "tsne",cols =c( 'lightgrey', '#EC7B1B', 'red'),pt.size = 0.5)
fibroblast <- FeaturePlot(sce,features = "RGS5",reduction = "tsne",cols =c( 'lightgrey', '#EC7B1B', 'red'),pt.size = 0.5)
B <- FeaturePlot(sce,features = "MS4A1",reduction = "tsne",cols =c( 'lightgrey', '#EC7B1B', 'red'),pt.size = 0.5)
plasma <- FeaturePlot(sce,features = "MZB1",reduction = "tsne",cols =c( 'lightgrey', '#EC7B1B', 'red'),pt.size = 0.5)
Tcell <- FeaturePlot(sce,features = "CD3D",reduction = "tsne",cols =c( 'lightgrey', '#EC7B1B', 'red'),pt.size = 0.5)  
NK<- FeaturePlot(sce,features = "KLRB1",reduction = "tsne",cols =c( 'lightgrey', '#EC7B1B', 'red'),pt.size = 0.5)
myeloid<- FeaturePlot(sce,features = "CD68",reduction = "tsne",cols =c( 'lightgrey', '#EC7B1B', 'red'),pt.size = 0.5)
DC1 <- FeaturePlot(sce,features = "CD1C",reduction = "tsne",cols =c( 'lightgrey', '#EC7B1B', 'red'),pt.size = 0.5)
DC2 <- FeaturePlot(sce,features = "IDO1",reduction = "tsne",cols =c( 'lightgrey', '#EC7B1B', 'red'),pt.size = 0.5)
hypatocyte <- FeaturePlot(sce,features = "SERPINA1",reduction = "tsne",cols =c( 'lightgrey', '#EC7B1B', 'red'),pt.size = 0.5)

library(ggsci)
library(ggplot2)
library(patchwork)
library(MySeuratWrappers)
celltype=data.frame(ClusterID=0:27,
                    celltype=0:27)   
celltype[celltype$ClusterID %in% c(0,3,4,10,19,26,14),2]='T/NK' 
celltype[celltype$ClusterID %in% c(12,8,2,17,20,7,21,16),2]='Hepatocyte'   

celltype[celltype$ClusterID %in% c(25,6,13),2]='Endothelial' 

celltype[celltype$ClusterID %in% c(1,5,15,27,22,23,24),2]='Myeloid' 
celltype[celltype$ClusterID %in% c(9),2]="Fibroblast" 

celltype[celltype$ClusterID %in% c(11),2]='B cell' 
celltype[celltype$ClusterID %in% c(18),2]='Plasma'



for(i in 1:nrow(celltype)){
  sce@meta.data[which(sce@meta.data$RNA_snn_res.0.8 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}

table(sce@meta.data$celltype)

genes_to_check = c('CD79A', 'MS4A1' ,  #B
                   'PECAM1', 'VWF', ## endo
                   "RGS5","COL1A1","PDGFRB",#hsc
                   "ALB","SERPINA1","APOH",#hcc
                   'CD68', 'CD163', 'CD14','C1QA',  'C1QB', ##macro
                   'IGHG1', 'MZB1',    #plasma cells#
                   'PTPRC', 'CD3D', 'CD3E') #T 

genes_to_check=unique(genes_to_check) 
p <- DotPlot(sce, features = unique(genes_to_check),
             assay='RNA' ,group.by = 'celltype' )  +scale_color_gradientn(values = seq(0,1,0.1),colours = c("white","grey","blue","red"))+
  labs(y = "Cell type", x = "Features")+ coord_flip()

ggsave(p,filename = "../03_result/dotplot.pdf",width = 10.86,height = 7.51)

Idents(sce)=sce$celltype
saveRDS(sce,"../01_data/sce_celltype.Rdata")

##Fig2C####
Bcell <- subset(sce,idents = c("B cell","Plasma"))
library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggpubr)
Bcell <-readRDS("../01_data/Bcell.rdata")
scRNA <- CreateSeuratObject(Bcell@assays$RNA@counts)

Sample <- substr(colnames(scRNA),1,6)
patient<- substr(Sample,1,5)
tissue <- ifelse(substr(Sample,6,6)=="N","Normal","Tumor")
group <- data.frame(Samples=Sample,Tissue=tissue,Patient=patient)
rownames(group) <- colnames(scRNA)
scRNA <- AddMetaData(scRNA,metadata =group)


table(scRNA@meta.data$Samples)
table(scRNA@meta.data$Tissue)

###
sce <- scRNA
sce <- NormalizeData(sce, 
                     normalization.method = "LogNormalize",
                     scale.factor = 1e4) 
sce <- FindVariableFeatures(sce,nfeatures = 3000)
sce <- ScaleData(sce)
sce <- RunPCA(sce, features = VariableFeatures(object = sce))
ElbowPlot(sce,ndims = 50)

#####
library(harmony)
seuratObj <- RunHarmony(sce, "orig.ident",max.iter.harmony = 20)
names(seuratObj@reductions)
ElbowPlot(seuratObj,ndims = 50,reduction = "harmony")

ndims <- 1:30

seuratObj <- RunUMAP(seuratObj,  dims = ndims, 
                     reduction = "harmony")

seuratObj <- RunTSNE(seuratObj,  dims = ndims, 
                     reduction = "harmony")

DimPlot(seuratObj,reduction = "umap",label=F,pt.size = 0.5) 

sce=seuratObj
sce <- FindNeighbors(sce, reduction = "harmony",
                     dims = ndims) 
sce <- FindClusters(sce,resolution = 0.5)

FeaturePlot(sce,features = "MS4A1",pt.size = 1.5,reduction = "tsne")
ggsave(filename = "../03_result/B cell/MS4A1.pdf",width = 6.01,height = 4.8)

FeaturePlot(sce,features = "IGHG1",pt.size = 1.5,reduction = "tsne")
ggsave(filename = "../03_result/B cell/IGHG1.pdf",width = 6.01,height = 4.8)

celltype=data.frame(ClusterID=0:5,
                    celltype=0:5)   

celltype[celltype$ClusterID %in% c(0,2,4),2]='B cell'   

celltype[celltype$ClusterID %in% c(1,3,5),2]='Plasma' 

for(i in 1:nrow(celltype)){
  sce@meta.data[which(sce@meta.data$RNA_snn_res.0.5 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}

table(sce@meta.data$celltype)

Idents(sce)=sce$celltype

####
Bcell <- sce
table(Idents(Bcell), Bcell$Samples)
Cellratio <- prop.table(table(Idents(Bcell), Bcell$Samples), margin = 2)
Cellratio <- as.data.frame(Cellratio)
Cellratio$Tissue <- rep(c("Normal","Normal","Tumor","Tumor"),4)
B <- Cellratio[Cellratio$Var1=="B cells",]
plasma <- Cellratio[Cellratio$Var1=="Plasma cells",]

mycolor <- c("#6C9D88","#6980AC")

p1 <- ggplot(B, aes(x = Tissue, y = Freq)) +
  geom_boxplot(aes(fill = Tissue), show.legend = F, width = 0.3) +  #
  scale_fill_manual(values = c("#6C9D88","#6980AC")) + 
  geom_point(size = 1.5) +#
  theme(panel.grid = element_blank(), 
        axis.line = element_line(colour = 'black', size = 0.5), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 15, hjust = 0.5), 
        plot.subtitle = element_text(size = 15, hjust = 0.5), 
        axis.text = element_text(size = 15, color = 'black'),
        axis.title = element_text(size = 15, color = 'black')) +
  labs(x = '', y = 'Percentages', title = 'Tumor vs Normal')+
  stat_compare_means(method = "wilcox.test",paired = F, comparisons=list(c("Tumor", "Normal")))+ggtitle(label = "B cells")

ggsave(p1,filename = "../PDF/Fig2/Bcell.per.pdf",width = 4.04,height = 4.34)


p2 <- ggplot(plasma, aes(x = Tissue, y = Freq)) +
  geom_boxplot(aes(fill = Tissue), show.legend = F, width = 0.3) +  #
  scale_fill_manual(values = c("#6C9D88","#6980AC")) + 
  geom_point(size = 1.5) +#
  theme(panel.grid = element_blank(), 
        axis.line = element_line(colour = 'black', size = 0.5), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 15, hjust = 0.5), 
        plot.subtitle = element_text(size = 15, hjust = 0.5), 
        axis.text = element_text(size = 15, color = 'black'),
        axis.title = element_text(size = 15, color = 'black')) +
  labs(x = '', y = 'Percentages', title = 'Tumor vs Normal')+
  stat_compare_means(method = "wilcox.test",paired = F, comparisons=list(c("Tumor", "Normal")))+ggtitle(label = "B cells")

ggsave(p2,filename = "../PDF/Fig2/Plasmacell.per.pdf",width = 4.04,height = 4.34)







































