screen -r  220142

export LSF_DOCKER_VOLUMES='/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/:/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/'
bsub -Is -G compute-gjrandolph -q general-interactive -R 'rusage[mem=128GB]' -a 'docker(jichanghan/scrna_seurat4_monocle2_wgcna_0213:latest)' /bin/bash

library(Seurat)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)   
library(sctransform)
library(harmony)
library(metap)
library(multtest)
library(biomaRt) 
library(monocle3)

rm(list=ls()) 


myinf1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/raw/Nan_PW_Normal"
PW_C.data <- Read10X(data.dir=myinf1)
PW_C <- CreateSeuratObject(counts = PW_C.data, project = "PW_C",min.cells=3,min.features=200)
PW_C<-RenameCells(object=PW_C,add.cell.id=unique(PW_C@meta.data$orig.ident))
PW_C<-RenameCells(object=PW_C,new.names=gsub('-1','',Cells(PW_C)))

myinf4 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/raw/JH_PW_Kid1/"
Kid1.data <- Read10X(data.dir=myinf4)
Kid1 <- CreateSeuratObject(counts = Kid1.data, project = "Kid1",min.cells=3,min.features=200)
Kid1<-RenameCells(object=Kid1,add.cell.id=unique(Kid1@meta.data$orig.ident))
Kid1<-RenameCells(object=Kid1,new.names=gsub('-1','',Cells(Kid1)))


myinf5 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/raw/JH_PW_Kid2/"
Kid2.data <- Read10X(data.dir=myinf5)
Kid2 <- CreateSeuratObject(counts = Kid2.data, project = "Kid2",min.cells=3,min.features=200)
Kid2<-RenameCells(object=Kid2,add.cell.id=unique(Kid2@meta.data$orig.ident))
Kid2<-RenameCells(object=Kid2,new.names=gsub('-1','',Cells(Kid2)))


myinf6 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/raw/JH_PW_Kid4/"
Kid4.data <- Read10X(data.dir=myinf6)
Kid4 <- CreateSeuratObject(counts = Kid4.data, project = "Kid4",min.cells=3,min.features=200)
Kid4<-RenameCells(object=Kid4,add.cell.id=unique(Kid4@meta.data$orig.ident))
Kid4<-RenameCells(object=Kid4,new.names=gsub('-1','',Cells(Kid4)))


myinf7 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/raw/JH_PW_Kid6/"
Kid6.data <- Read10X(data.dir=myinf7)
Kid6 <- CreateSeuratObject(counts = Kid6.data, project = "Kid6",min.cells=3,min.features=200)
Kid6<-RenameCells(object=Kid6,add.cell.id=unique(Kid6@meta.data$orig.ident))
Kid6<-RenameCells(object=Kid6,new.names=gsub('-1','',Cells(Kid6)))


myinf8 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/raw/JH_PW_Kid7/"
Kid7.data <- Read10X(data.dir=myinf8)
Kid7 <- CreateSeuratObject(counts = Kid7.data, project = "Kid7",min.cells=3,min.features=200)
Kid7<-RenameCells(object=Kid7,add.cell.id=unique(Kid7@meta.data$orig.ident))
Kid7<-RenameCells(object=Kid7,new.names=gsub('-1','',Cells(Kid7)))


myinf9 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/raw/JH_PW_Kid8/"
Kid8.data <- Read10X(data.dir=myinf9)
Kid8 <- CreateSeuratObject(counts = Kid8.data, project = "Kid8",min.cells=3,min.features=200)
Kid8<-RenameCells(object=Kid8,add.cell.id=unique(Kid8@meta.data$orig.ident))
Kid8<-RenameCells(object=Kid8,new.names=gsub('-1','',Cells(Kid8)))


myinf10 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/raw/JH_PW_Kid9/"
Kid9.data <- Read10X(data.dir=myinf10)
Kid9 <- CreateSeuratObject(counts = Kid9.data, project = "Kid9",min.cells=3,min.features=200)
Kid9<-RenameCells(object=Kid9,add.cell.id=unique(Kid9@meta.data$orig.ident))
Kid9<-RenameCells(object=Kid9,new.names=gsub('-1','',Cells(Kid9)))


myinf11 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/raw/JH_PW_Kid10/"
Kid10.data <- Read10X(data.dir=myinf11)
Kid10 <- CreateSeuratObject(counts = Kid10.data, project = "Kid10",min.cells=3,min.features=200)
Kid10<-RenameCells(object=Kid10,add.cell.id=unique(Kid10@meta.data$orig.ident))
Kid10<-RenameCells(object=Kid10,new.names=gsub('-1','',Cells(Kid10)))

myinf12 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/raw/JH_PW_Kid11/"
Kid11.data <- Read10X(data.dir=myinf12)
Kid11 <- CreateSeuratObject(counts = Kid11.data, project = "Kid11",min.cells=3,min.features=200)
Kid11<-RenameCells(object=Kid11,add.cell.id=unique(Kid11@meta.data$orig.ident))
Kid11<-RenameCells(object=Kid11,new.names=gsub('-1','',Cells(Kid11)))


myinf13 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/raw/JH_PW_Adu_Q/"
PW_Q.data <- Read10X(data.dir=myinf13)
PW_Q <- CreateSeuratObject(counts = PW_Q.data, project = "PW_Q",min.cells=3,min.features=200)
PW_Q<-RenameCells(object=PW_Q,add.cell.id=unique(PW_Q@meta.data$orig.ident))
PW_Q<-RenameCells(object=PW_Q,new.names=gsub('-1','',Cells(PW_Q)))

myinf14 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/raw/JH_PW_Adu_S/"
PW_S.data <- Read10X(data.dir=myinf14)
PW_S <- CreateSeuratObject(counts = PW_S.data, project = "PW_S",min.cells=3,min.features=200)
PW_S<-RenameCells(object=PW_S,add.cell.id=unique(PW_S@meta.data$orig.ident))
PW_S<-RenameCells(object=PW_S,new.names=gsub('-1','',Cells(PW_S)))


myinf15 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/raw/JH_PW_Adu_U/"
PW_U.data <- Read10X(data.dir=myinf15)
PW_U <- CreateSeuratObject(counts = PW_U.data, project = "PW_U",min.cells=3,min.features=200)
PW_U<-RenameCells(object=PW_U,add.cell.id=unique(PW_U@meta.data$orig.ident))
PW_U<-RenameCells(object=PW_U,new.names=gsub('-1','',Cells(PW_U)))



myinf16 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/raw/JH_PW_Adu_W/"
PW_W.data <- Read10X(data.dir=myinf16)
PW_W <- CreateSeuratObject(counts = PW_W.data, project = "PW_W",min.cells=3,min.features=200)
PW_W<-RenameCells(object=PW_W,add.cell.id=unique(PW_W@meta.data$orig.ident))
PW_W<-RenameCells(object=PW_W,new.names=gsub('-1','',Cells(PW_W)))


myinf17 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/raw/JH_PW_Adu1/sample_feature_bc_matrix/"
PW_Adu1.data <- Read10X(data.dir=myinf17)
PW_Adu1 <- CreateSeuratObject(counts = PW_Adu1.data, project = "PW_Adu1",min.cells=3,min.features=200)
PW_Adu1<-RenameCells(object=PW_Adu1,add.cell.id=unique(PW_Adu1@meta.data$orig.ident))
PW_Adu1<-RenameCells(object=PW_Adu1,new.names=gsub('-1','',Cells(PW_Adu1)))


myinf18 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/raw/JH_PW_Adu2/sample_feature_bc_matrix/"
PW_Adu2.data <- Read10X(data.dir=myinf18)
PW_Adu2 <- CreateSeuratObject(counts = PW_Adu2.data, project = "PW_Adu2",min.cells=3,min.features=200)
PW_Adu2<-RenameCells(object=PW_Adu2,add.cell.id=unique(PW_Adu2@meta.data$orig.ident))
PW_Adu2<-RenameCells(object=PW_Adu2,new.names=gsub('-1','',Cells(PW_Adu2)))





xx<- merge(PW_C, y = c(Kid1,Kid2,Kid4,Kid6,Kid7,Kid8,Kid9,Kid10,Kid11,PW_Q,PW_S,PW_U,PW_W,PW_Adu1,PW_Adu2), project = "original")

myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/Seruat_PW_PW_Adult_kit_Human__raw.RDS"
 saveRDS(xx, file = myoutf)


xx[["percent.mt"]] <- PercentageFeatureSet(xx, pattern = "^MT-")
myoutf <-paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/QC_control_nGene_UMI_mito_all.pdf")
    pdf(myoutf,width=20,height=5)
     VlnPlot(xx, features = c("nFeature_RNA", "nCount_RNA", "percent.mt",raster=FALSE), ncol = 3)
    dev.off()

     myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/QC_control_correlation_all.pdf"
      pdf(myoutf,width=10,height=5)
      plot1 <- FeatureScatter(xx, feature1 = "nCount_RNA", feature2 = "percent.mt",raster=FALSE)
     plot2 <- FeatureScatter(xx, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",raster=FALSE)
     CombinePlots(plots = list(plot1, plot2))
       dev.off()




PW_Adult_kit.list<-list(PW_C,Kid1,Kid2,Kid4,Kid6,Kid7,Kid8,Kid9,Kid10,Kid11,PW_Q,PW_S,PW_U,PW_W,PW_Adu1,PW_Adu2)

PW_Adult_kit.list <- lapply(X = PW_Adult_kit.list, FUN = function(x) {
    x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")
    x <- subset(x, subset = nCount_RNA<40000 & nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})


############################ Using RPCA to anchor############ 
 features <- SelectIntegrationFeatures(object.list = PW_Adult_kit.list,nfeatures = 3000)
 PW_Adult_kit.list <- lapply(X = PW_Adult_kit.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
 })

 immune.anchors <- FindIntegrationAnchors(object.list = PW_Adult_kit.list, anchor.features = features, reduction = "rpca")
 PW_Adult_kit.integrated <- IntegrateData(anchorset = immune.anchors)
 DefaultAssay(PW_Adult_kit.integrated) <- "integrated"

 PW_Adult_kit.integrated <- ScaleData(PW_Adult_kit.integrated, verbose = FALSE)
  PW_Adult_kit.integrated <- RunPCA(PW_Adult_kit.integrated, npcs = 30, verbose = FALSE)


 x<-c(1:15)


 eigs <- PW_Adult_kit.integrated[["pca"]]@stdev^2 #PCA的STBandard dev的平方#
  print(paste0("Variance captured by 35 PCs: ",sum(eigs[x] / sum(eigs)))) # 0.898#




# PW_Adult_kit.integrated do clustering and TSNE plot No regresSTBissue ---------------------------------------------


 #[5]Do clustering
 PW_Adult_kit.integrated<-FindNeighbors(PW_Adult_kit.integrated, dims = x)
  PW_Adult_kit.integrated <- FindClusters(PW_Adult_kit.integrated,resolution = seq(0.1,1.0,0.1))

PrintFindClustersParams(object = PW_Adult_kit.integrated)

   clusters_resolution<-sapply(grep("^integrated_snn_res",colnames(PW_Adult_kit.integrated@meta.data),value = TRUE),
                            function(x) length(unique(PW_Adult_kit.integrated@meta.data[,x])))
 clusters_resolution


 #PC=x PW_Adult_kit.integrated
#integrated_snn_res.0.1 integrated_snn_res.0.2 integrated_snn_res.0.3 
#                     9                     13                     17 
#integrated_snn_res.0.4 integrated_snn_res.0.5 integrated_snn_res.0.6 
#                    21                     22                     25 
#integrated_snn_res.0.7 integrated_snn_res.0.8 integrated_snn_res.0.9 
#                    28                     31                     31 
#  integrated_snn_res.1 
#                    33 


 #Set up a loop to look for markers

 tar <- c("integrated_snn_res.0.1","integrated_snn_res.0.2","integrated_snn_res.0.3","integrated_snn_res.0.4",
         "integrated_snn_res.0.5","integrated_snn_res.0.6","integrated_snn_res.0.7","integrated_snn_res.0.8",
         "integrated_snn_res.0.9","integrated_snn_res.1")




 PW_Adult_kit.integrated <- RunUMAP(PW_Adult_kit.integrated, dims= x)

 PW_Adult_kit.integrated<-RunTSNE(PW_Adult_kit.integrated,reduction = "pca",dims=x)

 for(i in 1:length(tar))
 {
  Idents(object = PW_Adult_kit.integrated) <- tar[i]
  myoutf1 <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/UMAP_SB6_",tar[i],"_undefined.pdf")
  pdf(myoutf1,7,5)
  print(DimPlot(PW_Adult_kit.integrated, reduction = "umap"))
  dev.off()
 }


 for(i in 1:length(tar))
 {
  Idents(object = PW_Adult_kit.integrated) <- tar[i]
  myoutf1 <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/TSNE_SB6_",tar[i],"_undefined.pdf")
  pdf(myoutf1,7,5)
  print(DimPlot(PW_Adult_kit.integrated, reduction = "tsne"))
  dev.off()
 }


 myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/TSNE_PW_Adult_kit.integrated_groupby_orig.ident.pdf"
  pdf(myoutf,width=8,height=5)
 DimPlot(object = PW_Adult_kit.integrated,group.by='orig.ident',reduction="umap",pt.size=0.6)
 dev.off()



 myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/Seruat_PW_Adult_kit.integrated__finished.RDS"
 saveRDS(PW_Adult_kit.integrated, file = myoutf)

 myinf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/Seruat_PW_Adult_kit.integrated__finished.RDS"
 PW_Adult_kit.integrated<-readRDS(myinf)


tar1 <- c("integrated_snn_res.0.6","integrated_snn_res.0.7","integrated_snn_res.0.8")

 for(i in 1 : length(tar1))
 {
  cat("\r",i)
  Idents(object = PW_Adult_kit.integrated) <- tar1[i]
  tmp <- PW_Adult_kit.integrated
  tmp.markers <- FindAllMarkers(object = tmp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/cluster_markers_",tar1[i],"_ScRNA.xls")
  write.table(tmp.markers,myoutf,sep="\t",quote=F)
  myoutf1 <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/heatmap_markers_",tar1[i],"_.pdf")
  pdf(myoutf1,30,35)
  top20 <- tmp.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
  print(DoHeatmap(tmp, features = top20$gene) + NoLegend())
  dev.off()
 }



#########################Find DEG of each cluster using findconserve markers#########
 tar<-"integrated_snn_res.0.3"
 Idents(object = PW_Adult_kit.integrated) <- tar
 tar1<-unique(as.numeric(PW_Adult_kit.integrated@meta.data$"integrated_snn_res.0.3"))
 tar1<-sort(tar1-1)
 DEG<-data.frame()
 for (i in 1:length(tar1))
 {
  markers <- FindConservedMarkers(PW_Adult_kit.integrated, ident.1 = tar1[i], grouping.var = "orig.ident", verbose = T,only.pos = T,logfc.threshold = 0.25,min.cells.group=1)
  markers[,"cluster"]<-tar1[i]
  markers[,'gene']<-row.names(markers)
  tar2<-grep('log2FC',colnames(markers))
  markers[,'average_Log2FC']<-apply(markers[,colnames(markers)[tar2]],1,mean)
  tag<-markers[,'max_pval']<=0.05
  markers<-markers[tag,]
  markers<-markers[order(markers[,'average_Log2FC'],decreasing=TRUE),]
  DEG<-rbind(DEG,markers)
 }

 myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/cluster_markers_res.0.3_conservemarkers.xls")
  write.table(DEG,myoutf,sep="\t",quote=F)








myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/Seruat_PW_Adult_kit.integrated__finished.RDS"
PW_Adult_kit.integrated<-readRDS(myinf)

 DefaultAssay(PW_Adult_kit.integrated) <- "RNA"

data_RNA<-PW_Adult_kit.integrated[['RNA']]@data
data_integrated<-PW_Adult_kit.integrated[['integrated']]@data
data_RNA<-data_RNA[Immune_general,]
data_integrated<-data_integrated[Immune_general,]
cells<-sample(colnames(data_RNA),20)

Immune_general<-c('CD3G','CD3E','CD8A','CD4','ITGAM','IL13','IL4','SELP','CD68','GATA6','CD1C','LYVE1','CCR2','TIMD4','CD19')



DefaultAssay(PW_Adult_kit.integrated) <- "integrated"
tar1 <- c("integrated_snn_res.0.6")

 for(i in 1 : length(tar1))
 {
  cat("\r",i)
  Idents(object = PW_Adult_kit.integrated) <- tar1[i]
  tmp <- PW_Adult_kit.integrated
  tmp.markers <- FindMarkers(object = tmp, only.pos = FALSE, ident.1='10',min.pct = 0.1, logfc.threshold = 0.1)
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/Cluster_10_DEG.xls")
  write.table(tmp.markers,myoutf,sep="\t",quote=F)
 }



 
DefaultAssay(PW_Adult_kit.integrated) <- "integrated"
tar1 <- c("integrated_snn_res.0.6")

 for(i in 1 : length(tar1))
 {
  cat("\r",i)
  Idents(object = PW_Adult_kit.integrated) <- tar1[i]
  tmp <- PW_Adult_kit.integrated
  tmp.markers <- FindMarkers(object = tmp, only.pos = FALSE, ident.1='1',ident.2='4',min.pct = 0.1, logfc.threshold = 0.1)
  myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.6/cluster1_vs_cluster4.xls"
  write.table(tmp.markers,myoutf,sep="\t",quote=F)
 }




Immune_general<-c('CD3G','CD3E','CD8A','CD4','ITGAM','ITGAX','MRC1','CD14','CD5','GATA6','CD1C','LYVE1','CCR2','TIMD4','CD19')


p<-FeaturePlot(object = PW_Adult_kit.integrated, features = Immune_general,cols= c('#deebf7', "#bd0026"),
               reduction = "umap",slot = "data",max.cutoff=4,pt.size=0.3,combine=F)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + theme(plot.title = element_text(size=40),
                          axis.title = element_text(family="Helvetica",face = "bold", color = "black",size=20), 
                           axis.text=element_text(family="Helvetica",face = "bold", color = "black",size=16),
                           axis.line = element_line(colour = 'black', size = 1),
                           axis.ticks = element_line(colour = "black", size = 1),
                           axis.ticks.length=unit(.3, "cm"))
}

myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/Immune_general_marker.pdf"
pdf(myoutf,width=50,height=30)
print(cowplot::plot_grid(plotlist = p,ncol=5))
dev.off()


T_general<-c('CD3G','CD3E','CD8A','CD4','ITGAM','ITGAX','NKG7','KLRG1','KLRB1','CD14','GATA6','CD1C','LYVE1','CCR2','TIMD4','CD19')


tar<-"integrated_snn_res.0.6"
 Idents(object = PW_Adult_kit.integrated) <- tar
myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/T_general_violin_res.0.6.pdf"
  pdf(myoutf,width=40,height=25)
 VlnPlot(object = PW_Adult_kit.integrated,assay='RNA',features = T_general,pt.size=0.2)
 dev.off()





tag<-PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.6=='7'

p<-FeaturePlot(object = PW_Adult_kit.integrated,cells=row.names(PW_Adult_kit.integrated@meta.data)[tag],
 features = c('GATA6', 'SELP'),blend = TRUE, blend.threshold = 0.000000001,pt.size = 1.5)

myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/GATA6_SELP_coexpression.pdf"
pdf(myoutf,width=40,height=12)
print(p)
dev.off()


myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/GATA6_SELP_coexpression.pdf"
pdf(myoutf,width=40,height=12)
FeaturePlot(PW_Adult_kit.integrated,cells=row.names(PW_Adult_kit.integrated@meta.data)[tag],
 features = c("GATA6", "SELP"), order = T, min.cutoff = 0,max.cutoff = 1, 
 blend = T, cols = c("darkblue", "green", "magenta"), blend.threshold = 0,pt.size = 1.5) &DarkTheme()
dev.off()







Idents(object = PW_Adult_kit.integrated) <- 'integrated_snn_res.0.6'
levels(PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.6)<-seq(0,24,1)
  myoutf1 <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.6/UMAP_recolored.pdf"
  pdf(myoutf1,7,5)
  print(DimPlot(PW_Adult_kit.integrated, reduction = "umap",cols=c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
            "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
            "#8dd3c7", "#bebada", "#fb8072", "#80b1d3", "#fdb462",
            "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5",
            "#ffed6f", "#a65628", "#f781bf", "#6d9eeb", "#ff9896")))
  dev.off()



#############Single cell Heatmap Plot Seurat##################

colors=c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
            "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
            "#8dd3c7", "#bebada", "#fb8072", "#80b1d3", "#fdb462",
            "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5",
            "#ffed6f", "#a65628", "#f781bf", "#6d9eeb", "#ff9896")


myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/Seruat_PW_Adult_kit.integrated__finished.RDS"
PW_Adult_kit.integrated<-readRDS(myinf)

Idents(object = PW_Adult_kit.integrated) <- 'integrated_snn_res.0.6'
  tmp <- PW_Adult_kit.integrated
  tmp.markers <- FindAllMarkers(object = PW_Adult_kit.integrated, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.3)


tar<-sort(unique(PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.6))
all_cell<-character()
for(i in 1:(length(tar)-3))
{
  cat("\r",i)
  len<-sum(PW_Normal@meta.data$RNA_snn_res.0.3==tar[i])
  tag<-sample(1:len,230)
  info<-PW_Normal@meta.data[PW_Normal$RNA_snn_res.0.3==tar[i],]
  barcode_use<-row.names(info)[tag]
  all_cell<-c(all_cell,barcode_use)
}
cluster15_16_barcode<-row.names(PW_Normal@meta.data)[PW_Normal@meta.data$RNA_snn_res.0.3%in%c('8','9','10')]

all_cell<-c(all_cell,cluster15_16_barcode)

myoutf1 <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PWC_Nan_For_Figure1/log_norm/PCx/heatmap_marker_genes_ordered_245cells_each_cluster.pdf")

tmp.markers <- tmp.markers[tmp.markers$p_val_adj<0.05,]

ts<-tmp.markers[order(factor(tmp.markers$cluster, levels = tar)),]

top20 <- ts %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

PW_Normal@meta.data[,"ordered_ident"]<-factor(PW_Normal@meta.data$RNA_snn_res.0.3,
                                                                    levels=seq(0,11,1))

pdf(myoutf1,13,18)
print(DoHeatmap(PW_Normal, cells=all_cell,features = top20$gene,group.by="RNA_snn_res.0.3",assay='RNA') + 
        scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')), 
                              mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')), 
                              midpoint = 0, guide = "colourbar", aesthetics = "fill")+NoLegend())
dev.off()










############Group#########################
PW_Adult_kit.integrated@meta.data[,'Group']<-rep('NA',nrow(PW_Adult_kit.integrated@meta.data))
tag1<-grep('PW',PW_Adult_kit.integrated@meta.data$orig.ident)
tag2<-grep('Kid',PW_Adult_kit.integrated@meta.data$orig.ident)
PW_Adult_kit.integrated@meta.data$Group[tag1]<-'Adult'
PW_Adult_kit.integrated@meta.data$Group[tag2]<-'Pediatric'


##############Check SELP, TIMD4 expression profile############################
myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/Seruat_PW_Adult_kit.integrated__finished.RDS"
PW_Adult_kit.integrated<-readRDS(myinf)

data<-as.matrix(GetAssayData(PW_Adult_kit.integrated,assay='RNA'))
PW_Adult_kit.integrated@meta.data[,'SELP']<-rep('0',nrow(PW_Adult_kit.integrated@meta.data))
PW_Adult_kit.integrated@meta.data[,'TIMD4']<-rep('0',nrow(PW_Adult_kit.integrated@meta.data))
PW_Adult_kit.integrated@meta.data[,'GATA6']<-rep('NA',nrow(PW_Adult_kit.integrated@meta.data))

tag1<-data['SELP',]>0
tag2<-data['TIMD4',]>0
tag3<-data['GATA6',]>0

sum(tag1)
sum(tag2)
sum(tag3)

sum(tag1&tag3)





name_1<-colnames(data)[tag1]
name_2<-colnames(data)[tag2]
name_3<-colnames(data)[tag3]

tag<-PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.6%in%c('11','2')
tag4<-PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.6%in%c('7')


PW_Adult_kit.integrated@meta.data[name_1,'SELP']<-1
PW_Adult_kit.integrated@meta.data[name_2,'TIMD4']<-1
PW_Adult_kit.integrated@meta.data$GATA6[tag]<-'GATA6-'
PW_Adult_kit.integrated@meta.data$GATA6[tag3&tag4&tag1]<-'GATA6+'


table(PW_Adult_kit.integrated@meta.data$SELP,PW_Adult_kit.integrated@meta.data$orig.ident)
table(PW_Adult_kit.integrated@meta.data$TIMD4,PW_Adult_kit.integrated@meta.data$orig.ident)
table(PW_Adult_kit.integrated@meta.data$GATA6,PW_Adult_kit.integrated@meta.data$orig.ident)




 
DefaultAssay(PW_Adult_kit.integrated) <- "integrated"
tar1 <- c("GATA6")

 for(i in 1 : length(tar1))
 {
  cat("\r",i)
  Idents(object = PW_Adult_kit.integrated) <- tar1[i]
  tmp <- PW_Adult_kit.integrated
  tmp.markers <- FindMarkers(object = tmp, only.pos = FALSE, ident.1='GATA6+',ident.2='GATA6-',min.pct = 0.001, logfc.threshold = 0.001)
  myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.6/GATA6vsGATA6-_DEG_GSEA.xls"
  write.table(tmp.markers,myoutf,sep="\t",quote=F)
 }




myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.6/GATA6vsGATA6-_DEG_GSEA_groupby_orig.ident.pdf"
  pdf(myoutf,width=8,height=5)
 DimPlot(object = PW_Adult_kit.integrated,reduction="umap",pt.size=0.6)
 dev.off()






  tmpinf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.6/GATA6vsGATA6-_DEG_GSEA.xls"
  res <- read.table(tmpinf,sep="\t",quote=NULL)
  res <- res[order(res$avg_log2FC,decreasing=T),]
  
  xx <- res$avg_log2FC
  names(xx) <- row.names(res)
  names(xx) <- toupper(names(xx))
  
  tag = duplicated(names(xx))
  res = res[tag==0,]
  row.names(res) <- toupper(row.names(res))
  res[,"Name"] <- row.names(res)
  res[,"metric"] <- res$avg_log2FC
  res <- res[,c("Name","metric")]
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.6/GATA6vsGATA6-.rnk")
  write.table(res,file=myoutf,quote=F,sep="\t",row.names=F)
}





dir <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_ADULT_Kid_normal_Human/PCx/res.0.5/GSEA/marker/integrated/"
files <- list.files(dir)
file_nam <- gsub(".xls","",files)
for(i in 1:length(files))
{
  cat("\r",i)
  tmpinf <- paste0(dir,files[i])
  res <- read.table(tmpinf,sep="\t",quote=NULL)
  #tag<-res$average_Log2FC>0
  #res<-res[tag,]
  res <- res[order(res$avg_log2FC,decreasing=T),]
  
  xx <- res$avg_log2FC
  names(xx) <- row.names(res)
  names(xx) <- toupper(names(xx))
  
  tag = duplicated(names(xx))
  res = res[tag==0,]
  row.names(res) <- toupper(row.names(res))
  res[,"Name"] <- row.names(res)
  res[,"metric"] <- res$avg_log2FC
  res <- res[,c("Name","metric")]
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_ADULT_Kid_normal_Human/PCx/res.0.5/GSEA/preranked_GSEA/",file_nam[i],".rnk")
  write.table(res,file=myoutf,quote=F,sep="\t",row.names=F)
}


########3Feature plot#################

myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_ADULT_Kid_normal_Human/Seruat_PW_Adult_kit.integrated__finished.RDS"
PW_Adult_kit.integrated<-readRDS(myinf)

 DefaultAssay(PW_Adult_kit.integrated) <- "RNA"

data_RNA<-PW_Adult_kit.integrated[['RNA']]@data
data_integrated<-PW_Adult_kit.integrated[['integrated']]@data
data_RNA<-data_RNA[Immune_general,]
data_integrated<-data_integrated[Immune_general,]
cells<-sample(colnames(data_RNA),20)


 DefaultAssay(PW_Adult_kit.integrated) <- "RNA"
Immune_general<-c('CD3G','CD3E','CD8A','CD4','ITGAM','ITGAX','MRC1','CD14','CD68','GATA6','CD1C','LYVE1','CCR2','TIMD4','CD19','TLN1')

p<-FeaturePlot(object = PW_Adult_kit.integrated, features = Immune_general,cols= c('#deebf7', "#bd0026"),
               reduction = "umap",slot = "data",max.cutoff=4,pt.size=0.3,combine=F)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + theme(plot.title = element_text(size=40),
                          axis.title = element_text(family="Helvetica",face = "bold", color = "black",size=20), 
                           axis.text=element_text(family="Helvetica",face = "bold", color = "black",size=16),
                           axis.line = element_line(colour = 'black', size = 1),
                           axis.ticks = element_line(colour = "black", size = 1),
                           axis.ticks.length=unit(.3, "cm"))
}
 
myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_ADULT_Kid_normal_Human/PCx/Immune_general_marker.pdf"
pdf(myoutf,width=50,height=30)
print(cowplot::plot_grid(plotlist = p,ncol=5))
dev.off()






tar1 <- c("integrated_snn_res.0.3","integrated_snn_res.0.4")

 for(i in 1 : length(tar1))
 {
  cat("\r",i)
  Idents(object = PW_Adult_kit.integrated) <- tar1[i]
  tmp <- PW_Adult_kit.integrated
  tmp.markers <- FindAllMarkers(object = tmp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_ADULT_Kid_normal_Human/PCx/cluster_markers_",tar1[i],"_ScRNA_RNAassay.xls")
  write.table(tmp.markers,myoutf,sep="\t",quote=F)
 }


myeloid_general<-c('SELL','IRF4','F5','CSF1R','VCAN','MERTK','MRC1','CD14','ITGAX','ITGAM','GATA6','CD1C','LYVE1','CCR2','TIMD4','SELP','CX3CR1','CD163',
'FLT3','MARCO','IL3RA','ICAM2','SIGLEC6','RETN','CD226',"CLEC4C","CLEC9A",'CD209','VSIG4','CD36','FCGR1A','MERTK','FCER1A','CD80','CD86','CLEC10A')

flow_panel<-c('FCGR1A','CD14','ITGAM','ITGAX','LYVE1','TIMD4','CD226','SELP','CCR2','HLA-DRA','CCR7','DPP4','MARCO','MRC1','FOLR2',
'LILRB4','CD1C','CLEC9A','SIGLEC6','C5AR1','IL3RA','CD274','FCER1A','SIRPA','CD163')

important<-c('GATA6','SELP','TIMD4','LYVE1','MARCO','CCR2','CD1C','CD14','ICAM2','IL3RA','MRC1','CD163')

paper<-c('GATA6','SELP','TIMD4','LYVE1','MARCO','MRC1','CCR2','CD1C','CD14')


p<-FeaturePlot(object = PW_Adult_kit.integrated, features = paper,cols= c('#deebf7', "#bd0026"),
               reduction = "umap",slot = "data",max.cutoff=4,pt.size=4,combine=F,order = TRUE)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + theme(plot.title = element_text(size=50),
                          axis.title = element_text(family="Helvetica",face = "bold", color = "black",size=40), 
                           axis.text=element_text(family="Helvetica",face = "bold", color = "black",size=30),
                           axis.line = element_line(colour = 'black', size = 3),
                           axis.ticks = element_line(colour = "black", size = 3),
                           axis.ticks.length=unit(1, "cm"))
}
 
myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_ADULT_Kid_normal_Human/PCx/paper_marker.pdf"
pdf(myoutf,width=50,height=40)
print(cowplot::plot_grid(plotlist = p,ncol=3))
dev.off()




p<-FeaturePlot(object = PW_Adult_kit.integrated, features = myeloid_general,cols= c('#deebf7', "#bd0026"),
               reduction = "umap",slot = "data",max.cutoff=4,pt.size=0.3,combine=F,order = TRUE)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + theme(plot.title = element_text(size=40),
                          axis.title = element_text(family="Helvetica",face = "bold", color = "black",size=20), 
                           axis.text=element_text(family="Helvetica",face = "bold", color = "black",size=16),
                           axis.line = element_line(colour = 'black', size = 1),
                           axis.ticks = element_line(colour = "black", size = 1),
                           axis.ticks.length=unit(.3, "cm"))
}
 
myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_ADULT_Kid_normal_Human/PCx/Myeloid_marker.pdf"
pdf(myoutf,width=80,height=45)
print(cowplot::plot_grid(plotlist = p,ncol=8))
dev.off()




p<-FeaturePlot(object = PW_Adult_kit.integrated, features = flow_panel,cols= c('#deebf7', "#bd0026"),
               reduction = "umap",slot = "data",max.cutoff=4,pt.size=0.3,combine=F,order = TRUE)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + theme(plot.title = element_text(size=40),
                          axis.title = element_text(family="Helvetica",face = "bold", color = "black",size=20), 
                           axis.text=element_text(family="Helvetica",face = "bold", color = "black",size=16),
                           axis.line = element_line(colour = 'black', size = 1),
                           axis.ticks = element_line(colour = "black", size = 1),
                           axis.ticks.length=unit(.3, "cm"))
}
 
myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_ADULT_Kid_normal_Human/PCx/Flow_panel_marker.pdf"
pdf(myoutf,width=50,height=50)
print(cowplot::plot_grid(plotlist = p,ncol=5))
dev.off()




p<-FeaturePlot(object = PW_Adult_kit.integrated, features = important,cols= c('#deebf7', "#bd0026"),
               reduction = "umap",slot = "data",max.cutoff=4,pt.size=0.3,combine=F,order = TRUE)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + theme(plot.title = element_text(size=50),
                          axis.title = element_text(family="Helvetica",face = "bold", color = "black",size=30), 
                           axis.text=element_text(family="Helvetica",face = "bold", color = "black",size=20),
                           axis.line = element_line(colour = 'black', size = 1),
                           axis.ticks = element_line(colour = "black", size = 1),
                           axis.ticks.length=unit(.3, "cm"))
}
 
myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_ADULT_Kid_normal_Human/PCx/important_marker.pdf"
pdf(myoutf,width=65,height=50)
print(cowplot::plot_grid(plotlist = p,ncol=4))
dev.off()




#############Cluster contribution by each different patient##########
   
   Idents(PW_Adult_kit.integrated)<-"integrated_snn_res.0.6"
 
 xx<-table(PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.3,PW_Adult_kit.integrated@meta.data$orig.ident)
 yy<-t(apply(xx,1,function(x){x/sum(x)})   )
 myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_ADULT_Kid_normal_Human/PCx/res.0.2/cluster_patient_distribution/patient_distribution_to_each_cluster.xls"
 write.table(yy,myoutf,sep="\t",quote=F)
 
 data<-matrix(0,ncol(yy)*nrow(yy),3)
 colnames(data)<-c("cluster","patient","value")

 data<-as.data.frame(data)
 data$cluster<-rep(row.names(yy),ncol(yy))

 data$patient<-as.character(sapply(colnames(yy),function(x){rep(x,nrow(yy))}))

 for (i in 1:nrow(data))
 {
  number<-yy[data$cluster[i],data$patient[i]]
  data$value[i]<-number
 }


 data$cluster<-factor(data$cluster,levels=unique(data$cluster))
 data$patient<-factor(data$patient,levels=c("PW_normal","PW_CRC1","PW_CRC2",'Kid1','Kid2','Kid4','Kid6','Kid7','Kid8','Kid9','Kid10'))

  col<-c("#2171b5","#cb181d","#74c476")

 myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_ADULT_Kid_normal_Human/PCx/res.0.3/cluster_patient_distribution/stackplot_patient_contribution_each_cluster_proportion.pdf"
 pdf(myoutf,width=8,height=8)
 ggplot() + geom_bar(aes(y = value, x = cluster, fill = patient), data = data,
                    stat="identity")+theme(axis.text.x = element_text(angle = 60, hjust = 1,size=10),
                                           axis.title.x = element_blank())

 dev.off()






myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_ADULT_Kid_normal_Human/Seruat_PW_Adult_kit.integrated__finished.RDS"
PW_Adult_kit.integrated<-readRDS(myinf)

tar<-unique(PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.6)
tar<-tar[order(tar)]
chi<-as.matrix(table(PW_Adult_kit.integrated@meta.data$orig.ident,PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.6))

chi_square<-chisq.test(chi)

chi_observed<-chi_square$observed
chi_expected<-chi_square$expected
ROE<-chi_observed/chi_expected

myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_ADULT_Kid_normal_Human/PCx/res.0.5/cluster_patient_distribution/ROE.xls"
write.table(ROE,myoutf,sep="\t",quote=F)



each_cluster_all<-apply(chi,2,sum)
chi_enrichment<-chi
for (i in 1:nrow(chi))
{
    for (k in 1:ncol(chi))
    {
      chi_enrichment[i,k]<-chi[i,k]*sum(chi)/sum(chi[i,])/sum(chi[,k])
    }
}

myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_ADULT_Kid_normal_Human/PCx/res.0.5/cluster_patient_distribution/hyper_enrichment.xls"
write.table(chi_enrichment,myoutf,sep="\t",quote=F)











myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_ADULT_Kid_normal_Human/Seruat_PW_Adult_kit.integrated__finished.RDS"
PW_Adult_kit.integrated<-readRDS(myinf)

tar<-unique(PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.6)
tar<-tar[order(tar)]
chi<-as.matrix(table(PW_Adult_kit.integrated@meta.data$orig.ident,PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.6))

chi_square<-chisq.test(chi)

chi_observed<-chi_square$observed
chi_expected<-chi_square$expected
ROE<-chi_observed/chi_expected

myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_ADULT_Kid_normal_Human/PCx/res.0.5/cluster_patient_distribution/ROE.xls"
write.table(ROE,myoutf,sep="\t",quote=F)



each_cluster_all<-apply(chi,2,sum)
chi_enrichment<-chi
for (i in 1:nrow(chi))
{
    for (k in 1:ncol(chi))
    {
      chi_enrichment[i,k]<-chi[i,k]*sum(chi)/sum(chi[i,])/sum(chi[,k])
    }
}

myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_ADULT_Kid_normal_Human/PCx/res.0.7/cluster_patient_distribution/hyper_enrichment.xls"
write.table(chi_enrichment,myoutf,sep="\t",quote=F)




#################Output cell of each orig.ident and adult versus kid#############
myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/Seruat_PW_Adult_kit.integrated__finished.RDS"
PW_Adult_kit.integrated<-readRDS(myinf)

Idents(PW_Adult_kit.integrated)<-"integrated_snn_res.0.6"
 tar<-unique(PW_Adult_kit.integrated@meta.data$Group)
 for (i in 1:length(tar))
 {
    tag<-PW_Adult_kit.integrated@meta.data$Group==tar[i]
    cell_ID<-row.names(PW_Adult_kit.integrated@meta.data)[tag]
    myoutf1 <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.6/each_sample/UMAP_SB6_",tar[i],"_undefined.pdf")
  pdf(myoutf1,7,5)
  print(DimPlot(PW_Adult_kit.integrated, reduction = "umap",pt.size=0.1,cells =cell_ID))
  dev.off()
 
 }
  
  



myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_ADULT_Kid_normal_Human/Seruat_PW_Adult_kit.integrated__finished.RDS"
PW_Adult_kit.integrated<-readRDS(myinf)

Idents(PW_Adult_kit.integrated)<-"integrated_snn_res.0.6"
cell_ID_adult<-row.names(PW_Adult_kit.integrated@meta.data)[as.numeric(PW_Adult_kit.integrated@meta.data$age)>18]
cell_ID_kid<-row.names(PW_Adult_kit.integrated@meta.data)[as.numeric(PW_Adult_kit.integrated@meta.data$age)<18]


    myoutf1 <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_ADULT_Kid_normal_Human/PCx/res.0.5/each_sample/All_Adult.pdf")
  pdf(myoutf1,7,5)
  print(DimPlot(PW_Adult_kit.integrated, reduction = "umap",pt.size=0.1,cells =cell_ID_adult))
  dev.off()
 

myoutf1 <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_ADULT_Kid_normal_Human/PCx/res.0.5/each_sample/All_Kid.pdf")
  pdf(myoutf1,7,5)
  print(DimPlot(PW_Adult_kit.integrated, reduction = "umap",pt.size=0.1,cells =cell_ID_kid))
  dev.off()
 
  
  

##########First get average gene expression of each cluster and then do z transformation to do heatmap for selected genes###########

myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/Seruat_PW_Adult_kit.integrated__finished.RDS"
exp<-readRDS(myinf)
DefaultAssay(exp) <- "RNA"

myoutf1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.6/expression_data/CD8_normalized_data.xls"
myoutf2 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.6/expression_data/CD8_normalized_cell_info.xls"


info <- exp@meta.data
res <- GetAssayData(object = exp,slot='data')
res <- as.matrix(res)
row.names(res) <- toupper(row.names(res))
avg <- apply(res,1, mean)
tag <- avg >0
res <- res[tag,]

res <- unique(res)

tag <- duplicated(row.names(res))

dup_gene <- row.names(res)[tag]

tag <- row.names(res) == dup_gene
#tag==logical(0) Thus ignore the res1,res STBep
#res1 <- res[tag==0,]

#res2 <- res[tag==1,]
#avg <- apply(res2,1, mean)

#res <- rbind(res1, res2[2,])
res <- as.data.frame(res)

write.table(res,myoutf1,sep="\t",quote=F)
write.table(info,myoutf2,sep="\t",quote=F)






myinf1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_ADULT_Kid_normal_Human/PCx/res.0.5/expression_data/CD8_normalized_data.xls"
myinf2 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_ADULT_Kid_normal_Human/PCx/res.0.5/expression_data/CD8_normalized_cell_info.xls"

data<-read.table(myinf1,sep="\t",quote=NULL)
info<-read.table(myinf2,sep="\t",quote=NULL)

DC_markers<-c('TCF4','IL3RA','CLEC4C','LILRB4','FCER1A','FCER1G','LILRA4',
            'IRF8','IRF4','ID2','BATF3','CLEC9A','XCR1','BTLA',
            'CD163','S100A9','S100A8','VCAN','FCN1','CD14','BST1','CD36',
            'FCGR2B','CD1C','HLA-DQA2','HLA-DQA1','HLA-DQB1','HLA-DQP1',
            'ZEB2','KLF4','ITGAX','ITGAM','CD2','SIRPA','CLEC4A','CLEC10A',
            'AXL','SIGLEC6','CX3CR1','SIGLEC1')

PM_markers<-c('GARNL3','CAV1','HGF','LAPTM5','GATA6','CD9','TGFB2','APOC1','F5','HDC','LYZ1',
               'FSCN1','ID1','PLXDC2','SELP','CXCL13','FABP5','ADGRE1','GAS6','ID2',
               'MRC1','FOLR2','LYVE1','MARCO','CD163','CD209','CD74','FPR1',
               'ANXA5','CBR2','CLEC4A1','MAF','ITGB1','CXCL13',
               'FABP7','CCR2','CSF1R','SELL','ITGAX','HLA-DQA2','HLA-DQA1','HLA-DQB1','HLA-DQP1','ITGAM',
               'RCN3','FCGRT','CCL24','IRF7','EMP1','LGALS1','HAL','RPS18')

Paper_markers<-c('GATA6','SELP','TIMD4','LYVE1','MARCO','MRC1','CCR2','CD1C','CD14','CLEC9A')

#cluster_markers<-c('S100A9','S100A8','VCAN','FCN1','CD14','MAFB','ITGAM','SPP1',
#                    'FCER1A','CD1C','CD1E','HLA-DPB1','HLA-DQB1','HLA-DRA','CLEC10A','FCGR2B','CLEC4A','TCF4','IL3RA','ITGAX','ZEB2',
#                    'LILRB4','CSF1R','GSN','GPR183','RGS1','IL1B','FOS','NR4A1',
#                    'C1QC','C1QA','C1QB','FN1','LYVE1','FOLR2','MARCO','CD163','FCGR3A','TIMD4','APOE','CD68',
#                    'IL32','CCL2','GZMA','CCL5','NKG7','GNLY','ID3','STAT1','CD36','CCL20','ITGA4','CD3G','CD8A',
#                    'CCL4','CXCL10','CCL3','CXCL9','CXCL11',
#                    'CXCR6','GATA3','RORA','TIGIT','CD69',
#                    'MKI67','TK1','CCNA2','TUBA1B',
#                    'IRF8','ID2','BATF3','CLEC9A','XCR1','BTLA',
#                     'NEAT1','CIITA','NFKBIZ','FOXP1','TLR2','IRF4'
#                    )



data_heatmap_paper<-data[Paper_markers,]
data_heatmap_paper<-na.omit(data_heatmap_paper)

tar <- unique(info$integrated_snn_res.0.6)
#tar <- as.numeric(tar)
tar <- tar[order(tar)]

for(i in 1 : length(tar))
{
  res <- matrix(0, nrow(data_heatmap_paper),4)
  row.names(res) <- row.names(data_heatmap_paper)
  colnames(res) <- c("C1_avg_expression","C1_median_expression","C2_avg_expression","C2_median_expression")
  
  tag <- info$integrated_snn_res.0.6 == tar[i]
  
  sam1 <- row.names(info)[tag==1]
  sam2 <- row.names(info)[tag==0]
  
  for(k in 1 : nrow(data_heatmap_paper))
  {
    cat("\r",i,"----->",k)
    xx1 <- as.numeric(data_heatmap_paper[k, sam1])
    xx2 <- as.numeric(data_heatmap_paper[k, sam2])
    
    res[k,"C1_avg_expression"] <- mean(xx1)
    res[k,"C2_avg_expression"] <- mean(xx2)
    res[k,"C1_median_expression"]<-median(xx1)
    res[k,"C2_median_expression"] <- median(xx2)
  }
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_ADULT_Kid_normal_Human/PCx/res.0.5/cluster_identity/pheatmap/Paper_Heatmap_",tar[i],"_vs_others_normalized_expression.xls")
  write.table(res,myoutf,sep="\t",quote=F)
}




res<-matrix(0,nrow(data_heatmap_paper),length(tar))
 row.names(res)<-row.names(data_heatmap_paper)
 colnames(res)<-tar[order(tar)]



for (i in 1:length(tar))
{
  myinf<-paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_ADULT_Kid_normal_Human/PCx/res.0.5/cluster_identity/pheatmap/Paper_Heatmap_",tar[i],"_vs_others_normalized_expression.xls")
  xx<-read.table(myinf,sep='\t')
  res[,i]<-xx$C1_avg_expression
}

res <- apply(res,1, function(arg) (arg-mean(arg))/sd(arg))
res <- t(res)


#res<-res[,c("2",'5','4','8','9','0','1','7','3','6')]

#ts<-res[sort(row.names(res)),]
breaklist_1<-seq(-2.5,0,0.2)
breaklist_2<-seq(0,2.6,0.2)
breaklist<-unique(c(breaklist_1,breaklist_2))

myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_ADULT_Kid_normal_Human/PCx/res.0.5/cluster_identity/pheatmap/heatmap/Paper_Heatmap_gene_each_cluster.pdf")
pdf(myoutf,width=20,height=5)
pheatmap(mat=res,
         color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(length(breaklist)),
         breaks = breaklist,
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         border_color      = NA,
         show_colnames     = T,
         show_rownames     = T,
         drop_levels       = TRUE,
         scale='none',
         
         fontsize=20)
dev.off()








data_DC<-data[DC_markers,]
data_DC<-na.omit(data_DC)

tar <- unique(info$integrated_snn_res.0.6)
#tar <- as.numeric(tar)
tar <- tar[order(tar)]

for(i in 1 : length(tar))
{
  res <- matrix(0, nrow(data_DC),4)
  row.names(res) <- row.names(data_DC)
  colnames(res) <- c("C1_avg_expression","C1_median_expression","C2_avg_expression","C2_median_expression")
  
  tag <- info$integrated_snn_res.0.6 == tar[i]
  
  sam1 <- row.names(info)[tag==1]
  sam2 <- row.names(info)[tag==0]
  
  for(k in 1 : nrow(data_DC))
  {
    cat("\r",i,"----->",k)
    xx1 <- as.numeric(data_DC[k, sam1])
    xx2 <- as.numeric(data_DC[k, sam2])
    
    res[k,"C1_avg_expression"] <- mean(xx1)
    res[k,"C2_avg_expression"] <- mean(xx2)
    res[k,"C1_median_expression"]<-median(xx1)
    res[k,"C2_median_expression"] <- median(xx2)
  }
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_ADULT_Kid_normal_Human/PCx/res.0.5/cluster_identity/pheatmap/Human_DC_gene_",tar[i],"_vs_others_normalized_expression.xls")
  write.table(res,myoutf,sep="\t",quote=F)
}




res<-matrix(0,nrow(data_DC),length(tar))
 row.names(res)<-row.names(data_DC)
 colnames(res)<-tar[order(tar)]



for (i in 1:length(tar))
{
  myinf<-paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_ADULT_Kid_normal_Human/PCx/res.0.5/cluster_identity/pheatmap/Human_DC_gene_",tar[i],"_vs_others_normalized_expression.xls")
  xx<-read.table(myinf,sep='\t')
  res[,i]<-xx$C1_avg_expression
}

res <- apply(res,1, function(arg) (arg-mean(arg))/sd(arg))
res <- t(res)


#res<-res[,c("2",'5','4','8','9','0','1','7','3','6')]

#ts<-res[sort(row.names(res)),]
breaklist_1<-seq(-2.5,0,0.2)
breaklist_2<-seq(0,2.6,0.2)
breaklist<-unique(c(breaklist_1,breaklist_2))

myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_ADULT_Kid_normal_Human/PCx/res.0.5/cluster_identity/pheatmap/heatmap/Human_DC_gene_each_cluster.pdf")
pdf(myoutf,width=25,height=15)
pheatmap(mat=res,
         color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(length(breaklist)),
         breaks = breaklist,
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         border_color      = NA,
         show_colnames     = T,
         show_rownames     = T,
         drop_levels       = TRUE,
         scale='none',
         
         fontsize=20)
dev.off()












data_PM<-data[PM_markers,]
data_PM<-na.omit(data_PM)

tar <- unique(info$integrated_snn_res.0.6)
#tar <- as.numeric(tar)
tar <- tar[order(tar)]

for(i in 1 : length(tar))
{
  res <- matrix(0, nrow(data_PM),4)
  row.names(res) <- row.names(data_PM)
  colnames(res) <- c("C1_avg_expression","C1_median_expression","C2_avg_expression","C2_median_expression")
  
  tag <- info$integrated_snn_res.0.6 == tar[i]
  
  sam1 <- row.names(info)[tag==1]
  sam2 <- row.names(info)[tag==0]
  
  for(k in 1 : nrow(data_PM))
  {
    cat("\r",i,"----->",k)
    xx1 <- as.numeric(data_PM[k, sam1])
    xx2 <- as.numeric(data_PM[k, sam2])
    
    res[k,"C1_avg_expression"] <- mean(xx1)
    res[k,"C2_avg_expression"] <- mean(xx2)
    res[k,"C1_median_expression"]<-median(xx1)
    res[k,"C2_median_expression"] <- median(xx2)
  }
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_ADULT_Kid_normal_Human/PCx/res.0.5/cluster_identity/pheatmap/Mouse_Peritoneal_MAC_gene_",tar[i],"_vs_others_normalized_expression.xls")
  write.table(res,myoutf,sep="\t",quote=F)
}




res<-matrix(0,nrow(data_PM),length(tar))
 row.names(res)<-row.names(data_PM)
 colnames(res)<-tar[order(tar)]



for (i in 1:length(tar))
{
  myinf<-paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_ADULT_Kid_normal_Human/PCx/res.0.5/cluster_identity/pheatmap/Mouse_Peritoneal_MAC_gene_",tar[i],"_vs_others_normalized_expression.xls")
  xx<-read.table(myinf,sep='\t')
  res[,i]<-xx$C1_avg_expression
}

res <- apply(res,1, function(arg) (arg-mean(arg))/sd(arg))
res <- t(res)


#res<-res[,c("2",'5','4','8','9','0','1','7','3','6')]

#ts<-res[sort(row.names(res)),]
breaklist_1<-seq(-2.5,0,0.2)
breaklist_2<-seq(0,2.6,0.2)
breaklist<-unique(c(breaklist_1,breaklist_2))

myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_ADULT_Kid_normal_Human/PCx/res.0.5/cluster_identity/pheatmap/heatmap/Mouse_Peritoneal_MAC_gene_each_cluster.pdf")
pdf(myoutf,width=15,height=15)
pheatmap(mat=res,
         color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(length(breaklist)),
         breaks = breaklist,
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         border_color      = NA,
         show_colnames     = T,
         show_rownames     = T,
         drop_levels       = TRUE,
         scale='none',
         
         fontsize=20)
dev.off()










data_each_cluster<-data[cluster_markers,]
data_each_cluster<-na.omit(data_each_cluster)

tar <- unique(info$integrated_snn_res.0.3)
#tar <- as.numeric(tar)
tar <- tar[order(tar)]

for(i in 1 : length(tar))
{
  res <- matrix(0, nrow(data_each_cluster),4)
  row.names(res) <- row.names(data_each_cluster)
  colnames(res) <- c("C1_avg_expression","C1_median_expression","C2_avg_expression","C2_median_expression")
  
  tag <- info$integrated_snn_res.0.3 == tar[i]
  
  sam1 <- row.names(info)[tag==1]
  sam2 <- row.names(info)[tag==0]
  
  for(k in 1 : nrow(data_each_cluster))
  {
    cat("\r",i,"----->",k)
    xx1 <- as.numeric(data_each_cluster[k, sam1])
    xx2 <- as.numeric(data_each_cluster[k, sam2])
    
    res[k,"C1_avg_expression"] <- mean(xx1)
    res[k,"C2_avg_expression"] <- mean(xx2)
    res[k,"C1_median_expression"]<-median(xx1)
    res[k,"C2_median_expression"] <- median(xx2)
  }
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_ADULT_Kid_normal_Human/PCx/res.0.3/cluster_identity/pheatmap/marker_each_cluster_",tar[i],"_vs_others_normalized_expression.xls")
  write.table(res,myoutf,sep="\t",quote=F)
}




res<-matrix(0,nrow(data_each_cluster),length(tar))
 row.names(res)<-row.names(data_each_cluster)
 colnames(res)<-tar[order(tar)]



for (i in 1:length(tar))
{
  myinf<-paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_ADULT_Kid_normal_Human/PCx/res.0.3/cluster_identity/pheatmap/marker_each_cluster_",tar[i],"_vs_others_normalized_expression.xls")
  xx<-read.table(myinf,sep='\t')
  res[,i]<-xx$C1_avg_expression
}

res <- apply(res,1, function(arg) (arg-mean(arg))/sd(arg))
res <- t(res)


#res<-res[,c("2",'5','4','8','9','0','1','7','3','6')]

#ts<-res[sort(row.names(res)),]
breaklist_1<-seq(-2.5,0,0.2)
breaklist_2<-seq(0,2.6,0.2)
breaklist<-unique(c(breaklist_1,breaklist_2))

myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_ADULT_Kid_normal_Human/PCx/res.0.3/cluster_identity/pheatmap/heatmap/marker_each_cluster.pdf")
pdf(myoutf,width=8,height=17)
pheatmap(mat=res,
         color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(length(breaklist)),
         breaks = breaklist,
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         border_color      = NA,
         show_colnames     = T,
         show_rownames     = T,
         drop_levels       = TRUE,
         scale='none',
         
         fontsize=14)
dev.off()
















PM_markers<-c('GARNL3','CAV1','HGF','LAPTM5','GATA6','CD9','TGFB2','APOC1','F5','HDC','LYZ1',
               'FSCN1','ID1','PLXDC2','SELP','CXCL13','FABP5','ADGRE1','GAS6','ID2',
               'MRC1','FOLR2','LYVE1','MARCO','CD163','CD209','CD74','FPR1',
               'ANXA5','CBR2','CLEC4A1','MAF','ITGB1','PF4',
               'FABP7','CCR2','CSF1R','SELL','ITGAX','HLA-DQA2','HLA-DQA1','HLA-DQB1','HLA-DQP1','ITGAM',
               'RCN3','FCGRT','CCL24','IRF7','EMP1','LGALS1','HAL','RPS18')



PM_markers<-c('VCAN','S100A8','CCR2','CD14','CD1C','MRC1','FABP5','FOLR2','LYVE1','MARCO','CD163',
               'CSF1R','TREM2','HLA-DQB1','TIMD4','KLRB1','NKG7','GZMB','TK1','MKI67','XCR1','CLEC9A','ZEB2','GATA6','F5','SELP','CXCL13','TCF4','IL3RA','CLEC4C','SIGLEC6','CCR7')


Paper_markers<-c('GATA6','SELP','TIMD4','LYVE1','MARCO','MRC1','CCR2','CD1C','CD14','CLEC9A')


myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_ADULT_Kid_normal_Human/PCx/res.0.5/cluster_identity/Dotplot/Dot_plot_paper_genes.pdf"
pdf(myoutf,width=13,height=5)
print(DotPlot(PW_Adult_kit.integrated, cols=c("#08519c","#ef3b2c"),group.by="integrated_snn_res.0.6",features = Paper_markers)+
        theme(axis.text.x = element_text(angle = 45,size=10,hjust = 1))+coord_flip())

dev.off()





myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_ADULT_Kid_normal_Human/PCx/res.0.5/cluster_identity/Dotplot/Dot_plot_MAC_heterogeneity.pdf"
pdf(myoutf,width=15,height=6)
print(DotPlot(PW_Adult_kit.integrated, cols=c("#08519c","#ef3b2c"),group.by="integrated_snn_res.0.6",features = PM_markers)+
        theme(axis.text.x = element_text(angle = 45,size=10,hjust = 1)))

dev.off()


DC_markers<-c('TCF4','IL3RA','CLEC4C','LILRB4','FCER1A','FCER1G','LILRA4',
            'IRF8','IRF4','ID2','BATF3','CLEC9A','XCR1','BTLA',
            'CD163','S100A9','S100A8','VCAN','FCN1','CD14','BST1','CD36',
            'FCGR2B','CD1C','HLA-DQA2','HLA-DQA1','HLA-DQB1','HLA-DQP1',
            'ZEB2','KLF4','ITGAX','ITGAM','CD2','SIRPA','CLEC4A','CLEC10A',
            'AXL','SIGLEC6','CX3CR1','SIGLEC1')


myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_ADULT_Kid_normal_Human/PCx/res.0.5/cluster_identity/Dotplot/Dot_plot_DC_heterogeneity.pdf"
pdf(myoutf,width=15,height=6)
print(DotPlot(PW_Adult_kit.integrated, cols=c("#08519c","#ef3b2c"),group.by="integrated_snn_res.0.6",features = DC_markers)+
        theme(axis.text.x = element_text(angle = 45,size=10,hjust = 1)))

dev.off()





Idents(object = PW_Adult_kit.integrated) <- "integrated_snn_res.0.6"
  tmp <- PW_Adult_kit.integrated
  tmp.markers <- FindAllMarkers(object = tmp, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.3,assay='integrated')


tar<-sort(unique(PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.6))
all_cell<-character()
for(i in 1:(length(tar)-2))
{
  cat("\r",i)
  len<-sum(PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.6==tar[i])
  tag<-sample(1:len,245)
  info<-PW_Adult_kit.integrated@meta.data[PW_Adult_kit.integrated$integrated_snn_res.0.6==tar[i],]
  barcode_use<-row.names(info)[tag]
  all_cell<-c(all_cell,barcode_use)
}
cluster15_16_barcode<-row.names(PW_Adult_kit.integrated@meta.data)[PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.6%in%c('15','16')]

all_cell<-c(all_cell,cluster15_16_barcode)

myoutf1 <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_ADULT_Kid_normal_Human/PCx/res.0.5/cluster_identity/heatmap_marker_genes_ordered_245cells_each_cluster.pdf")

tmp.markers <- tmp.markers[tmp.markers$p_val_adj<0.05,]

ts<-tmp.markers[order(factor(tmp.markers$cluster, levels = tar)),]

top20 <- ts %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

PW_Adult_kit.integrated@meta.data[,"ordered_ident"]<-factor(PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.6,
                                                                    levels=seq(0,16,1))

pdf(myoutf1,12,10)
print(DoHeatmap(PW_Adult_kit.integrated, cells=all_cell,features = top20$gene,group.by="integrated_snn_res.0.6",assay='integrated') + 
        scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')), 
                              mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')), 
                              midpoint = 0, guide = "colourbar", aesthetics = "fill")+NoLegend())
dev.off()






####################Perform GSEA########

screen -r 101052


myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_ADULT_Kid_normal_Human/Seruat_PW_Adult_kit.integrated__finished.RDS"
PW_Adult_kit.integrated<-readRDS(myinf)
#DefaultAssay(PW_Adult_kit.integrated) <- "RNA"
#tar<-"integrated_snn_res.0.6"
# Idents(object = PW_Adult_kit.integrated) <- tar
# tar1<-unique(as.numeric(PW_Adult_kit.integrated@meta.data$"integrated_snn_res.0.6"))
# tar1<-sort(tar1-1)
# for (i in 1:length(tar1))
# {
#  markers <- FindConservedMarkers(PW_Adult_kit.integrated, ident.1 = tar1[i], grouping.var = "orig.ident", verbose = T,logfc.threshold = 0.001,min.cells.group=1,min.pct=0.001)
#  markers[,"cluster"]<-tar1[i]
#  markers[,'gene']<-row.names(markers)
#  tar2<-grep('log2FC',colnames(markers))
#  markers[,'average_Log2FC']<-apply(markers[,colnames(markers)[tar2]],1,mean)
#  markers<-markers[order(markers[,'average_Log2FC'],decreasing=TRUE),]
#  myoutf<-paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_ADULT_Kid_normal_Human/PCx/res.0.5/GSEA/marker/Cluster_",tar1[i],"differentially_expressed_markers.xls")
#  write.table(markers,myoutf,quote=F,sep="\t")
# }




DefaultAssay(PW_Adult_kit.integrated) <- "integrated"
tar<-"integrated_snn_res.0.6"
 Idents(object = PW_Adult_kit.integrated) <- tar
 tar1<-as.character(unique(PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.6))

for(i in 1:length(tar1))
{
  cluster.markers <- FindMarkers(PW_Adult_kit.integrated, ident.1 =tar1[i],logfc.threshold = 0.001,min.pct = 0.001)
  myoutf<-paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_ADULT_Kid_normal_Human/PCx/res.0.5/GSEA/marker/integrated/Cluster_",tar1[i],"differentially_expressed_markers.xls")
  write.table(cluster.markers,myoutf,quote=F,sep="\t")
}


  tmp <- PW_Adult_kit.integrated
  tmp.markers <- FindAllMarkers(object = tmp, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.3,assay='integrated')
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_ADULT_Kid_normal_Human/PCx/cluster_markers_",tar1[i],"_ScRNA.xls")
  write.table(tmp.markers,myoutf,sep="\t",quote=F)
  
  cell.list <- WhichCells(tmp, downsample = 200)
  myoutf1 <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_ADULT_Kid_normal_Human/PCx/heatmap_markers_",tar1[i],"_.pdf")
  pdf(myoutf1,30,35)
  top20 <- tmp.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
  print(DoHeatmap(tmp, cells=cell.list,features = top20$gene,assay='integrated') + NoLegend())
  dev.off()
 




dir <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/PM_KA_manuscript/Human_bulkRNA/analysis/marker/"
files <- list.files(dir)
file_nam <- gsub(".xls","",files)
for(i in 1:length(files))
{
  cat("\r",i)
  tmpinf <- paste0(dir,files[i])
  res <- read.table(tmpinf,sep="\t",quote=NULL)
  res <- res[order(res$Log2FC,decreasing=T),]
  
  xx <- res$Log2FC
  names(xx) <- row.names(res)
  names(xx) <- toupper(names(xx))
  
  tag = duplicated(names(xx))
  res = res[tag==0,]
  row.names(res) <- toupper(row.names(res))
  res[,"Name"] <- row.names(res)
  res[,"metric"] <- res$Log2FC
  res <- res[,c("Name","metric")]
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/PM_KA_manuscript/Human_bulkRNA/analysis/GSEA/",file_nam[i],".rnk")
  write.table(res,file=myoutf,quote=F,sep="\t",row.names=F)
}





dir <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_ADULT_Kid_normal_Human/PCx/res.0.5/GSEA/marker/integrated/"
files <- list.files(dir)
file_nam <- gsub(".xls","",files)
for(i in 1:length(files))
{
  cat("\r",i)
  tmpinf <- paste0(dir,files[i])
  res <- read.table(tmpinf,sep="\t",quote=NULL)
  #tag<-res$average_Log2FC>0
  #res<-res[tag,]
  res <- res[order(res$avg_log2FC,decreasing=T),]
  
  xx <- res$avg_log2FC
  names(xx) <- row.names(res)
  names(xx) <- toupper(names(xx))
  
  tag = duplicated(names(xx))
  res = res[tag==0,]
  row.names(res) <- toupper(row.names(res))
  res[,"Name"] <- row.names(res)
  res[,"metric"] <- res$avg_log2FC
  res <- res[,c("Name","metric")]
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_ADULT_Kid_normal_Human/PCx/res.0.5/GSEA/preranked_GSEA/",file_nam[i],".rnk")
  write.table(res,file=myoutf,quote=F,sep="\t",row.names=F)
}








################heatmap single cell Seurat.
myoutf1 <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.6/cluster_identity/heatmap_marker_genes_different_phenotype_Tcell.pdf")


Idents(PW_Adult_kit.integrated)<-"integrated_snn_res.0.6"
 
tmp.markers <- FindAllMarkers(object = PW_Adult_kit.integrated, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.25,return.thresh=0.05)


top30 <- tmp.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)


pdf(myoutf1,35,35)
print(DoHeatmap(PW_Adult_kit.integrated, features = top50$gene) + 
        scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')), 
                              mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')), 
                              midpoint = 0, guide = "colourbar", aesthetics = "fill"))
dev.off()

tar<-sort(unique(PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.6))
all_cell<-character()
for(i in 1:length(tar))
{
  cat("\r",i)
  len<-sum(PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.3==tar[i])
  tag<-sample(1:len,242)
  info<-PW_Adult_kit.integrated@meta.data[PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.3==tar[i],]
  barcode_use<-row.names(info)[tag]
  all_cell<-combine(all_cell,barcode_use)
}



myoutf1 <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_ADULT_Kid_normal_Human/PCx/res.0.3/cluster_identity/heatmap_marker_genes_ordered_242cells_each_cluster.pdf")


Idents(PW_Adult_kit.integrated)<-"integrated_snn_res.0.3"

tmp.markers <- FindAllMarkers(object = PW_Adult_kit.integrated, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.15,return.thresh=0.05)

ts<-tmp.markers[order(factor(tmp.markers$cluster, levels = tar)),]

top30 <- ts %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)

PW_Adult_kit.integrated@meta.data[,"ordered_ident"]<-factor(PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.3,levels=c('0','1','3','6','7',"2",'5','4','8','9'))

pdf(myoutf1,12,10)
print(DoHeatmap(PW_Adult_kit.integrated, cells=all_cell,features = top30$gene,group.by="ordered_ident") + 
        scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')), 
                              mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')), 
                              midpoint = 0, guide = "colourbar", aesthetics = "fill"))
dev.off()

myoutf2 <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_ADULT_Kid_normal_Human/PCx/res.0.3/cluster_identity/heatmap_marker_top30.xls")

write.table(top30,myoutf2,sep='\t')







tmp.markers <- FindAllMarkers(object = PW_Adult_kit.integrated,only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_ADULT_Kid_normal_Human/PCx/res.0.3/cluster_identity/markers_each_cluster.xls")
write.table(tmp.markers,myoutf,sep="\t",quote=F)

tmp.markers <- FindAllMarkers(object = PW_Adult_kit.integrated,only.pos = T, logfc.threshold = 0.01)
myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_ADULT_Kid_normal_Human/PCx/res.0.3/cluster_identity/markers_each_cluster_0.01.xls")
write.table(tmp.markers,myoutf,sep="\t",quote=F)


myinf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_ADULT_Kid_normal_Human/PCx/res.0.3/cluster_identity/markers_each_cluster_0.01.xls")
markers_all<-read.table(myinf,sep="\t")
marker0<-subset(markers_all,cluster=="0")$gene
marker1<-subset(markers_all,cluster=="1")$gene

###Do venn diagram show overlap of differentially expressed genes between two Trm clusters

myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_ADULT_Kid_normal_Human/PCx/res.0.3/cluster_identity/Venn_Two_Trm_Overlap.pdf")
pdf(myoutf,width=10,height=10)
draw.pairwise.venn(area1 = 109, area2 = 81,  cross.area = 45,category = c("C1-FOS","C2-IFNG"),
                 lwd = rep(2, 2), lty = rep("solid", 2), col =
                   c("#54278f","#fc4e2a"), fill = c("#54278f","#fc4e2a"), alpha = rep(0.5, 2))

dev.off()



myoutf <- paste0("/dartfs-hpc/rc/lab/T/TurkLab/Jichang/bulkTCR/Fred/tumor635.pdf")
pdf(myoutf,width=10,height=10)
draw.pairwise.venn(area1 = 2240, area2 = 2901,  cross.area = 1738,category = c("RNAseq","TCRseq"),
                   lwd = rep(2, 2), lty = rep("solid", 2), col =
                     c("#54278f","#fc4e2a"), fill = c("#54278f","#fc4e2a"), alpha = rep(0.5, 2))

dev.off()














#############Cluster contribution by each different patient##########
Myeloid_PW_Normal_CRC_Kid_B6_BC.integrated<-subset(PW_Normal_CRC_Kid_B6_BC.integrated,subset = CD8A<0.1&CD79A<0.1&CD19<0.1&CD3G<0.1&CD3E<0.1&NCAM1<0.1)
Lymphocyte_PW_Normal_CRC_Kid_B6_BC.integrated<-subset(PW_Normal_CRC_Kid_B6_BC.integrated,subset = CD8A>0.1|CD79A>0.1|CD19>0.1|CD3G>0.1|CD3E>0.1|NCAM1>0.1)

 table(Myeloid_PW_Normal_CRC_Kid_B6_BC.integrated@meta.data$orig.ident)

   ##CF_B6IF1   CF_B6N1  CF_BCIF1   CF_BCN1      Kid1     Kid10      Kid2      Kid4 
   #  4780      7456      5648      7792      2234      2505      1875      1663 
   #  Kid6      Kid7      Kid8      Kid9   PW_CRC1   PW_CRC2 PW_normal Xin_B6_N1 
   #  1959      4368       577      2437       972      3807      6059     14930 
#  Xin_B6_N3 
   # 11436 
table(Lymphocyte_PW_Normal_CRC_Kid_B6_BC.integrated@meta.data$orig.ident)
 #CF_B6IF1   CF_B6N1  CF_BCIF1   CF_BCN1      Kid1     Kid10      Kid2      Kid4 
 #      22       153        39       156      3827      3567      4145      1852 
 #    Kid6      Kid7      Kid8      Kid9   PW_CRC1   PW_CRC2 PW_normal Xin_B6_N1 
 #    1743      2707       338      1478      6742      3733      6113      9935 
#in_B6_N3 
#     5408
 
   Idents(PW_Adult_kit.integrated)<-"integrated_snn_res.0.3"
 
 xx<-table(PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.3,PW_Adult_kit.integrated@meta.data$orig.ident)
 yy<-t(apply(xx,1,function(x){x/sum(x)})   )
 myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/cluster_tissue_distribution/patient_distribution_to_each_cluster.xls"
 write.table(yy,myoutf,sep="\t",quote=F)
 
 data<-matrix(0,ncol(yy)*nrow(yy),3)
 colnames(data)<-c("cluster","patient","value")

 data<-as.data.frame(data)
 data$cluster<-rep(row.names(yy),ncol(yy))

 data$patient<-as.character(sapply(colnames(yy),function(x){rep(x,nrow(yy))}))

 for (i in 1:nrow(data))
 {
  number<-yy[data$cluster[i],data$patient[i]]
  data$value[i]<-number
 }


 data$cluster<-factor(data$cluster,levels=unique(data$cluster))
 data$patient<-factor(data$patient,levels=c("PW_normal","PW_CRC1","PW_CRC2"))

  col<-c("#2171b5","#cb181d","#74c476")

 myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/cluster_patient_distribution/ stackplot_patient_contribution_each_cluster_proportion.pdf"
 pdf(myoutf,width=8,height=8)
 ggplot() + geom_bar(aes(y = value, x = cluster, fill = patient), data = data,
                    stat="identity")+theme(axis.text.x = element_text(angle = 60, hjust = 1,size=10),
                                           axis.title.x = element_blank())+scale_fill_manual(values = col)

 dev.off()


#########Check resolution 0.4


Idents(PW_Adult_kit.integrated)<-"integrated_snn_res.0.4"
 

tmp.markers <- FindMarkers(object = PW_Adult_kit.integrated, ident.1 = 10, ident.2 = 2,only.pos = F, min.pct = 0.15, logfc.threshold = 0.1)
   myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.4/C10_C2.xls")
   write.table(tmp.markers,myoutf,sep="\t",quote=F)
   

tmp.markers <- FindMarkers(object = PW_Adult_kit.integrated, ident.1 = 6, ident.2 = 1,only.pos = F, min.pct = 0.15, logfc.threshold = 0.1)
   myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.4/C6_C1.xls")
   write.table(tmp.markers,myoutf,sep="\t",quote=F)
   


  

########Do heatmap (bulk expression each cluster) show identity of each cluster########

tmp.markers <- FindAllMarkers(object = PW_Adult_kit.integrated, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.25,return.thresh=0.05)

top20 <- tmp.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/Seruat_PW_Adult_kit.integrated__finished.RDS"
exp <- readRDS(myinf)

myoutf1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/expression_data/CD8_normalized_data.xls"
myoutf2 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/expression_data/CD8_normalized_cell_info.xls"


info <- exp@meta.data
res <- GetAssayData(object = exp,assay='integrated',slot='data')
res <- as.matrix(res)
row.names(res) <- toupper(row.names(res))
avg <- apply(res,1, mean)
tag <- avg >0
res <- res[tag,]

res <- unique(res)

tag <- duplicated(row.names(res))

dup_gene <- row.names(res)[tag]

tag <- row.names(res) == dup_gene
#tag==logical(0) Thus ignore the res1,res STBep
#res1 <- res[tag==0,]

#res2 <- res[tag==1,]
#avg <- apply(res2,1, mean)

#res <- rbind(res1, res2[2,])
res <- as.data.frame(res)

write.table(res,myoutf1,sep="\t",quote=F)
write.table(info,myoutf2,sep="\t",quote=F)


myinf1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/expression_data/CD8_normalized_data.xls"
myinf2 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/expression_data/CD8_normalized_cell_info.xls"

data<-read.table(myinf1,sep="\t",quote=NULL)
info<-read.table(myinf2,sep="\t",quote=NULL)

var_gene<-VariableFeatures(object = exp)


gene<-c("CD69","NR4A1","IFNG","TNF","GZMB","GZMA","RUNX3","RGS1","JUN","JUNB","FOS","FOSB","JUND",
        "SELL","CCR7","TCF7","LEF1","KLF2","S1PR1","CX3CR1","KLRG1","NKG7","FGFBP2","PRF1","GNLY",
        "TOX","CXCL13","PDCD1","CTLA4","LAG3","TIGIT","HAVCR2","DUSP1","NR4A2","EOMES","TBX21","CCR5","CXCR6",
        "CCL3","CCL4","CCL5","SLC4A10","RORC","CCR6","IL7R","TRAV1-2","KLRB1",'DPP4',"GZMK","CD27","MALAT1","LPAR6",
        "FXYD7","COTL1","CLU","IL15RA","IL2RA","IL2RG","IL2","FAS","ITGAL","IL32","ID3","BCL6","STAT3",
        "FOXO1","HOPX","ID2","KLF4","RUNX2","FASLG","RUNX1","FOXP1","KLF12","IKZF5","KLF7","GTF3A","HOXB7","ZNF318")





data<-data[row.names(data)%in%var_gene,]

tar <- unique(info$integrated_snn_res.0.3)
#tar <- as.numeric(tar)
tar <- tar[order(tar)]

for(i in 1 : length(tar))
{
  res <- matrix(0, nrow(data),4)
  row.names(res) <- row.names(data)
  colnames(res) <- c("C1_avg_expression","C1_median_expression","C2_avg_expression","C2_median_expression")
  
  tag <- info$integrated_snn_res.0.3 == tar[i]
  
  sam1 <- row.names(info)[tag==1]
  sam2 <- row.names(info)[tag==0]
  
  for(k in 1 : nrow(data))
  {
    cat("\r",i,"----->",k)
    xx1 <- as.numeric(data[k, sam1])
    xx2 <- as.numeric(data[k, sam2])
    
    res[k,"C1_avg_expression"] <- mean(xx1)
    res[k,"C2_avg_expression"] <- mean(xx2)
    res[k,"C1_median_expression"]<-median(xx1)
    res[k,"C2_median_expression"] <- median(xx2)
  }
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/Trm_signature/For_selected_var_genes_cluster_",tar[i],"_vs_others_normalized_expression.xls")
  write.table(res,myoutf,sep="\t",quote=F)
}





data<-data[row.names(data)%in%gene,]

tar <- unique(info$integrated_snn_res.0.3)
#tar <- as.numeric(tar)
tar <- tar[order(tar)]

for(i in 1 : length(tar))
{
  res <- matrix(0, nrow(data),4)
  row.names(res) <- row.names(data)
  colnames(res) <- c("C1_avg_expression","C1_median_expression","C2_avg_expression","C2_median_expression")
  
  tag <- info$integrated_snn_res.0.3 == tar[i]
  
  sam1 <- row.names(info)[tag==1]
  sam2 <- row.names(info)[tag==0]
  
  for(k in 1 : nrow(data))
  {
    cat("\r",i,"----->",k)
    xx1 <- as.numeric(data[k, sam1])
    xx2 <- as.numeric(data[k, sam2])
    
    res[k,"C1_avg_expression"] <- mean(xx1)
    res[k,"C2_avg_expression"] <- mean(xx2)
    res[k,"C1_median_expression"]<-median(xx1)
    res[k,"C2_median_expression"] <- median(xx2)
  }
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/Trm_signature/For_selected_genes_cluster_",tar[i],"_vs_others_normalized_expression.xls")
  write.table(res,myoutf,sep="\t",quote=F)
}


res<-matrix(0,nrow(data),length(tar))
row.names(res)<-row.names(data)
colnames(res)<-tar[order(tar)]

for (i in 1:length(tar))
{
  myinf<-paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/Trm_signature/For_selected_genes_cluster_",tar[i],"_vs_others_normalized_expression.xls")
  xx<-read.table(myinf,sep='\t')
  res[,i]<-xx$C1_avg_expression
}

res <- apply(res,1, function(arg) (arg-mean(arg))/sd(arg))
res <- t(res)


res<-res[,c("2",'5','4','8','9','0','1','7','3','6')]

ts<-res[sort(row.names(res)),]



gene<-c("CD69","NR4A1","IFNG","TNF","GZMB","RUNX3","RGS1","JUN","JUNB","FOS","FOSB","JUND",
        "SELL","CCR7","TCF7","LEF1","KLF2","S1PR1","CX3CR1","KLRG1","NKG7","FGFBP2","PRF1","GNLY",
        "TOX","CXCL13","PDCD1","CTLA4","LAG3","TIGIT","HAVCR2","DUSP1","NR4A2","EOMES","TBX21","CCR5","CXCR6",
        "CCL3","CCL4","CCL5","SLC4A10","RORC","CCR6","IL7R","TRAV1-2","KLRB1",'DPP4',"GZMK","CD27","MALAT1","LPAR6",
        "FXYD7","COTL1","CLU")




Trm_gene<-c("CD69","NR4A1","NR4A2","DUSP1","RGS1")
res_Trm_gene<-res[row.names(res)%in%Trm_gene,]


Cytokine<-c("IFNG","TNF","CCL3","CCL4","CCL5","CXCR6")
res_Cyto<-res[row.names(res)%in%Cytokine,]


JUN<-c("JUN","JUNB","FOS","FOSB","JUND")
res_JUN<-res[row.names(res)%in%JUN,]

Tcir<-c("SELL","CCR7","TCF7","LEF1","KLF2","S1PR1")
res_Tcir<-res[row.names(res)%in%Tcir,]

Teff<-c("CX3CR1","KLRG1","NKG7","PRF1","GZMB","GZMA","IL15RA")
res_Teff<-res[row.names(res)%in%Teff,]

Texh<-c("TOX","CXCL13","PDCD1","CTLA4","LAG3","TIGIT","HAVCR2")
res_Texh<-res[row.names(res)%in%Texh,]

MAIT<-c("SLC4A10","RORC","TRAV1-2","KLRB1",'DPP4',"IL15RA")
res_MAIT<-res[row.names(res)%in%MAIT,]

Tscm<-c("CD27","GZMK","IL7R","MALAT1")
res_Tscm<-res[row.names(res)%in%Tscm,]

Tun<-c("LPAR6","FXYD7","COTL1","CLU")
res_Tun<-res[row.names(res)%in%Tun,]

breaklist_1<-seq(-1,0,0.1)
breaklist_2<-seq(0,2,0.2)
breaklist<-unique(combine(breaklist_1,breaklist_2))



myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/cluster_identity/heatmap_identify_everycluster_all_selected_genes.pdf")
pdf(myoutf,width=20,height=30)
pheatmap(mat=ts,
         color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(length(breaklist)),
         breaks = breaklist,
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         border_color      = NA,
         show_colnames     = T,
         show_rownames     = T,
         drop_levels       = TRUE,
         scale='none',
         
         fontsize=20)
dev.off()



myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/cluster_identity/heatmap_identify_everycluster_MAIT.pdf")
pdf(myoutf,width=15,height=3)
pheatmap(mat=res_MAIT,
         color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(length(breaklist)),
         breaks = breaklist,
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         border_color      = NA,
         show_colnames     = T,
         show_rownames     = T,
         drop_levels       = TRUE,
         scale='none',
         
         fontsize=20)
dev.off()

myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/cluster_identity/heatmap_identify_everycluster_Tscm.pdf")
pdf(myoutf,width=15,height=2.4)
pheatmap(mat=res_Tscm,
         color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(length(breaklist)),
         breaks = breaklist,
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         border_color      = NA,
         show_colnames     = T,
         show_rownames     = T,
         drop_levels       = TRUE,
         scale='none',
         
         fontsize=20)
dev.off()

myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/cluster_identity/heatmap_identify_everycluster_Tun.pdf")
pdf(myoutf,width=15,height=2.4)
pheatmap(mat=res_Tun,
         color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(length(breaklist)),
         breaks = breaklist,
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         border_color      = NA,
         show_colnames     = T,
         show_rownames     = T,
         drop_levels       = TRUE,
         scale='none',
         
         fontsize=20)
dev.off()


myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/cluster_identity/heatmap_identify_everycluster_Trm.pdf")
pdf(myoutf,width=15,height=3)
pheatmap(mat=res_Trm_gene,
         color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(length(breaklist)),
         breaks = breaklist,
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         border_color      = NA,
         show_colnames     = T,
         show_rownames     = T,
         drop_levels       = TRUE,
         scale='none',
         
         fontsize=20)
dev.off()

myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/cluster_identity/heatmap_identify_everycluster_cytokine.pdf")
pdf(myoutf,width=15,height=3.6)
pheatmap(mat=res_Cyto,
         color = colorRampPalette(rev(brewer.pal(n = 9, name =
                                                   "RdYlBu")))(length(breaklist)),
         breaks = breaklist,
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         border_color      = NA,
         show_colnames     = T,
         show_rownames     = T,
         drop_levels       = TRUE,
         scale='none',
         
         fontsize=20)
dev.off()

myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/cluster_identity/heatmap_identify_everycluster_JUN.pdf")
pdf(myoutf,width=15,height=3)
pheatmap(mat=res_JUN,
         color = colorRampPalette(rev(brewer.pal(n = 9, name =
                                                   "RdYlBu")))(length(breaklist)),
         breaks = breaklist,
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         border_color      = NA,
         show_colnames     = T,
         show_rownames     = T,
         drop_levels       = TRUE,
         scale='none',
         
         fontsize=20)
dev.off()

myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/cluster_identity/heatmap_identify_everycluster_Tcir.pdf")
pdf(myoutf,width=15,height=3.6)
pheatmap(mat=res_Tcir,
         color = colorRampPalette(rev(brewer.pal(n = 9, name =
                                                   "RdYlBu")))(length(breaklist)),
         breaks = breaklist,
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         border_color      = NA,
         show_colnames     = T,
         show_rownames     = T,
         drop_levels       = TRUE,
         scale='none',
         
         fontsize=20)
dev.off()

myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/cluster_identity/heatmap_identify_everycluster_Teff.pdf")
pdf(myoutf,width=15,height=3.6)
pheatmap(mat=res_Teff,
         color = colorRampPalette(rev(brewer.pal(n = 9, name =
                                                   "RdYlBu")))(length(breaklist)),
         breaks = breaklist,
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         border_color      = NA,
         show_colnames     = T,
         show_rownames     = T,
         drop_levels       = TRUE,
         scale='none',
         
         fontsize=20)
dev.off()


myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/cluster_identity/heatmap_identify_everycluster_Texh.pdf")
pdf(myoutf,width=15,height=4.2)
pheatmap(mat=res_Texh,
         color = colorRampPalette(rev(brewer.pal(n = 9, name =
                                                   "RdYlBu")))(length(breaklist)),
         breaks = breaklist,
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         border_color      = NA,
         show_colnames     = T,
         show_rownames     = T,
         drop_levels       = TRUE,
         scale='none',
         
         fontsize=20)
dev.off()


































































##################Subset myeloid cells only for downstream analysis#########################
 myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/Seruat_PW_Adult_kit.integrated__finished.RDS"
PW_Adult_kit.integrated<-readRDS(myinf)
 
 tar<-"integrated_snn_res.0.4"
 Idents(object = PW_Adult_kit.integrated) <- tar
 subcluster<-c('0','2','6','4','8','11','13','18','19','14')
 tag<-PW_Adult_kit.integrated@meta.data[,tar]%in%subcluster
 tag1<-GetAssayData(PW_Adult_kit.integrated,assay='RNA', slot='data')['CD3E',]>0.05
 cell_ID<-row.names(PW_Adult_kit.integrated@meta.data)[tag&!tag1]
 
 myinf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/Seruat_PW_PW_Adult_kit_Human__raw.RDS"
 xx<-readRDS(myinf)

 
 Myeloid_PW_Adult_kit.integrated<-subset(x=xx,cells=cell_ID)
table(Myeloid_PW_Adult_kit.integrated@meta.data$orig.ident)
   #Kid1   Kid10   Kid11    Kid2    Kid4    Kid6    Kid7    Kid8    Kid9 PW_Adu1 
   #2213    2622    2347    2021    1716    1675    4445     472    2422    3008 
#  PW_Adu2    PW_C    PW_Q    PW_S    PW_U    PW_W 
  # 2568    5545    2865    3644    1225    4278 

 Myeloid_PW_Adult_kit.integrated.list <- SplitObject(Myeloid_PW_Adult_kit.integrated, split.by = "orig.ident")


 Myeloid_PW_Adult_kit.integrated.list <- lapply(X = Myeloid_PW_Adult_kit.integrated.list, FUN = function(x) {
    x <- NormalizeData(x,assay='RNA')
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 5000,assay='RNA')
 })


############################ Using RPCA to anchor############ 
 features <- SelectIntegrationFeatures(object.list = Myeloid_PW_Adult_kit.integrated.list,nfeatures = 3000)
 Myeloid_PW_Adult_kit.integrated.list <- lapply(X = Myeloid_PW_Adult_kit.integrated.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = T)
    x <- RunPCA(x, features = features, verbose = T)
 })

 immune.anchors <- FindIntegrationAnchors(object.list = Myeloid_PW_Adult_kit.integrated.list, anchor.features = features, reduction = "rpca")
 Myeloid_PW_Adult_kit.integrated <- IntegrateData(anchorset = immune.anchors)
 
 DefaultAssay(Myeloid_PW_Adult_kit.integrated) <- "integrated"



 Myeloid_PW_Adult_kit.integrated <- ScaleData(Myeloid_PW_Adult_kit.integrated, verbose = T)
  Myeloid_PW_Adult_kit.integrated <- RunPCA(Myeloid_PW_Adult_kit.integrated, npcs = 30, verbose = T)


 x<-c(1:15)


 eigs <- Myeloid_PW_Adult_kit.integrated[["pca"]]@stdev^2 #PCA的STBandard dev的平方#
  print(paste0("Variance captured by 35 PCs: ",sum(eigs[x] / sum(eigs)))) #0.840#




# Myeloid_PW_Adult_kit.integrated do clustering and TSNE plot No regresSTBissue ---------------------------------------------


 #[5]Do clustering
 Myeloid_PW_Adult_kit.integrated<-FindNeighbors(Myeloid_PW_Adult_kit.integrated, dims = x)
  Myeloid_PW_Adult_kit.integrated <- FindClusters(Myeloid_PW_Adult_kit.integrated,resolution = seq(0.1,1.0,0.1))

PrintFindClustersParams(object = Myeloid_PW_Adult_kit.integrated)

   clusters_resolution<-sapply(grep("^integrated_snn_res",colnames(Myeloid_PW_Adult_kit.integrated@meta.data),value = TRUE),
                            function(x) length(unique(Myeloid_PW_Adult_kit.integrated@meta.data[,x])))
 clusters_resolution


 #PC=x Myeloid_PW_Adult_kit.integrated
#integrated_snn_res.0.1 integrated_snn_res.0.2 integrated_snn_res.0.3 
#                     7                     12                     14 
#integrated_snn_res.0.4 integrated_snn_res.0.5 integrated_snn_res.0.6 
#                    17                     18                     21 
#integrated_snn_res.0.7 integrated_snn_res.0.8 integrated_snn_res.0.9 
#                    24                     24                     30 
#  integrated_snn_res.1 
#                    33 


 #Set up a loop to look for markers

 tar <- c("integrated_snn_res.0.1","integrated_snn_res.0.2","integrated_snn_res.0.3","integrated_snn_res.0.4",
         "integrated_snn_res.0.5","integrated_snn_res.0.6","integrated_snn_res.0.7","integrated_snn_res.0.8",
         "integrated_snn_res.0.9","integrated_snn_res.1")




 Myeloid_PW_Adult_kit.integrated <- RunUMAP(Myeloid_PW_Adult_kit.integrated, dims= x)

 #Myeloid_PW_Adult_kit.integrated<-RunTSNE(Myeloid_PW_Adult_kit.integrated,reduction = "pca",dims=x)

 for(i in 1:length(tar))
 {
  Idents(object = Myeloid_PW_Adult_kit.integrated) <- tar[i]
  myoutf1 <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/UMAP_SB6_",tar[i],"_undefined.pdf")
  pdf(myoutf1,7,5)
  print(DimPlot(Myeloid_PW_Adult_kit.integrated, reduction = "umap",pt.size=0.2))
  dev.off()
 }


 #for(i in 1:length(tar))
 #{
 # Idents(object = Myeloid_PW_Adult_kit.integrated) <- tar[i]
 # myoutf1 <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/TSNE_SB6_",tar[i],"_undefined.pdf")
 # pdf(myoutf1,7,5)
 # print(DimPlot(Myeloid_PW_Adult_kit.integrated, reduction = "tsne"))
 # dev.off()
 #}


 myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/TSNE_Myeloid_PW_Adult_kit.integrated_groupby_orig.ident.pdf"
  pdf(myoutf,width=8,height=5)
 DimPlot(object = Myeloid_PW_Adult_kit.integrated,group.by='orig.ident',reduction="umap",pt.size=0.2)
 dev.off()

gender_info<-c('F','M','F','M','F','M','M','M','F','F','F','F','M','F','F','M')
age_info<-c('24','7','10','1.3','1.2','14','15','8','17','7','29','49','30','36','36','39')
tar<-unique(Myeloid_PW_Adult_kit.integrated@meta.data$orig.ident)
Myeloid_PW_Adult_kit.integrated@meta.data[,'gender']<-rep(nrow(Myeloid_PW_Adult_kit.integrated@meta.data),'Unknown')
Myeloid_PW_Adult_kit.integrated@meta.data[,'age']<-rep(nrow(Myeloid_PW_Adult_kit.integrated@meta.data),'Unknown')


Myeloid_PW_Adult_kit.integrated@meta.data[,'group']<-rep('Unassigned',nrow(Myeloid_PW_Adult_kit.integrated@meta.data))
tag1<-as.numeric(Myeloid_PW_Adult_kit.integrated@meta.data$age)>18&Myeloid_PW_Adult_kit.integrated@meta.data$gender=='M'
tag2<-as.numeric(Myeloid_PW_Adult_kit.integrated@meta.data$age)>18&Myeloid_PW_Adult_kit.integrated@meta.data$gender=='F'
tag3<-as.numeric(Myeloid_PW_Adult_kit.integrated@meta.data$age)<=18&Myeloid_PW_Adult_kit.integrated@meta.data$gender=='M'
tag4<-as.numeric(Myeloid_PW_Adult_kit.integrated@meta.data$age)<=18&Myeloid_PW_Adult_kit.integrated@meta.data$gender=='F'
Myeloid_PW_Adult_kit.integrated@meta.data$group[tag1]<-'Adult_Male'
Myeloid_PW_Adult_kit.integrated@meta.data$group[tag2]<-'Adult_Female'
Myeloid_PW_Adult_kit.integrated@meta.data$group[tag3]<-'Kid_Male'
Myeloid_PW_Adult_kit.integrated@meta.data$group[tag4]<-'Kid_Female'

for (i in 1:length(tar))
{
  tag<-Myeloid_PW_Adult_kit.integrated@meta.data$orig.ident==tar[i]
  Myeloid_PW_Adult_kit.integrated@meta.data[tag,'gender']<-gender_info[i]
  Myeloid_PW_Adult_kit.integrated@meta.data[tag,'age']<-age_info[i]
}


 myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/Seruat_Myeloid_PW_Adult_kit.integrated__finished.RDS"
 saveRDS(Myeloid_PW_Adult_kit.integrated, file = myoutf)



DefaultAssay(Myeloid_PW_Adult_kit.integrated) <- "integrated"
tar1 <- c("integrated_snn_res.0.3","integrated_snn_res.0.4","integrated_snn_res.0.7","integrated_snn_res.0.6","integrated_snn_res.0.5")

 for(i in 1 : length(tar1))
 {
  cat("\r",i)
  Idents(object = Myeloid_PW_Adult_kit.integrated) <- tar1[i]
  tmp <- Myeloid_PW_Adult_kit.integrated
  tmp.markers <- FindAllMarkers(object = tmp, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.3,assay='integrated')
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/cluster_markers_",tar1[i],"_ScRNA.xls")
  write.table(tmp.markers,myoutf,sep="\t",quote=F)
  
  cell.list <- WhichCells(tmp, downsample = 200)
  myoutf1 <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/heatmap_markers_",tar1[i],"_.pdf")
  pdf(myoutf1,30,35)
  top20 <- tmp.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
  print(DoHeatmap(tmp, cells=cell.list,features = top20$gene,assay='integrated') + NoLegend())
  dev.off()
 }






#########################Find DEG of each cluster using findconserve markers#########
 tar<-"integrated_snn_res.0.3"
 Idents(object = Myeloid_PW_Adult_kit.integrated) <- tar
 tar1<-unique(as.numeric(Myeloid_PW_Adult_kit.integrated@meta.data$"integrated_snn_res.0.3"))
 tar1<-sort(tar1-1)
 DEG<-data.frame()
 for (i in 1:length(tar1))
 {
  markers <- FindConservedMarkers(Myeloid_PW_Adult_kit.integrated, ident.1 = tar1[i], grouping.var = "orig.ident", verbose = T,only.pos = T,logfc.threshold = 0.25,min.cells.group=1)
  markers[,"cluster"]<-tar1[i]
  markers[,'gene']<-row.names(markers)
  tar2<-grep('log2FC',colnames(markers))
  markers[,'average_Log2FC']<-apply(markers[,colnames(markers)[tar2]],1,mean)
  tag<-markers[,'max_pval']<=0.05
  markers<-markers[tag,]
  markers<-markers[order(markers[,'average_Log2FC'],decreasing=TRUE),]
  DEG<-rbind(DEG,markers)
 }





myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/Seruat_Myeloid_PW_Adult_kit.integrated__finished.RDS"
Myeloid_PW_Adult_kit.integrated<-readRDS(myinf)

 DefaultAssay(Myeloid_PW_Adult_kit.integrated) <- "RNA"

data_RNA<-Myeloid_PW_Adult_kit.integrated[['RNA']]@data
data_integrated<-Myeloid_PW_Adult_kit.integrated[['integrated']]@data
data_RNA<-data_RNA[Immune_general,]
data_integrated<-data_integrated[Immune_general,]
cells<-sample(colnames(data_RNA),20)


 DefaultAssay(Myeloid_PW_Adult_kit.integrated) <- "RNA"
Immune_general<-c('CD3G','CD3E','CD8A','CD4','ITGAM','ITGAX','MRC1','CD14','CD68','GATA6','CD1C','LYVE1','CCR2','TIMD4','CD19','TLN1')

p<-FeaturePlot(object = Myeloid_PW_Adult_kit.integrated, features = Immune_general,cols= c('#deebf7', "#bd0026"),
               reduction = "umap",slot = "data",max.cutoff=4,pt.size=0.3,combine=F)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + theme(plot.title = element_text(size=40),
                          axis.title = element_text(family="Helvetica",face = "bold", color = "black",size=20), 
                           axis.text=element_text(family="Helvetica",face = "bold", color = "black",size=16),
                           axis.line = element_line(colour = 'black', size = 1),
                           axis.ticks = element_line(colour = "black", size = 1),
                           axis.ticks.length=unit(.3, "cm"))
}
 
myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/Immune_general_marker.pdf"
pdf(myoutf,width=50,height=30)
print(cowplot::plot_grid(plotlist = p,ncol=5))
dev.off()






tar1 <- c("integrated_snn_res.0.3","integrated_snn_res.0.4")

 for(i in 1 : length(tar1))
 {
  cat("\r",i)
  Idents(object = Myeloid_PW_Adult_kit.integrated) <- tar1[i]
  tmp <- Myeloid_PW_Adult_kit.integrated
  tmp.markers <- FindAllMarkers(object = tmp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/cluster_markers_",tar1[i],"_ScRNA_RNAassay.xls")
  write.table(tmp.markers,myoutf,sep="\t",quote=F)
 }


myeloid_general<-c('SELL','IRF4','F5','CSF1R','VCAN','MERTK','MRC1','CD14','ITGAX','ITGAM','GATA6','CD1C','LYVE1','CCR2','TIMD4','SELP','CX3CR1','CD163',
'FLT3','MARCO','IL3RA','ICAM2','SIGLEC6','RETN','CD226',"CLEC4C","CLEC9A",'CD209','VSIG4','CD36','FCGR1A','MERTK','FCER1A','CD80','CD86','CLEC10A')

flow_panel<-c('FCGR1A','CD14','ITGAM','ITGAX','LYVE1','TIMD4','CD226','SELP','CCR2','HLA-DRA','CCR7','DPP4','MARCO','MRC1','FOLR2',
'LILRB4','CD1C','CLEC9A','SIGLEC6','C5AR1','IL3RA','CD274','FCER1A','SIRPA','CD163')

important<-c('GATA6','SELP','TIMD4','LYVE1','MARCO','CCR2','CD1C','CD14','ICAM2','IL3RA','MRC1','CD163')

paper<-c('GATA6','SELP','TIMD4','LYVE1','MARCO','MRC1','CCR2','CD1C','CD14')


p<-FeaturePlot(object = Myeloid_PW_Adult_kit.integrated, features = paper,cols= c('#deebf7', "#bd0026"),
               reduction = "umap",slot = "data",max.cutoff=4,pt.size=4,combine=F,order = TRUE)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + theme(plot.title = element_text(size=50),
                          axis.title = element_text(family="Helvetica",face = "bold", color = "black",size=40), 
                           axis.text=element_text(family="Helvetica",face = "bold", color = "black",size=30),
                           axis.line = element_line(colour = 'black', size = 3),
                           axis.ticks = element_line(colour = "black", size = 3),
                           axis.ticks.length=unit(1, "cm"))
}
 
myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/paper_marker.pdf"
pdf(myoutf,width=50,height=40)
print(cowplot::plot_grid(plotlist = p,ncol=3))
dev.off()




p<-FeaturePlot(object = Myeloid_PW_Adult_kit.integrated, features = myeloid_general,cols= c('#deebf7', "#bd0026"),
               reduction = "umap",slot = "data",max.cutoff=4,pt.size=0.3,combine=F,order = TRUE)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + theme(plot.title = element_text(size=40),
                          axis.title = element_text(family="Helvetica",face = "bold", color = "black",size=20), 
                           axis.text=element_text(family="Helvetica",face = "bold", color = "black",size=16),
                           axis.line = element_line(colour = 'black', size = 1),
                           axis.ticks = element_line(colour = "black", size = 1),
                           axis.ticks.length=unit(.3, "cm"))
}
 
myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/Myeloid_marker.pdf"
pdf(myoutf,width=80,height=45)
print(cowplot::plot_grid(plotlist = p,ncol=8))
dev.off()




p<-FeaturePlot(object = Myeloid_PW_Adult_kit.integrated, features = flow_panel,cols= c('#deebf7', "#bd0026"),
               reduction = "umap",slot = "data",max.cutoff=4,pt.size=0.3,combine=F,order = TRUE)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + theme(plot.title = element_text(size=40),
                          axis.title = element_text(family="Helvetica",face = "bold", color = "black",size=20), 
                           axis.text=element_text(family="Helvetica",face = "bold", color = "black",size=16),
                           axis.line = element_line(colour = 'black', size = 1),
                           axis.ticks = element_line(colour = "black", size = 1),
                           axis.ticks.length=unit(.3, "cm"))
}
 
myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/Flow_panel_marker.pdf"
pdf(myoutf,width=50,height=50)
print(cowplot::plot_grid(plotlist = p,ncol=5))
dev.off()




p<-FeaturePlot(object = Myeloid_PW_Adult_kit.integrated, features = important,cols= c('#deebf7', "#bd0026"),
               reduction = "umap",slot = "data",max.cutoff=4,pt.size=0.3,combine=F,order = TRUE)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + theme(plot.title = element_text(size=50),
                          axis.title = element_text(family="Helvetica",face = "bold", color = "black",size=30), 
                           axis.text=element_text(family="Helvetica",face = "bold", color = "black",size=20),
                           axis.line = element_line(colour = 'black', size = 1),
                           axis.ticks = element_line(colour = "black", size = 1),
                           axis.ticks.length=unit(.3, "cm"))
}
 
myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/important_marker.pdf"
pdf(myoutf,width=65,height=50)
print(cowplot::plot_grid(plotlist = p,ncol=4))
dev.off()








####Plot adult and children myeloid cell markers seperately seperating sex

tag1<-as.numeric(Myeloid_PW_Adult_kit.integrated@meta.data$age)>18&Myeloid_PW_Adult_kit.integrated@meta.data$gender=='M'
tag2<-as.numeric(Myeloid_PW_Adult_kit.integrated@meta.data$age)>18&Myeloid_PW_Adult_kit.integrated@meta.data$gender=='F'
tag3<-as.numeric(Myeloid_PW_Adult_kit.integrated@meta.data$age)<=18&Myeloid_PW_Adult_kit.integrated@meta.data$gender=='M'
tag4<-as.numeric(Myeloid_PW_Adult_kit.integrated@meta.data$age)<=18&Myeloid_PW_Adult_kit.integrated@meta.data$gender=='F'


p<-FeaturePlot(object = Myeloid_PW_Adult_kit.integrated, cells =row.names(Myeloid_PW_Adult_kit.integrated@meta.data)[tag1], features = myeloid_general,cols= c('#deebf7', "#bd0026"),
               reduction = "umap",slot = "data",max.cutoff=4,pt.size=0.3,combine=F,order = TRUE)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + theme(plot.title = element_text(size=40),
                          axis.title = element_text(family="Helvetica",face = "bold", color = "black",size=20), 
                           axis.text=element_text(family="Helvetica",face = "bold", color = "black",size=16),
                           axis.line = element_line(colour = 'black', size = 1),
                           axis.ticks = element_line(colour = "black", size = 1),
                           axis.ticks.length=unit(.3, "cm"))
}
 
myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/Myeloid_marker_adult_male.pdf"
pdf(myoutf,width=50,height=50)
print(cowplot::plot_grid(plotlist = p,ncol=5))
dev.off()



p<-FeaturePlot(object = Myeloid_PW_Adult_kit.integrated, cells =row.names(Myeloid_PW_Adult_kit.integrated@meta.data)[tag2], features = myeloid_general,cols= c('#deebf7', "#bd0026"),
               reduction = "umap",slot = "data",max.cutoff=4,pt.size=0.3,combine=F,order = TRUE)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + theme(plot.title = element_text(size=40),
                          axis.title = element_text(family="Helvetica",face = "bold", color = "black",size=20), 
                           axis.text=element_text(family="Helvetica",face = "bold", color = "black",size=16),
                           axis.line = element_line(colour = 'black', size = 1),
                           axis.ticks = element_line(colour = "black", size = 1),
                           axis.ticks.length=unit(.3, "cm"))
}
 
myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/Myeloid_marker_adult_female.pdf"
pdf(myoutf,width=50,height=50)
print(cowplot::plot_grid(plotlist = p,ncol=5))
dev.off()




p<-FeaturePlot(object = Myeloid_PW_Adult_kit.integrated, cells =row.names(Myeloid_PW_Adult_kit.integrated@meta.data)[tag3], features = myeloid_general,cols= c('#deebf7', "#bd0026"),
               reduction = "umap",slot = "data",max.cutoff=4,pt.size=0.3,combine=F,order = TRUE)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + theme(plot.title = element_text(size=40),
                          axis.title = element_text(family="Helvetica",face = "bold", color = "black",size=20), 
                           axis.text=element_text(family="Helvetica",face = "bold", color = "black",size=16),
                           axis.line = element_line(colour = 'black', size = 1),
                           axis.ticks = element_line(colour = "black", size = 1),
                           axis.ticks.length=unit(.3, "cm"))
}
 
myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/Myeloid_marker_kid_male.pdf"
pdf(myoutf,width=50,height=50)
print(cowplot::plot_grid(plotlist = p,ncol=5))
dev.off()



p<-FeaturePlot(object = Myeloid_PW_Adult_kit.integrated, cells =row.names(Myeloid_PW_Adult_kit.integrated@meta.data)[tag4], features = myeloid_general,cols= c('#deebf7', "#bd0026"),
               reduction = "umap",slot = "data",max.cutoff=4,pt.size=0.3,combine=F,order = TRUE)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + theme(plot.title = element_text(size=40),
                          axis.title = element_text(family="Helvetica",face = "bold", color = "black",size=20), 
                           axis.text=element_text(family="Helvetica",face = "bold", color = "black",size=16),
                           axis.line = element_line(colour = 'black', size = 1),
                           axis.ticks = element_line(colour = "black", size = 1),
                           axis.ticks.length=unit(.3, "cm"))
}
 
myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/Myeloid_marker_kid_female.pdf"
pdf(myoutf,width=50,height=50)
print(cowplot::plot_grid(plotlist = p,ncol=5))
dev.off()






#############Cluster contribution by each different patient##########
   
   Idents(Myeloid_PW_Adult_kit.integrated)<-"integrated_snn_res.0.5"
 
 xx<-table(Myeloid_PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.3,Myeloid_PW_Adult_kit.integrated@meta.data$orig.ident)
 yy<-t(apply(xx,1,function(x){x/sum(x)})   )
 myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.2/cluster_patient_distribution/patient_distribution_to_each_cluster.xls"
 write.table(yy,myoutf,sep="\t",quote=F)
 
 data<-matrix(0,ncol(yy)*nrow(yy),3)
 colnames(data)<-c("cluster","patient","value")

 data<-as.data.frame(data)
 data$cluster<-rep(row.names(yy),ncol(yy))

 data$patient<-as.character(sapply(colnames(yy),function(x){rep(x,nrow(yy))}))

 for (i in 1:nrow(data))
 {
  number<-yy[data$cluster[i],data$patient[i]]
  data$value[i]<-number
 }


 data$cluster<-factor(data$cluster,levels=unique(data$cluster))
 data$patient<-factor(data$patient,levels=c("PW_normal","PW_CRC1","PW_CRC2",'Kid1','Kid2','Kid4','Kid6','Kid7','Kid8','Kid9','Kid10'))

  col<-c("#2171b5","#cb181d","#74c476")

 myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.3/cluster_patient_distribution/stackplot_patient_contribution_each_cluster_proportion.pdf"
 pdf(myoutf,width=8,height=8)
 ggplot() + geom_bar(aes(y = value, x = cluster, fill = patient), data = data,
                    stat="identity")+theme(axis.text.x = element_text(angle = 60, hjust = 1,size=10),
                                           axis.title.x = element_blank())

 dev.off()






myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/Seruat_Myeloid_PW_Adult_kit.integrated__finished.RDS"
Myeloid_PW_Adult_kit.integrated<-readRDS(myinf)

tar<-unique(Myeloid_PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.5)
tar<-tar[order(tar)]
chi<-as.matrix(table(Myeloid_PW_Adult_kit.integrated@meta.data$orig.ident,Myeloid_PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.5))

chi_square<-chisq.test(chi)

chi_observed<-chi_square$observed
chi_expected<-chi_square$expected
ROE<-chi_observed/chi_expected

myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.5/cluster_patient_distribution/ROE.xls"
write.table(ROE,myoutf,sep="\t",quote=F)



each_cluster_all<-apply(chi,2,sum)
chi_enrichment<-chi
for (i in 1:nrow(chi))
{
    for (k in 1:ncol(chi))
    {
      chi_enrichment[i,k]<-chi[i,k]*sum(chi)/sum(chi[i,])/sum(chi[,k])
    }
}

myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.5/cluster_patient_distribution/hyper_enrichment.xls"
write.table(chi_enrichment,myoutf,sep="\t",quote=F)











myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/Seruat_Myeloid_PW_Adult_kit.integrated__finished.RDS"
Myeloid_PW_Adult_kit.integrated<-readRDS(myinf)

tar<-unique(Myeloid_PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.5)
tar<-tar[order(tar)]
chi<-as.matrix(table(Myeloid_PW_Adult_kit.integrated@meta.data$orig.ident,Myeloid_PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.5))

chi_square<-chisq.test(chi)

chi_observed<-chi_square$observed
chi_expected<-chi_square$expected
ROE<-chi_observed/chi_expected

myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.5/cluster_patient_distribution/ROE.xls"
write.table(ROE,myoutf,sep="\t",quote=F)



each_cluster_all<-apply(chi,2,sum)
chi_enrichment<-chi
for (i in 1:nrow(chi))
{
    for (k in 1:ncol(chi))
    {
      chi_enrichment[i,k]<-chi[i,k]*sum(chi)/sum(chi[i,])/sum(chi[,k])
    }
}

myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.7/cluster_patient_distribution/hyper_enrichment.xls"
write.table(chi_enrichment,myoutf,sep="\t",quote=F)




#################Output cell of each orig.ident and adult versus kid#############
myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/Seruat_Myeloid_PW_Adult_kit.integrated__finished.RDS"
Myeloid_PW_Adult_kit.integrated<-readRDS(myinf)

Idents(Myeloid_PW_Adult_kit.integrated)<-"integrated_snn_res.0.5"
 tar<-unique(Myeloid_PW_Adult_kit.integrated@meta.data$orig.ident)
 for (i in 1:length(tar))
 {
    tag<-Myeloid_PW_Adult_kit.integrated@meta.data$orig.ident==tar[i]
    cell_ID<-row.names(Myeloid_PW_Adult_kit.integrated@meta.data)[tag]
    myoutf1 <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.5/each_sample/UMAP_SB6_",tar[i],"_undefined.pdf")
  pdf(myoutf1,7,5)
  print(DimPlot(Myeloid_PW_Adult_kit.integrated, reduction = "umap",pt.size=0.1,cells =cell_ID))
  dev.off()
 
 }
  
  



myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/Seruat_Myeloid_PW_Adult_kit.integrated__finished.RDS"
Myeloid_PW_Adult_kit.integrated<-readRDS(myinf)

Idents(Myeloid_PW_Adult_kit.integrated)<-"integrated_snn_res.0.5"
cell_ID_adult<-row.names(Myeloid_PW_Adult_kit.integrated@meta.data)[as.numeric(Myeloid_PW_Adult_kit.integrated@meta.data$age)>18]
cell_ID_kid<-row.names(Myeloid_PW_Adult_kit.integrated@meta.data)[as.numeric(Myeloid_PW_Adult_kit.integrated@meta.data$age)<18]


    myoutf1 <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.5/each_sample/All_Adult.pdf")
  pdf(myoutf1,7,5)
  print(DimPlot(Myeloid_PW_Adult_kit.integrated, reduction = "umap",pt.size=0.1,cells =cell_ID_adult))
  dev.off()
 

myoutf1 <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.5/each_sample/All_Kid.pdf")
  pdf(myoutf1,7,5)
  print(DimPlot(Myeloid_PW_Adult_kit.integrated, reduction = "umap",pt.size=0.1,cells =cell_ID_kid))
  dev.off()
 
  
  







################Find differentially expressed genes between adult, kid, male and female even in the same cluster###############################

tar <- unique(Myeloid_PW_Adult_kit.integrated@meta.data$group)
tar1<-combn(tar,2)
tar2<-tar1[,c(2,3,4,5)]

Idents(object = Myeloid_PW_Adult_kit.integrated) <- 'integrated_snn_res.0.5'
 DefaultAssay(Myeloid_PW_Adult_kit.integrated) <- "integrated"




 for(i in 1 : ncol(tar2))
 {
  cat("\r",i)
  tmp <- Myeloid_PW_Adult_kit.integrated
  tmp.markers <- FindMarkers(object = tmp, ident.1=tar2[,i][1],ident.2 =tar2[,i][2],min.pct = 0.1, logfc.threshold = 0.1)
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/DEG_",tar2[,i][1],'_',tar2[,i][2],".xls")
  write.table(tmp.markers,myoutf,sep="\t",quote=F)
 }

tar3<-as.character(unique(Myeloid_PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.5))


for(i in 1 : ncol(tar2))
 {
  cat("\r",i)
  markers_all<-data.frame()
  for (j in 1:length(tar3))
  {
  tmp <- subset(Myeloid_PW_Adult_kit.integrated,idents=tar3[j])
  Idents(object = tmp) <- 'group'
  tmp.markers <- FindMarkers(object = tmp, ident.1=tar2[,i][1],ident.2 =tar2[,i][2],min.pct = 0.1, logfc.threshold = 0.1)
  tmp.markers[,'cluster']<-rep(tar3[j],nrow(tmp.markers))
  markers_all<-rbind(markers_all,tmp.markers)
  }
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/DEG_",tar2[,i][1],'_',tar2[,i][2],".xls")
  write.table(markers_all,myoutf,sep="\t",quote=F)
 }











Myeloid_PW_Adult_kit.integrated@meta.data[,'age_group_only']<-rep('unknown',nrow(Myeloid_PW_Adult_kit.integrated@meta.data))
tag1<-as.numeric(Myeloid_PW_Adult_kit.integrated@meta.data$age)>=18
tag2<-as.numeric(Myeloid_PW_Adult_kit.integrated@meta.data$age)<18
Myeloid_PW_Adult_kit.integrated@meta.data$age_group_only[tag1]<-'Adult'
Myeloid_PW_Adult_kit.integrated@meta.data$age_group_only[tag2]<-'Kid'


  markers_all<-data.frame()
  for (j in 1:length(tar3))
  {
  tmp <- subset(Myeloid_PW_Adult_kit.integrated,idents=tar3[j])
  Idents(object = tmp) <- 'age_group_only'
  tmp.markers <- FindMarkers(object = tmp, ident.1='Adult',ident.2 ='Kid',min.pct = 0.1, logfc.threshold = 0.1)
  tmp.markers[,'cluster']<-rep(tar3[j],nrow(tmp.markers))
  markers_all<-rbind(markers_all,tmp.markers)
  }
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/DEG_Adult_vs_Kid.xls")
  write.table(markers_all,myoutf,sep="\t",quote=F)







################Cell cycle scoring###################
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/Seruat_Myeloid_PW_Adult_kit.integrated__finished.RDS"
Myeloid_PW_Adult_kit.integrated<-readRDS(myinf)

DefaultAssay(Myeloid_PW_Adult_kit.integrated) <- "RNA"
Myeloid_PW_Adult_kit.integrated<-CellCycleScoring(Myeloid_PW_Adult_kit.integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
cell_cycle_score<-t(Myeloid_PW_Adult_kit.integrated@meta.data[,c("S.Score",'G2M.Score')])
adt_assay <- CreateAssayObject(counts = cell_cycle_score)
Myeloid_PW_Adult_kit.integrated[['cell_cycle']]<-adt_assay
DefaultAssay(Myeloid_PW_Adult_kit.integrated) <- "cell_cycle"
p<-FeaturePlot(object = Myeloid_PW_Adult_kit.integrated, features = c('S.Score','G2M.Score'),cols= c('#deebf7', "#bd0026"),
               reduction = "umap",slot = "data",max.cutoff=1.5,pt.size=0.3,combine=F)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + theme(plot.title = element_text(size=40),
                          axis.title = element_text(family="Helvetica",face = "bold", color = "black",size=20), 
                           axis.text=element_text(family="Helvetica",face = "bold", color = "black",size=16),
                           axis.line = element_line(colour = 'black', size = 1),
                           axis.ticks = element_line(colour = "black", size = 1),
                           axis.ticks.length=unit(.3, "cm"))
}
 
myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.5/cell_cycle/Cell_Cycle_Score.pdf"
pdf(myoutf,width=25,height=10)
print(cowplot::plot_grid(plotlist = p,ncol=2))
dev.off()

myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.5/cell_cycle/By_cell_cycle_Myeloid_PW_Adult_kit.integrated.pdf"
  pdf(myoutf,width=8,height=5)
 DimPlot(object = Myeloid_PW_Adult_kit.integrated,group.by='Phase',reduction="umap",pt.size=0.3)
 dev.off()

myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/Seruat_Myeloid_PW_Adult_kit.integrated__finished.RDS"
 saveRDS(Myeloid_PW_Adult_kit.integrated, file = myoutf)

Idents(Myeloid_PW_Adult_kit.integrated)<-"integrated_snn_res.0.5"
 
myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.5/cell_cycle/Cell_cycle_violin.pdf"
  pdf(myoutf,width=13,height=5)
 VlnPlot(object = Myeloid_PW_Adult_kit.integrated,assay='cell_cycle',features = c('S.Score','G2M.Score'),pt.size=0.3)
 dev.off()




Myeloid_marker<-c('XCR1','IL3RA','CLEC4C','SIGLEC6','ITGAM','ITGAX','MRC1','CD14','CD68','GATA6','CD1C','LYVE1','CCR2','TIMD4','CLEC9A')


p<-FeaturePlot(object = Myeloid_PW_Adult_kit.integrated, features = Myeloid_marker,cols= c('#deebf7', "#bd0026"),
               reduction = "umap",slot = "data",max.cutoff=4,pt.size=0.3,combine=F)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + theme(plot.title = element_text(size=40),
                          axis.title = element_text(family="Helvetica",face = "bold", color = "black",size=20), 
                           axis.text=element_text(family="Helvetica",face = "bold", color = "black",size=16),
                           axis.line = element_line(colour = 'black', size = 1),
                           axis.ticks = element_line(colour = "black", size = 1),
                           axis.ticks.length=unit(.3, "cm"))
}
 
myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/Myeloid_marker.pdf"
pdf(myoutf,width=50,height=30)
print(cowplot::plot_grid(plotlist = p,ncol=5))
dev.off()
















##########First get average gene expression of each cluster and then do z transformation to do heatmap for selected genes###########

myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/Seruat_Myeloid_PW_Adult_kit.integrated__finished.RDS"
exp<-readRDS(myinf)
DefaultAssay(exp) <- "RNA"

myoutf1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.5/expression_data/CD8_normalized_data.xls"
myoutf2 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.5/expression_data/CD8_normalized_cell_info.xls"


info <- exp@meta.data
res <- GetAssayData(object = exp,slot='data')
res <- as.matrix(res)
row.names(res) <- toupper(row.names(res))
avg <- apply(res,1, mean)
tag <- avg >0
res <- res[tag,]

res <- unique(res)

tag <- duplicated(row.names(res))

dup_gene <- row.names(res)[tag]

tag <- row.names(res) == dup_gene
#tag==logical(0) Thus ignore the res1,res STBep
#res1 <- res[tag==0,]

#res2 <- res[tag==1,]
#avg <- apply(res2,1, mean)

#res <- rbind(res1, res2[2,])
res <- as.data.frame(res)

write.table(res,myoutf1,sep="\t",quote=F)
write.table(info,myoutf2,sep="\t",quote=F)






myinf1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.5/expression_data/CD8_normalized_data.xls"
myinf2 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.5/expression_data/CD8_normalized_cell_info.xls"

data<-read.table(myinf1,sep="\t",quote=NULL)
info<-read.table(myinf2,sep="\t",quote=NULL)

DC_markers<-c('TCF4','IL3RA','CLEC4C','LILRB4','FCER1A','FCER1G','LILRA4',
            'IRF8','IRF4','ID2','BATF3','CLEC9A','XCR1','BTLA',
            'CD163','S100A9','S100A8','VCAN','FCN1','CD14','BST1','CD36',
            'FCGR2B','CD1C','HLA-DQA2','HLA-DQA1','HLA-DQB1','HLA-DQP1',
            'ZEB2','KLF4','ITGAX','ITGAM','CD2','SIRPA','CLEC4A','CLEC10A',
            'AXL','SIGLEC6','CX3CR1','SIGLEC1')

PM_markers<-c('GARNL3','CAV1','HGF','LAPTM5','GATA6','CD9','TGFB2','APOC1','F5','HDC','LYZ1',
               'FSCN1','ID1','PLXDC2','SELP','CXCL13','FABP5','ADGRE1','GAS6','ID2',
               'MRC1','FOLR2','LYVE1','MARCO','CD163','CD209','CD74','FPR1',
               'ANXA5','CBR2','CLEC4A1','MAF','ITGB1','CXCL13',
               'FABP7','CCR2','CSF1R','SELL','ITGAX','HLA-DQA2','HLA-DQA1','HLA-DQB1','HLA-DQP1','ITGAM',
               'RCN3','FCGRT','CCL24','IRF7','EMP1','LGALS1','HAL','RPS18')

Paper_markers<-c('GATA6','SELP','TIMD4','LYVE1','MARCO','MRC1','CCR2','CD1C','CD14','CLEC9A')

#cluster_markers<-c('S100A9','S100A8','VCAN','FCN1','CD14','MAFB','ITGAM','SPP1',
#                    'FCER1A','CD1C','CD1E','HLA-DPB1','HLA-DQB1','HLA-DRA','CLEC10A','FCGR2B','CLEC4A','TCF4','IL3RA','ITGAX','ZEB2',
#                    'LILRB4','CSF1R','GSN','GPR183','RGS1','IL1B','FOS','NR4A1',
#                    'C1QC','C1QA','C1QB','FN1','LYVE1','FOLR2','MARCO','CD163','FCGR3A','TIMD4','APOE','CD68',
#                    'IL32','CCL2','GZMA','CCL5','NKG7','GNLY','ID3','STAT1','CD36','CCL20','ITGA4','CD3G','CD8A',
#                    'CCL4','CXCL10','CCL3','CXCL9','CXCL11',
#                    'CXCR6','GATA3','RORA','TIGIT','CD69',
#                    'MKI67','TK1','CCNA2','TUBA1B',
#                    'IRF8','ID2','BATF3','CLEC9A','XCR1','BTLA',
#                     'NEAT1','CIITA','NFKBIZ','FOXP1','TLR2','IRF4'
#                    )



data_heatmap_paper<-data[Paper_markers,]
data_heatmap_paper<-na.omit(data_heatmap_paper)

tar <- unique(info$integrated_snn_res.0.5)
#tar <- as.numeric(tar)
tar <- tar[order(tar)]

for(i in 1 : length(tar))
{
  res <- matrix(0, nrow(data_heatmap_paper),4)
  row.names(res) <- row.names(data_heatmap_paper)
  colnames(res) <- c("C1_avg_expression","C1_median_expression","C2_avg_expression","C2_median_expression")
  
  tag <- info$integrated_snn_res.0.5 == tar[i]
  
  sam1 <- row.names(info)[tag==1]
  sam2 <- row.names(info)[tag==0]
  
  for(k in 1 : nrow(data_heatmap_paper))
  {
    cat("\r",i,"----->",k)
    xx1 <- as.numeric(data_heatmap_paper[k, sam1])
    xx2 <- as.numeric(data_heatmap_paper[k, sam2])
    
    res[k,"C1_avg_expression"] <- mean(xx1)
    res[k,"C2_avg_expression"] <- mean(xx2)
    res[k,"C1_median_expression"]<-median(xx1)
    res[k,"C2_median_expression"] <- median(xx2)
  }
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.5/cluster_identity/pheatmap/Paper_Heatmap_",tar[i],"_vs_others_normalized_expression.xls")
  write.table(res,myoutf,sep="\t",quote=F)
}




res<-matrix(0,nrow(data_heatmap_paper),length(tar))
 row.names(res)<-row.names(data_heatmap_paper)
 colnames(res)<-tar[order(tar)]



for (i in 1:length(tar))
{
  myinf<-paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.5/cluster_identity/pheatmap/Paper_Heatmap_",tar[i],"_vs_others_normalized_expression.xls")
  xx<-read.table(myinf,sep='\t')
  res[,i]<-xx$C1_avg_expression
}

res <- apply(res,1, function(arg) (arg-mean(arg))/sd(arg))
res <- t(res)


#res<-res[,c("2",'5','4','8','9','0','1','7','3','6')]

#ts<-res[sort(row.names(res)),]
breaklist_1<-seq(-2.5,0,0.2)
breaklist_2<-seq(0,2.6,0.2)
breaklist<-unique(c(breaklist_1,breaklist_2))

myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.5/cluster_identity/pheatmap/heatmap/Paper_Heatmap_gene_each_cluster.pdf")
pdf(myoutf,width=20,height=5)
pheatmap(mat=res,
         color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(length(breaklist)),
         breaks = breaklist,
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         border_color      = NA,
         show_colnames     = T,
         show_rownames     = T,
         drop_levels       = TRUE,
         scale='none',
         
         fontsize=20)
dev.off()








data_DC<-data[DC_markers,]
data_DC<-na.omit(data_DC)

tar <- unique(info$integrated_snn_res.0.5)
#tar <- as.numeric(tar)
tar <- tar[order(tar)]

for(i in 1 : length(tar))
{
  res <- matrix(0, nrow(data_DC),4)
  row.names(res) <- row.names(data_DC)
  colnames(res) <- c("C1_avg_expression","C1_median_expression","C2_avg_expression","C2_median_expression")
  
  tag <- info$integrated_snn_res.0.5 == tar[i]
  
  sam1 <- row.names(info)[tag==1]
  sam2 <- row.names(info)[tag==0]
  
  for(k in 1 : nrow(data_DC))
  {
    cat("\r",i,"----->",k)
    xx1 <- as.numeric(data_DC[k, sam1])
    xx2 <- as.numeric(data_DC[k, sam2])
    
    res[k,"C1_avg_expression"] <- mean(xx1)
    res[k,"C2_avg_expression"] <- mean(xx2)
    res[k,"C1_median_expression"]<-median(xx1)
    res[k,"C2_median_expression"] <- median(xx2)
  }
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.5/cluster_identity/pheatmap/Human_DC_gene_",tar[i],"_vs_others_normalized_expression.xls")
  write.table(res,myoutf,sep="\t",quote=F)
}




res<-matrix(0,nrow(data_DC),length(tar))
 row.names(res)<-row.names(data_DC)
 colnames(res)<-tar[order(tar)]



for (i in 1:length(tar))
{
  myinf<-paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.5/cluster_identity/pheatmap/Human_DC_gene_",tar[i],"_vs_others_normalized_expression.xls")
  xx<-read.table(myinf,sep='\t')
  res[,i]<-xx$C1_avg_expression
}

res <- apply(res,1, function(arg) (arg-mean(arg))/sd(arg))
res <- t(res)


#res<-res[,c("2",'5','4','8','9','0','1','7','3','6')]

#ts<-res[sort(row.names(res)),]
breaklist_1<-seq(-2.5,0,0.2)
breaklist_2<-seq(0,2.6,0.2)
breaklist<-unique(c(breaklist_1,breaklist_2))

myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.5/cluster_identity/pheatmap/heatmap/Human_DC_gene_each_cluster.pdf")
pdf(myoutf,width=25,height=15)
pheatmap(mat=res,
         color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(length(breaklist)),
         breaks = breaklist,
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         border_color      = NA,
         show_colnames     = T,
         show_rownames     = T,
         drop_levels       = TRUE,
         scale='none',
         
         fontsize=20)
dev.off()












data_PM<-data[PM_markers,]
data_PM<-na.omit(data_PM)

tar <- unique(info$integrated_snn_res.0.5)
#tar <- as.numeric(tar)
tar <- tar[order(tar)]

for(i in 1 : length(tar))
{
  res <- matrix(0, nrow(data_PM),4)
  row.names(res) <- row.names(data_PM)
  colnames(res) <- c("C1_avg_expression","C1_median_expression","C2_avg_expression","C2_median_expression")
  
  tag <- info$integrated_snn_res.0.5 == tar[i]
  
  sam1 <- row.names(info)[tag==1]
  sam2 <- row.names(info)[tag==0]
  
  for(k in 1 : nrow(data_PM))
  {
    cat("\r",i,"----->",k)
    xx1 <- as.numeric(data_PM[k, sam1])
    xx2 <- as.numeric(data_PM[k, sam2])
    
    res[k,"C1_avg_expression"] <- mean(xx1)
    res[k,"C2_avg_expression"] <- mean(xx2)
    res[k,"C1_median_expression"]<-median(xx1)
    res[k,"C2_median_expression"] <- median(xx2)
  }
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.5/cluster_identity/pheatmap/Mouse_Peritoneal_MAC_gene_",tar[i],"_vs_others_normalized_expression.xls")
  write.table(res,myoutf,sep="\t",quote=F)
}




res<-matrix(0,nrow(data_PM),length(tar))
 row.names(res)<-row.names(data_PM)
 colnames(res)<-tar[order(tar)]



for (i in 1:length(tar))
{
  myinf<-paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.5/cluster_identity/pheatmap/Mouse_Peritoneal_MAC_gene_",tar[i],"_vs_others_normalized_expression.xls")
  xx<-read.table(myinf,sep='\t')
  res[,i]<-xx$C1_avg_expression
}

res <- apply(res,1, function(arg) (arg-mean(arg))/sd(arg))
res <- t(res)


#res<-res[,c("2",'5','4','8','9','0','1','7','3','6')]

#ts<-res[sort(row.names(res)),]
breaklist_1<-seq(-2.5,0,0.2)
breaklist_2<-seq(0,2.6,0.2)
breaklist<-unique(c(breaklist_1,breaklist_2))

myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.5/cluster_identity/pheatmap/heatmap/Mouse_Peritoneal_MAC_gene_each_cluster.pdf")
pdf(myoutf,width=15,height=15)
pheatmap(mat=res,
         color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(length(breaklist)),
         breaks = breaklist,
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         border_color      = NA,
         show_colnames     = T,
         show_rownames     = T,
         drop_levels       = TRUE,
         scale='none',
         
         fontsize=20)
dev.off()










data_each_cluster<-data[cluster_markers,]
data_each_cluster<-na.omit(data_each_cluster)

tar <- unique(info$integrated_snn_res.0.3)
#tar <- as.numeric(tar)
tar <- tar[order(tar)]

for(i in 1 : length(tar))
{
  res <- matrix(0, nrow(data_each_cluster),4)
  row.names(res) <- row.names(data_each_cluster)
  colnames(res) <- c("C1_avg_expression","C1_median_expression","C2_avg_expression","C2_median_expression")
  
  tag <- info$integrated_snn_res.0.3 == tar[i]
  
  sam1 <- row.names(info)[tag==1]
  sam2 <- row.names(info)[tag==0]
  
  for(k in 1 : nrow(data_each_cluster))
  {
    cat("\r",i,"----->",k)
    xx1 <- as.numeric(data_each_cluster[k, sam1])
    xx2 <- as.numeric(data_each_cluster[k, sam2])
    
    res[k,"C1_avg_expression"] <- mean(xx1)
    res[k,"C2_avg_expression"] <- mean(xx2)
    res[k,"C1_median_expression"]<-median(xx1)
    res[k,"C2_median_expression"] <- median(xx2)
  }
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.3/cluster_identity/pheatmap/marker_each_cluster_",tar[i],"_vs_others_normalized_expression.xls")
  write.table(res,myoutf,sep="\t",quote=F)
}




res<-matrix(0,nrow(data_each_cluster),length(tar))
 row.names(res)<-row.names(data_each_cluster)
 colnames(res)<-tar[order(tar)]



for (i in 1:length(tar))
{
  myinf<-paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.3/cluster_identity/pheatmap/marker_each_cluster_",tar[i],"_vs_others_normalized_expression.xls")
  xx<-read.table(myinf,sep='\t')
  res[,i]<-xx$C1_avg_expression
}

res <- apply(res,1, function(arg) (arg-mean(arg))/sd(arg))
res <- t(res)


#res<-res[,c("2",'5','4','8','9','0','1','7','3','6')]

#ts<-res[sort(row.names(res)),]
breaklist_1<-seq(-2.5,0,0.2)
breaklist_2<-seq(0,2.6,0.2)
breaklist<-unique(c(breaklist_1,breaklist_2))

myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.3/cluster_identity/pheatmap/heatmap/marker_each_cluster.pdf")
pdf(myoutf,width=8,height=17)
pheatmap(mat=res,
         color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(length(breaklist)),
         breaks = breaklist,
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         border_color      = NA,
         show_colnames     = T,
         show_rownames     = T,
         drop_levels       = TRUE,
         scale='none',
         
         fontsize=14)
dev.off()
















PM_markers<-c('GARNL3','CAV1','HGF','LAPTM5','GATA6','CD9','TGFB2','APOC1','F5','HDC','LYZ1',
               'FSCN1','ID1','PLXDC2','SELP','CXCL13','FABP5','ADGRE1','GAS6','ID2',
               'MRC1','FOLR2','LYVE1','MARCO','CD163','CD209','CD74','FPR1',
               'ANXA5','CBR2','CLEC4A1','MAF','ITGB1','PF4',
               'FABP7','CCR2','CSF1R','SELL','ITGAX','HLA-DQA2','HLA-DQA1','HLA-DQB1','HLA-DQP1','ITGAM',
               'RCN3','FCGRT','CCL24','IRF7','EMP1','LGALS1','HAL','RPS18')



PM_markers<-c('VCAN','S100A8','CCR2','CD14','CD1C','MRC1','FABP5','FOLR2','LYVE1','MARCO','CD163',
               'CSF1R','TREM2','HLA-DQB1','TIMD4','KLRB1','NKG7','GZMB','TK1','MKI67','XCR1','CLEC9A','ZEB2','GATA6','F5','SELP','CXCL13','TCF4','IL3RA','CLEC4C','SIGLEC6','CCR7')


Paper_markers<-c('GATA6','SELP','TIMD4','LYVE1','MARCO','MRC1','CCR2','CD1C','CD14','CLEC9A')


myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.5/cluster_identity/Dotplot/Dot_plot_paper_genes.pdf"
pdf(myoutf,width=13,height=5)
print(DotPlot(Myeloid_PW_Adult_kit.integrated, cols=c("#08519c","#ef3b2c"),group.by="integrated_snn_res.0.5",features = Paper_markers)+
        theme(axis.text.x = element_text(angle = 45,size=10,hjust = 1))+coord_flip())

dev.off()





myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.5/cluster_identity/Dotplot/Dot_plot_MAC_heterogeneity.pdf"
pdf(myoutf,width=15,height=6)
print(DotPlot(Myeloid_PW_Adult_kit.integrated, cols=c("#08519c","#ef3b2c"),group.by="integrated_snn_res.0.5",features = PM_markers)+
        theme(axis.text.x = element_text(angle = 45,size=10,hjust = 1)))

dev.off()


DC_markers<-c('TCF4','IL3RA','CLEC4C','LILRB4','FCER1A','FCER1G','LILRA4',
            'IRF8','IRF4','ID2','BATF3','CLEC9A','XCR1','BTLA',
            'CD163','S100A9','S100A8','VCAN','FCN1','CD14','BST1','CD36',
            'FCGR2B','CD1C','HLA-DQA2','HLA-DQA1','HLA-DQB1','HLA-DQP1',
            'ZEB2','KLF4','ITGAX','ITGAM','CD2','SIRPA','CLEC4A','CLEC10A',
            'AXL','SIGLEC6','CX3CR1','SIGLEC1')


myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.5/cluster_identity/Dotplot/Dot_plot_DC_heterogeneity.pdf"
pdf(myoutf,width=15,height=6)
print(DotPlot(Myeloid_PW_Adult_kit.integrated, cols=c("#08519c","#ef3b2c"),group.by="integrated_snn_res.0.5",features = DC_markers)+
        theme(axis.text.x = element_text(angle = 45,size=10,hjust = 1)))

dev.off()





Idents(object = Myeloid_PW_Adult_kit.integrated) <- "integrated_snn_res.0.5"
  tmp <- Myeloid_PW_Adult_kit.integrated
  tmp.markers <- FindAllMarkers(object = tmp, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.3,assay='integrated')


tar<-sort(unique(Myeloid_PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.5))
all_cell<-character()
for(i in 1:(length(tar)-2))
{
  cat("\r",i)
  len<-sum(Myeloid_PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.5==tar[i])
  tag<-sample(1:len,245)
  info<-Myeloid_PW_Adult_kit.integrated@meta.data[Myeloid_PW_Adult_kit.integrated$integrated_snn_res.0.5==tar[i],]
  barcode_use<-row.names(info)[tag]
  all_cell<-c(all_cell,barcode_use)
}
cluster15_16_barcode<-row.names(Myeloid_PW_Adult_kit.integrated@meta.data)[Myeloid_PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.5%in%c('15','16')]

all_cell<-c(all_cell,cluster15_16_barcode)

myoutf1 <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.5/cluster_identity/heatmap_marker_genes_ordered_245cells_each_cluster.pdf")

tmp.markers <- tmp.markers[tmp.markers$p_val_adj<0.05,]

ts<-tmp.markers[order(factor(tmp.markers$cluster, levels = tar)),]

top20 <- ts %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

Myeloid_PW_Adult_kit.integrated@meta.data[,"ordered_ident"]<-factor(Myeloid_PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.5,
                                                                    levels=seq(0,16,1))

pdf(myoutf1,12,10)
print(DoHeatmap(Myeloid_PW_Adult_kit.integrated, cells=all_cell,features = top20$gene,group.by="integrated_snn_res.0.5",assay='integrated') + 
        scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')), 
                              mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')), 
                              midpoint = 0, guide = "colourbar", aesthetics = "fill")+NoLegend())
dev.off()






##########Check GATA6################
GATA6_exp<-as.matrix(GetAssayData(Myeloid_PW_Adult_kit.integrated,assay='RNA',slot='data')[c('GATA6','SELP'),])
apply(GATA6_exp,1,summary)

#              GATA6        SELP
#Min.    0.000000000 0.000000000
#1st Qu. 0.000000000 0.000000000
#Median  0.000000000 0.000000000
#Mean    0.009995213 0.007311954
#3rd Qu. 0.000000000 0.000000000
#Max.    2.931832980 2.733619776

tag1<-GATA6_exp[1,]>mean(GATA6_exp[1,])
tag2<-GATA6_exp[2,]>mean(GATA6_exp[2,])

GATA6_pos<-colnames(GATA6_exp)[tag1]
SELP_pos<-colnames(GATA6_exp)[tag2]

Myeloid_PW_Adult_kit.integrated@meta.data[,'GATA6']<-rep("Neg",nrow(Myeloid_PW_Adult_kit.integrated@meta.data))
Myeloid_PW_Adult_kit.integrated@meta.data[,'SELP']<-rep("Neg",nrow(Myeloid_PW_Adult_kit.integrated@meta.data))
Myeloid_PW_Adult_kit.integrated@meta.data[GATA6_pos,'GATA6']<-'Pos'
Myeloid_PW_Adult_kit.integrated@meta.data[SELP_pos,'SELP']<-'Pos'

GATA6_dis<-table(Myeloid_PW_Adult_kit.integrated@meta.data[,'GATA6'],Myeloid_PW_Adult_kit.integrated@meta.data[,'orig.ident'])
SELP_dis<-table(Myeloid_PW_Adult_kit.integrated@meta.data[,'SELP'],Myeloid_PW_Adult_kit.integrated@meta.data[,'orig.ident'])



chi_square<-chisq.test(GATA6_dis)

chi_observed<-chi_square$observed
chi_expected<-chi_square$expected
ROE_GATA6<-chi_observed/chi_expected


chi_square<-chisq.test(SELP_dis)

chi_observed<-chi_square$observed
chi_expected<-chi_square$expected
ROE_SELP<-chi_observed/chi_expected



GATA6_enrichment<-GATA6_dis
for (i in 1:nrow(GATA6_dis))
{
    for (k in 1:ncol(GATA6_dis))
    {
      GATA6_enrichment[i,k]<-GATA6_dis[i,k]*sum(GATA6_dis)/sum(GATA6_dis[i,])/sum(GATA6_dis[,k])
    }
}



myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.7/cluster_patient_distribution/GATA6_enrichment.xls"
write.table(ROE_GATA6,myoutf,sep="\t",quote=F)

myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.7/cluster_patient_distribution/SELP_enrichment.xls"
write.table(ROE_SELP,myoutf,sep="\t",quote=F)





##############Find doublets in the dataset#############
 
 myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/Seruat_Myeloid_PW_Adult_kit.integrated__finished.RDS"
 Myeloid_PW_Adult_kit.integrated<-readRDS(myinf)
 cell_ID<-row.names(Myeloid_PW_Adult_kit.integrated@meta.data)

 myinf1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/raw/Nan_PW_Normal"
 PW_normal.data <- Read10X(data.dir=myinf1)
 PW_normal <- CreateSeuratObject(counts = PW_normal.data, project = "PW_normal",min.cells=3,min.features=200)
 PW_normal<-RenameCells(object=PW_normal,add.cell.id=unique(PW_normal@meta.data$orig.ident))
 PW_normal<-RenameCells(object=PW_normal,new.names=gsub('-1','',Cells(PW_normal)))


  myinf2 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/raw/Nan_PW_CRC1/sample_feature_bc_matrix"
 PW_CRC1.data <- Read10X(data.dir=myinf2)
 PW_CRC1 <- CreateSeuratObject(counts = PW_CRC1.data, project = "PW_CRC1",min.cells=3,min.features=200)
 PW_CRC1<-RenameCells(object=PW_CRC1,add.cell.id=unique(PW_CRC1@meta.data$orig.ident))
 PW_CRC1<-RenameCells(object=PW_CRC1,new.names=gsub('-1','',Cells(PW_CRC1)))


 myinf3 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/raw/Nan_PW_CRC2/sample_feature_bc_matrix"
 PW_CRC2.data <- Read10X(data.dir=myinf3)
 PW_CRC2 <- CreateSeuratObject(counts = PW_CRC2.data, project = "PW_CRC2",min.cells=3,min.features=200)
 PW_CRC2<-RenameCells(object=PW_CRC2,add.cell.id=unique(PW_CRC2@meta.data$orig.ident))
 PW_CRC2<-RenameCells(object=PW_CRC2,new.names=gsub('-1','',Cells(PW_CRC2)))

 xx<- merge(PW_normal, y = c(PW_CRC1, PW_CRC2), project = "PW_Adult_kit")


 Myeloid_PW_Adult_kit.doublet<-subset(x=xx,cells=cell_ID)
 
 Myeloid_PW_Adult_kit.doublet.list <- SplitObject(Myeloid_PW_Adult_kit.doublet, split.by = "orig.ident")

 Myeloid_PW_Adult_kit.doublet.list <- lapply(X = Myeloid_PW_Adult_kit.doublet.list, FUN = function(x) {
    x <- NormalizeData(x,assay='RNA')
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000,assay='RNA')
    x <- ScaleData(x, verbose = T)
    x <- RunPCA(x,verbose = T)
    x<-FindNeighbors(x, dims = 1:15)
   x <- FindClusters(x,resolution = 0.3)
 
 })









 Myeloid_PW_Adult_kit.doublet.list <- lapply(X = Myeloid_PW_Adult_kit.doublet.list, FUN = function(x) {
    yy <- paramSweep_v3(x, PCs = 1:15, sct = FALSE)
    zz<-summarizeSweep(yy, GT = FALSE)
    tt<-find.pK(zz)
    pK<-tt %>% filter(BCmetric == max(BCmetric)) %>% select (pK)
    pK<- as.numeric (as.character(pK[[1]]))
    homotypic.prop <- modelHomotypic(x@meta.data$"integrated_snn_res.0.3")
    nExp_poi <- round(0.055*nrow(x@meta.data)) 
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    x <- doubletFinder_v3(x, PCs = 1:15, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
 })


 Myeloid_PW_Adult_kit.doublet<-merge(
  Myeloid_PW_Adult_kit.doublet.list[[1]],
  y =c(Myeloid_PW_Adult_kit.doublet.list[[2]], Myeloid_PW_Adult_kit.doublet.list[[3]]),
  add.cell.ids = NULL,
  project = "Myeloid_PW_Adult_kit.doublet"
)

tag<-Myeloid_PW_Adult_kit.doublet@meta.data[,"DF.classifications_0.25_0.17_448"]=='Doublet'
doublet_ID<-row.names(Myeloid_PW_Adult_kit.doublet@meta.data[tag,])

Myeloid_PW_Adult_kit.integrated@meta.data[,"doublet"]<-rep('Singlet',nrow(Myeloid_PW_Adult_kit.integrated@meta.data))
Myeloid_PW_Adult_kit.integrated@meta.data[doublet_ID,"doublet"]<-"Doublet"

myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/doublet_finder.pdf"
  pdf(myoutf,width=8,height=5)
 DimPlot(object = Myeloid_PW_Adult_kit.integrated,group.by='doublet',reduction="umap",pt.size=0.3,cells.highlight=doublet_ID,cols.highlight = "#e6550d",
  sizes.highlight = 0.3)
 dev.off()














####################Perform GSEA########

screen -r 101052


myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/Seruat_Myeloid_PW_Adult_kit.integrated__finished.RDS"
Myeloid_PW_Adult_kit.integrated<-readRDS(myinf)
#DefaultAssay(Myeloid_PW_Adult_kit.integrated) <- "RNA"
#tar<-"integrated_snn_res.0.5"
# Idents(object = Myeloid_PW_Adult_kit.integrated) <- tar
# tar1<-unique(as.numeric(Myeloid_PW_Adult_kit.integrated@meta.data$"integrated_snn_res.0.5"))
# tar1<-sort(tar1-1)
# for (i in 1:length(tar1))
# {
#  markers <- FindConservedMarkers(Myeloid_PW_Adult_kit.integrated, ident.1 = tar1[i], grouping.var = "orig.ident", verbose = T,logfc.threshold = 0.001,min.cells.group=1,min.pct=0.001)
#  markers[,"cluster"]<-tar1[i]
#  markers[,'gene']<-row.names(markers)
#  tar2<-grep('log2FC',colnames(markers))
#  markers[,'average_Log2FC']<-apply(markers[,colnames(markers)[tar2]],1,mean)
#  markers<-markers[order(markers[,'average_Log2FC'],decreasing=TRUE),]
#  myoutf<-paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.5/GSEA/marker/Cluster_",tar1[i],"differentially_expressed_markers.xls")
#  write.table(markers,myoutf,quote=F,sep="\t")
# }




DefaultAssay(Myeloid_PW_Adult_kit.integrated) <- "integrated"
tar<-"integrated_snn_res.0.5"
 Idents(object = Myeloid_PW_Adult_kit.integrated) <- tar
 tar1<-as.character(unique(Myeloid_PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.5))

for(i in 1:length(tar1))
{
  cluster.markers <- FindMarkers(Myeloid_PW_Adult_kit.integrated, ident.1 =tar1[i],logfc.threshold = 0.001,min.pct = 0.001)
  myoutf<-paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.5/GSEA/marker/integrated/Cluster_",tar1[i],"differentially_expressed_markers.xls")
  write.table(cluster.markers,myoutf,quote=F,sep="\t")
}


  tmp <- Myeloid_PW_Adult_kit.integrated
  tmp.markers <- FindAllMarkers(object = tmp, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.3,assay='integrated')
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/cluster_markers_",tar1[i],"_ScRNA.xls")
  write.table(tmp.markers,myoutf,sep="\t",quote=F)
  
  cell.list <- WhichCells(tmp, downsample = 200)
  myoutf1 <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/heatmap_markers_",tar1[i],"_.pdf")
  pdf(myoutf1,30,35)
  top20 <- tmp.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
  print(DoHeatmap(tmp, cells=cell.list,features = top20$gene,assay='integrated') + NoLegend())
  dev.off()
 




dir <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/PM_KA_manuscript/Human_bulkRNA/analysis/marker/"
files <- list.files(dir)
file_nam <- gsub(".xls","",files)
for(i in 1:length(files))
{
  cat("\r",i)
  tmpinf <- paste0(dir,files[i])
  res <- read.table(tmpinf,sep="\t",quote=NULL)
  res <- res[order(res$Log2FC,decreasing=T),]
  
  xx <- res$Log2FC
  names(xx) <- row.names(res)
  names(xx) <- toupper(names(xx))
  
  tag = duplicated(names(xx))
  res = res[tag==0,]
  row.names(res) <- toupper(row.names(res))
  res[,"Name"] <- row.names(res)
  res[,"metric"] <- res$Log2FC
  res <- res[,c("Name","metric")]
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/PM_KA_manuscript/Human_bulkRNA/analysis/GSEA/",file_nam[i],".rnk")
  write.table(res,file=myoutf,quote=F,sep="\t",row.names=F)
}





dir <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.5/GSEA/marker/integrated/"
files <- list.files(dir)
file_nam <- gsub(".xls","",files)
for(i in 1:length(files))
{
  cat("\r",i)
  tmpinf <- paste0(dir,files[i])
  res <- read.table(tmpinf,sep="\t",quote=NULL)
  #tag<-res$average_Log2FC>0
  #res<-res[tag,]
  res <- res[order(res$avg_log2FC,decreasing=T),]
  
  xx <- res$avg_log2FC
  names(xx) <- row.names(res)
  names(xx) <- toupper(names(xx))
  
  tag = duplicated(names(xx))
  res = res[tag==0,]
  row.names(res) <- toupper(row.names(res))
  res[,"Name"] <- row.names(res)
  res[,"metric"] <- res$avg_log2FC
  res <- res[,c("Name","metric")]
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.5/GSEA/preranked_GSEA/",file_nam[i],".rnk")
  write.table(res,file=myoutf,quote=F,sep="\t",row.names=F)
}























##########First get average gene expression of each cluster and then do z transformation to do heatmap for selected genes###########

myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/Seruat_Myeloid_PW_Adult_kit.integrated__finished.RDS"
exp <- readRDS(myinf)
DefaultAssay(exp) <- "integrated"

myoutf1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.3/expression_data/CD8_normalized_data.xls"
myoutf2 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.3/expression_data/CD8_normalized_cell_info.xls"


info <- exp@meta.data
res <- GetAssayData(object = exp)
res <- as.matrix(res)
row.names(res) <- toupper(row.names(res))
avg <- apply(res,1, mean)
tag <- avg >0
res <- res[tag,]

res <- unique(res)

tag <- duplicated(row.names(res))

dup_gene <- row.names(res)[tag]

tag <- row.names(res) == dup_gene
#tag==logical(0) Thus ignore the res1,res STBep
#res1 <- res[tag==0,]

#res2 <- res[tag==1,]
#avg <- apply(res2,1, mean)

#res <- rbind(res1, res2[2,])
res <- as.data.frame(res)

write.table(res,myoutf1,sep="\t",quote=F)
write.table(info,myoutf2,sep="\t",quote=F)

myinf1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.3/expression_data/CD8_normalized_data.xls"
myinf2 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.3/expression_data/CD8_normalized_cell_info.xls"

data<-read.table(myinf1,sep="\t",quote=NULL)
info<-read.table(myinf2,sep="\t",quote=NULL)


dir <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/PM_KA_manuscript/Human_bulkRNA/analysis/marker/"
files <- list.files(dir)
file_nam <- gsub(".txt","",files)
for(i in 1:length(files))
{
  cat("\r",i)
  tmpinf <- paste0(dir,files[i])
  res1 <- read.table(tmpinf,sep="\t",quote=NULL)
  #tag<-res$average_Log2FC>0
  #res<-res[tag,]
  res1 <- res1[order(res1$Log2FC,decreasing=T),]
  FC2_sig<-row.names(res1[res1$Log2FC>2,])
  FC1_sig<-row.names(res1[res1$Log2FC>1,])
  top100_sig<-row.names(res1[1:100,])
  
  data_FC2<-data[row.names(data)%in%FC2_sig,]
  data_FC1<-data[row.names(data)%in%FC1_sig,]
  data_top100<-data[row.names(data)%in%top100_sig,]
  tar <- unique(info$integrated_snn_res.0.3)
 #tar <- as.numeric(tar)
 tar <- tar[order(tar)]

 for(m in 1 : length(tar))
 {
  res <- matrix(0, nrow(data_FC2),4)
  row.names(res) <- row.names(data_FC2)
  colnames(res) <- c("C1_avg_expression","C1_median_expression","C2_avg_expression","C2_median_expression")
  
  tag <- info$integrated_snn_res.0.3 == tar[m]
  
  sam1 <- row.names(info)[tag==1]
  sam2 <- row.names(info)[tag==0]
  
  for(k in 1 : nrow(data_FC2))
  {
    cat("\r",m,"----->",k)
    xx1 <- as.numeric(data_FC2[k, sam1])
    xx2 <- as.numeric(data_FC2[k, sam2])
    
    res[k,"C1_avg_expression"] <- mean(xx1)
    res[k,"C2_avg_expression"] <- mean(xx2)
    res[k,"C1_median_expression"]<-median(xx1)
    res[k,"C2_median_expression"] <- median(xx2)
  }
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.3/scRNA_bulk/",file_nam[i],"_FC2_genes_cluster_",tar[m],"_vs_others_normalized_expression.xls")
  write.table(res,myoutf,sep="\t",quote=F)
 }

 for(m in 1 : length(tar))
 {
  res <- matrix(0, nrow(data_FC1),4)
  row.names(res) <- row.names(data_FC1)
  colnames(res) <- c("C1_avg_expression","C1_median_expression","C2_avg_expression","C2_median_expression")
  
  tag <- info$integrated_snn_res.0.3 == tar[m]
  
  sam1 <- row.names(info)[tag==1]
  sam2 <- row.names(info)[tag==0]
  
  for(k in 1 : nrow(data_FC1))
  {
    cat("\r",m,"----->",k)
    xx1 <- as.numeric(data_FC1[k, sam1])
    xx2 <- as.numeric(data_FC1[k, sam2])
    
    res[k,"C1_avg_expression"] <- mean(xx1)
    res[k,"C2_avg_expression"] <- mean(xx2)
    res[k,"C1_median_expression"]<-median(xx1)
    res[k,"C2_median_expression"] <- median(xx2)
  }
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.3/scRNA_bulk/",file_nam[i],"_FC1_genes_cluster_",tar[m],"_vs_others_normalized_expression.xls")
  write.table(res,myoutf,sep="\t",quote=F)
}

for(m in 1 : length(tar))
 {
  res <- matrix(0, nrow(data_top100),4)
  row.names(res) <- row.names(data_top100)
  colnames(res) <- c("C1_avg_expression","C1_median_expression","C2_avg_expression","C2_median_expression")
  
  tag <- info$integrated_snn_res.0.3 == tar[m]
  
  sam1 <- row.names(info)[tag==1]
  sam2 <- row.names(info)[tag==0]
  
  for(k in 1 : nrow(data_top100))
  {
    cat("\r",m,"----->",k)
    xx1 <- as.numeric(data_top100[k, sam1])
    xx2 <- as.numeric(data_top100[k, sam2])
    
    res[k,"C1_avg_expression"] <- mean(xx1)
    res[k,"C2_avg_expression"] <- mean(xx2)
    res[k,"C1_median_expression"]<-median(xx1)
    res[k,"C2_median_expression"] <- median(xx2)
  }
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.3/scRNA_bulk/",file_nam[i],"_top100_genes_cluster_",tar[m],"_vs_others_normalized_expression.xls")
  write.table(res,myoutf,sep="\t",quote=F)
}

}



dir <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/PM_KA_manuscript/Human_bulkRNA/analysis/marker/"
files <- list.files(dir)
file_nam <- gsub(".txt","",files)
for (k in 1:length(file_nam))
{
  myinf<-paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.3/scRNA_bulk/",file_nam[k],"_FC2_genes_cluster_0_vs_others_normalized_expression.xls")
  xx<-read.table(myinf,sep='\t')
  
  res<-matrix(0,nrow(xx),length(tar))
 row.names(res)<-row.names(xx)
 colnames(res)<-tar[order(tar)]
  for (i in 1:length(tar))
{
  myinf<-paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.3/scRNA_bulk/",file_nam[k],"_FC2_genes_cluster_",tar[i],"_vs_others_normalized_expression.xls")
  xx<-read.table(myinf,sep='\t')
  res[,i]<-xx$C1_avg_expression
}
res <- apply(res,1, function(arg) (arg-mean(arg))/sd(arg))
res <- t(res)


#res<-res[,c("2",'5','4','8','9','0','1','7','3','6')]

#ts<-res[sort(row.names(res)),]
breaklist_1<-seq(-2.5,0,0.2)
breaklist_2<-seq(0,2.6,0.2)
breaklist<-unique(c(breaklist_1,breaklist_2))

myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.3/scRNA_bulk/heatmaps/",file_nam[k],"_FC2_heatmap_everycluster.pdf")
pdf(myoutf,width=15,height=30)
pheatmap(mat=res,
         color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(length(breaklist)),
         breaks = breaklist,
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         border_color      = NA,
         show_colnames     = T,
         show_rownames     = T,
         drop_levels       = TRUE,
         scale='none',
         
         fontsize=20)
dev.off()


}






myinf1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.3/expression_data/CD8_normalized_data.xls"
myinf2 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.3/expression_data/CD8_normalized_cell_info.xls"

data<-read.table(myinf1,sep="\t",quote=NULL)
info<-read.table(myinf2,sep="\t",quote=NULL)

myinf3<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/PM_KA_manuscript/Mouse_PM/Meta_data/GATA6KO_vs_WT_DEG.txt"
myinf4<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/PM_KA_manuscript/Mouse_PM/Meta_data/GATA6KO_vs_WT_DEG.txt"



data3<-read.table(myinf3, sep="\t", header=T)
data4<-read.table(myinf4, sep="\t", header=T)

GATA6KO_sig<-toupper(data3[data3$Log2FC>1,'X'])
LPMWT_sig<-toupper(data3[data3$Log2FC<=-1,'X'])


##Mouse_to_human <- function(x){

   #human <- useEnsembl(biomart = "ensembl", 
   #                dataset = "hsapiens_gene_ensembl", 
   #                mirror = "uswest")
   #mouse <- useEnsembl(biomart = "ensembl", 
   #                dataset = "mmusculus_gene_ensembl", 
   #                mirror = "uswest")


#human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
#mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

#   genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
#
 #  humanx <- unique(genesV2[, 2])

 #  return(humanx)
#}

data_GATA6KO<-data[GATA6KO_sig,]
data_GATA6KO<-na.omit(data_GATA6KO)

tar <- unique(info$integrated_snn_res.0.3)
#tar <- as.numeric(tar)
tar <- tar[order(tar)]

for(i in 1 : length(tar))
{
  res <- matrix(0, nrow(data_GATA6KO),4)
  row.names(res) <- row.names(data_GATA6KO)
  colnames(res) <- c("C1_avg_expression","C1_median_expression","C2_avg_expression","C2_median_expression")
  
  tag <- info$integrated_snn_res.0.3 == tar[i]
  
  sam1 <- row.names(info)[tag==1]
  sam2 <- row.names(info)[tag==0]
  
  for(k in 1 : nrow(data_GATA6KO))
  {
    cat("\r",i,"----->",k)
    xx1 <- as.numeric(data_GATA6KO[k, sam1])
    xx2 <- as.numeric(data_GATA6KO[k, sam2])
    
    res[k,"C1_avg_expression"] <- mean(xx1)
    res[k,"C2_avg_expression"] <- mean(xx2)
    res[k,"C1_median_expression"]<-median(xx1)
    res[k,"C2_median_expression"] <- median(xx2)
  }
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.3/scRNA_bulk/GATA6KO_sig_",tar[i],"_vs_others_normalized_expression.xls")
  write.table(res,myoutf,sep="\t",quote=F)
}




data_LPMWT<-data[LPMWT_sig,]
data_LPMWT<-na.omit(data_LPMWT)

tar <- unique(info$integrated_snn_res.0.3)
#tar <- as.numeric(tar)
tar <- tar[order(tar)]

for(i in 1 : length(tar))
{
  res <- matrix(0, nrow(data_LPMWT),4)
  row.names(res) <- row.names(data_LPMWT)
  colnames(res) <- c("C1_avg_expression","C1_median_expression","C2_avg_expression","C2_median_expression")
  
  tag <- info$integrated_snn_res.0.3 == tar[i]
  
  sam1 <- row.names(info)[tag==1]
  sam2 <- row.names(info)[tag==0]
  
  for(k in 1 : nrow(data_LPMWT))
  {
    cat("\r",i,"----->",k)
    xx1 <- as.numeric(data_LPMWT[k, sam1])
    xx2 <- as.numeric(data_LPMWT[k, sam2])
    
    res[k,"C1_avg_expression"] <- mean(xx1)
    res[k,"C2_avg_expression"] <- mean(xx2)
    res[k,"C1_median_expression"]<-median(xx1)
    res[k,"C2_median_expression"] <- median(xx2)
  }
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.3/scRNA_bulk/data_LPMWT_sig_",tar[i],"_vs_others_normalized_expression.xls")
  write.table(res,myoutf,sep="\t",quote=F)
}






for (i in 1:length(tar))
{
  myinf<-paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.3/scRNA_bulk/GATA6KO_sig_",tar[i],"_vs_others_normalized_expression.xls")
  xx<-read.table(myinf,sep='\t')
  res[,i]<-xx$C1_avg_expression
}

res <- apply(res,1, function(arg) (arg-mean(arg))/sd(arg))
res <- t(res)


#res<-res[,c("2",'5','4','8','9','0','1','7','3','6')]

#ts<-res[sort(row.names(res)),]
breaklist_1<-seq(-2.5,0,0.2)
breaklist_2<-seq(0,2.6,0.2)
breaklist<-unique(c(breaklist_1,breaklist_2))

myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.3/scRNA_bulk/heatmaps/GATA6KO_signature.pdf")
pdf(myoutf,width=15,height=30)
pheatmap(mat=res,
         color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(length(breaklist)),
         breaks = breaklist,
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         border_color      = NA,
         show_colnames     = T,
         show_rownames     = T,
         drop_levels       = TRUE,
         scale='none',
         
         fontsize=20)
dev.off()



for (i in 1:length(tar))
{
  myinf<-paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.3/scRNA_bulk/data_LPMWT_sig_",tar[i],"_vs_others_normalized_expression.xls")
  xx<-read.table(myinf,sep='\t')
  res[,i]<-xx$C1_avg_expression
}

res <- apply(res,1, function(arg) (arg-mean(arg))/sd(arg))
res <- t(res)


#res<-res[,c("2",'5','4','8','9','0','1','7','3','6')]

#ts<-res[sort(row.names(res)),]
breaklist_1<-seq(-2.5,0,0.2)
breaklist_2<-seq(0,2.6,0.2)
breaklist<-unique(c(breaklist_1,breaklist_2))

myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.3/scRNA_bulk/heatmaps/LPMWT_signature.pdf")
pdf(myoutf,width=15,height=30)
pheatmap(mat=res,
         color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(length(breaklist)),
         breaks = breaklist,
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         border_color      = NA,
         show_colnames     = T,
         show_rownames     = T,
         drop_levels       = TRUE,
         scale='none',
         
         fontsize=20)
dev.off()












################heatmap single cell Seurat.
myoutf1 <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.3/cluster_identity/heatmap_marker_genes_different_phenotype_Tcell.pdf")


Idents(Myeloid_PW_Adult_kit.integrated)<-"integrated_snn_res.0.3"
 
tmp.markers <- FindAllMarkers(object = Myeloid_PW_Adult_kit.integrated, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.25,return.thresh=0.05)


top50 <- tmp.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)


pdf(myoutf1,35,35)
print(DoHeatmap(Myeloid_PW_Adult_kit.integrated, features = top50$gene) + 
        scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')), 
                              mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')), 
                              midpoint = 0, guide = "colourbar", aesthetics = "fill"))
dev.off()

tar<-c('0','1','3','6','7',"2",'5','4','8','9')
all_cell<-character()
for(i in 1:length(tar))
{
  cat("\r",i)
  len<-sum(Myeloid_PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.3==tar[i])
  tag<-sample(1:len,242)
  info<-Myeloid_PW_Adult_kit.integrated@meta.data[Myeloid_PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.3==tar[i],]
  barcode_use<-row.names(info)[tag]
  all_cell<-combine(all_cell,barcode_use)
}



myoutf1 <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.3/cluster_identity/heatmap_marker_genes_ordered_242cells_each_cluster.pdf")


Idents(Myeloid_PW_Adult_kit.integrated)<-"integrated_snn_res.0.3"

tmp.markers <- FindAllMarkers(object = Myeloid_PW_Adult_kit.integrated, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.15,return.thresh=0.05)

ts<-tmp.markers[order(factor(tmp.markers$cluster, levels = tar)),]

top30 <- ts %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)

Myeloid_PW_Adult_kit.integrated@meta.data[,"ordered_ident"]<-factor(Myeloid_PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.3,levels=c('0','1','3','6','7',"2",'5','4','8','9'))

pdf(myoutf1,12,10)
print(DoHeatmap(Myeloid_PW_Adult_kit.integrated, cells=all_cell,features = top30$gene,group.by="ordered_ident") + 
        scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')), 
                              mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')), 
                              midpoint = 0, guide = "colourbar", aesthetics = "fill"))
dev.off()

myoutf2 <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.3/cluster_identity/heatmap_marker_top30.xls")

write.table(top30,myoutf2,sep='\t')







tmp.markers <- FindAllMarkers(object = Myeloid_PW_Adult_kit.integrated,only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.3/cluster_identity/markers_each_cluster.xls")
write.table(tmp.markers,myoutf,sep="\t",quote=F)

tmp.markers <- FindAllMarkers(object = Myeloid_PW_Adult_kit.integrated,only.pos = T, logfc.threshold = 0.01)
myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.3/cluster_identity/markers_each_cluster_0.01.xls")
write.table(tmp.markers,myoutf,sep="\t",quote=F)


myinf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.3/cluster_identity/markers_each_cluster_0.01.xls")
markers_all<-read.table(myinf,sep="\t")
marker0<-subset(markers_all,cluster=="0")$gene
marker1<-subset(markers_all,cluster=="1")$gene

###Do venn diagram show overlap of differentially expressed genes between two Trm clusters

myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/PCx/res.0.3/cluster_identity/Venn_Two_Trm_Overlap.pdf")
pdf(myoutf,width=10,height=10)
draw.pairwise.venn(area1 = 109, area2 = 81,  cross.area = 45,category = c("C1-FOS","C2-IFNG"),
                 lwd = rep(2, 2), lty = rep("solid", 2), col =
                   c("#54278f","#fc4e2a"), fill = c("#54278f","#fc4e2a"), alpha = rep(0.5, 2))

dev.off()



myoutf <- paste0("/dartfs-hpc/rc/lab/T/TurkLab/Jichang/bulkTCR/Fred/tumor635.pdf")
pdf(myoutf,width=10,height=10)
draw.pairwise.venn(area1 = 2240, area2 = 2901,  cross.area = 1738,category = c("RNAseq","TCRseq"),
                   lwd = rep(2, 2), lty = rep("solid", 2), col =
                     c("#54278f","#fc4e2a"), fill = c("#54278f","#fc4e2a"), alpha = rep(0.5, 2))

dev.off()















































































































##############@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
########Subset T cell and analyze
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

 myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/Seruat_PW_Adult_kit.integrated__finished.RDS"
 PW_Adult_kit.integrated<-readRDS(myinf)
 
 tar<-"integrated_snn_res.0.6"
 Idents(object = PW_Adult_kit.integrated) <- tar
 subcluster<-c('3','5','6','8','9','10','12','13','14','16','19','20','24')
 tag<-PW_Adult_kit.integrated@meta.data[,tar]%in%subcluster
 tag1<-GetAssayData(PW_Adult_kit.integrated,assay='RNA', slot='data')['CD14',]>0.05
 cell_ID<-row.names(PW_Adult_kit.integrated@meta.data)[tag&!tag1]
 
 myinf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/Seruat_PW_PW_Adult_kit_Human__raw.RDS"
 xx<-readRDS(myinf)

 
 Tcell_PW_Adult_kit_Human.integrated<-subset(x=xx,cells=cell_ID)
 table(Tcell_PW_Adult_kit_Human.integrated@meta.data$orig.ident)
#    Kid1   Kid10   Kid11    Kid2    Kid4    Kid6    Kid7    Kid8    Kid9 PW_Adu1 
#   3351    3106    1599    3750    1801    1950    2330     452    1399    2750 
#PW_Adu2    PW_C    PW_Q    PW_S    PW_U    PW_W 
#   1161    6029    3556     516     860    2994 

 Tcell_PW_Adult_kit_Human.integrated.list <- SplitObject(Tcell_PW_Adult_kit_Human.integrated, split.by = "orig.ident")


 Tcell_PW_Adult_kit_Human.integrated.list <- lapply(X = Tcell_PW_Adult_kit_Human.integrated.list, FUN = function(x) {
    x <- NormalizeData(x,assay='RNA')
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 5000,assay='RNA')
 })


############################ Using RPCA to anchor############ 
 features <- SelectIntegrationFeatures(object.list = Tcell_PW_Adult_kit_Human.integrated.list,nfeatures = 4000)
 Tcell_PW_Adult_kit_Human.integrated.list <- lapply(X = Tcell_PW_Adult_kit_Human.integrated.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = T)
    x <- RunPCA(x, features = features, verbose = T)
 })

 immune.anchors <- FindIntegrationAnchors(object.list = Tcell_PW_Adult_kit_Human.integrated.list, anchor.features = features, reduction = "rpca")
 Tcell_PW_Adult_kit_Human.integrated <- IntegrateData(anchorset = immune.anchors)
 
 DefaultAssay(Tcell_PW_Adult_kit_Human.integrated) <- "integrated"



 Tcell_PW_Adult_kit_Human.integrated <- ScaleData(Tcell_PW_Adult_kit_Human.integrated, verbose = T)
  Tcell_PW_Adult_kit_Human.integrated <- RunPCA(Tcell_PW_Adult_kit_Human.integrated, npcs = 30, verbose = T)


 x<-c(1:15)


 eigs <- Tcell_PW_Adult_kit_Human.integrated[["pca"]]@stdev^2 #PCA的STBandard dev的平方#
  print(paste0("Variance captured by 35 PCs: ",sum(eigs[x] / sum(eigs)))) #0.856#




# Tcell_PW_Adult_kit_Human.integrated do clustering and TSNE plot No regresSTBissue ---------------------------------------------


 #[5]Do clustering
 Tcell_PW_Adult_kit_Human.integrated<-FindNeighbors(Tcell_PW_Adult_kit_Human.integrated, dims = x)
  Tcell_PW_Adult_kit_Human.integrated <- FindClusters(Tcell_PW_Adult_kit_Human.integrated,resolution = seq(0.1,1.0,0.1))

PrintFindClustersParams(object = Tcell_PW_Adult_kit_Human.integrated)

   clusters_resolution<-sapply(grep("^integrated_snn_res",colnames(Tcell_PW_Adult_kit_Human.integrated@meta.data),value = TRUE),
                            function(x) length(unique(Tcell_PW_Adult_kit_Human.integrated@meta.data[,x])))
 clusters_resolution


 #PC=x Tcell_PW_Adult_kit_Human.integrated
#integrated_snn_res.0.1 integrated_snn_res.0.2 integrated_snn_res.0.3 
#                     9                     13                     16 
#integrated_snn_res.0.4 integrated_snn_res.0.5 integrated_snn_res.0.6 
#                    18                     18                     20 
#integrated_snn_res.0.7 integrated_snn_res.0.8 integrated_snn_res.0.9 
#                    20                     21                     22 
#  integrated_snn_res.1 
#                    24  


 #Set up a loop to look for markers

 tar <- c("integrated_snn_res.0.1","integrated_snn_res.0.2","integrated_snn_res.0.3","integrated_snn_res.0.4",
         "integrated_snn_res.0.5","integrated_snn_res.0.6","integrated_snn_res.0.7","integrated_snn_res.0.8",
         "integrated_snn_res.0.9","integrated_snn_res.1")




 Tcell_PW_Adult_kit_Human.integrated <- RunUMAP(Tcell_PW_Adult_kit_Human.integrated, dims= x)

 #Tcell_PW_Adult_kit_Human.integrated<-RunTSNE(Tcell_PW_Adult_kit_Human.integrated,reduction = "pca",dims=x)

 for(i in 1:length(tar))
 {
  Idents(object = Tcell_PW_Adult_kit_Human.integrated) <- tar[i]
  myoutf1 <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/UMAP_SB6_",tar[i],"_undefined.pdf")
  pdf(myoutf1,7,5)
  print(DimPlot(Tcell_PW_Adult_kit_Human.integrated, reduction = "umap",pt.size=0.2))
  dev.off()
 }


 #for(i in 1:length(tar))
 #{
 # Idents(object = Tcell_PW_Adult_kit_Human.integrated) <- tar[i]
 # myoutf1 <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/TSNE_SB6_",tar[i],"_undefined.pdf")
 # pdf(myoutf1,7,5)
 # print(DimPlot(Tcell_PW_Adult_kit_Human.integrated, reduction = "tsne"))
 # dev.off()
 #}


 myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/TSNE_Tcell_PW_Adult_kit_Human.integrated_groupby_orig.ident.pdf"
  pdf(myoutf,width=8,height=5)
 DimPlot(object = Tcell_PW_Adult_kit_Human.integrated,group.by='orig.ident',reduction="umap",pt.size=0.2)
 dev.off()

gender_info<-c('F','M','F','M','F','M','M','M','F','F','F','F','M','F','F','M')
age_info<-c('24','7','10','1.3','1.2','14','15','8','17','7','29','49','30','36','36','39')
tar<-unique(Tcell_PW_Adult_kit_Human.integrated@meta.data$orig.ident)
Tcell_PW_Adult_kit_Human.integrated@meta.data[,'gender']<-rep(nrow(Tcell_PW_Adult_kit_Human.integrated@meta.data),'Unknown')
Tcell_PW_Adult_kit_Human.integrated@meta.data[,'age']<-rep(nrow(Tcell_PW_Adult_kit_Human.integrated@meta.data),'Unknown')


Tcell_PW_Adult_kit_Human.integrated@meta.data[,'group']<-rep('Unassigned',nrow(Tcell_PW_Adult_kit_Human.integrated@meta.data))
tag1<-as.numeric(Tcell_PW_Adult_kit_Human.integrated@meta.data$age)>18&Tcell_PW_Adult_kit_Human.integrated@meta.data$gender=='M'
tag2<-as.numeric(Tcell_PW_Adult_kit_Human.integrated@meta.data$age)>18&Tcell_PW_Adult_kit_Human.integrated@meta.data$gender=='F'
tag3<-as.numeric(Tcell_PW_Adult_kit_Human.integrated@meta.data$age)<=18&Tcell_PW_Adult_kit_Human.integrated@meta.data$gender=='M'
tag4<-as.numeric(Tcell_PW_Adult_kit_Human.integrated@meta.data$age)<=18&Tcell_PW_Adult_kit_Human.integrated@meta.data$gender=='F'
Tcell_PW_Adult_kit_Human.integrated@meta.data$group[tag1]<-'Adult_Male'
Tcell_PW_Adult_kit_Human.integrated@meta.data$group[tag2]<-'Adult_Female'
Tcell_PW_Adult_kit_Human.integrated@meta.data$group[tag3]<-'Kid_Male'
Tcell_PW_Adult_kit_Human.integrated@meta.data$group[tag4]<-'Kid_Female'

for (i in 1:length(tar))
{
  tag<-Tcell_PW_Adult_kit_Human.integrated@meta.data$orig.ident==tar[i]
  Tcell_PW_Adult_kit_Human.integrated@meta.data[tag,'gender']<-gender_info[i]
  Tcell_PW_Adult_kit_Human.integrated@meta.data[tag,'age']<-age_info[i]
}


 myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/Seruat_Tcell_PW_Adult_kit_Human.integrated__finished.RDS"
 saveRDS(Tcell_PW_Adult_kit_Human.integrated, file = myoutf)



DefaultAssay(Tcell_PW_Adult_kit_Human.integrated) <- "integrated"
tar1 <- c("integrated_snn_res.0.6","integrated_snn_res.0.5","integrated_snn_res.0.4","integrated_snn_res.0.7","integrated_snn_res.0.8")

 for(i in 1 : length(tar1))
 {
  cat("\r",i)
  Idents(object = Tcell_PW_Adult_kit_Human.integrated) <- tar1[i]
  tmp <- Tcell_PW_Adult_kit_Human.integrated
  tmp.markers <- FindAllMarkers(object = tmp, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.3,assay='integrated')
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/cluster_markers_",tar1[i],"_ScRNA.xls")
  write.table(tmp.markers,myoutf,sep="\t",quote=F)
  
  cell.list <- WhichCells(tmp, downsample = 200)
  myoutf1 <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/heatmap_markers_",tar1[i],"_.pdf")
  pdf(myoutf1,30,35)
  top20 <- tmp.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
  print(DoHeatmap(tmp, cells=cell.list,features = top20$gene,assay='integrated') + NoLegend())
  dev.off()
 }



#########################Find DEG of each cluster using findconserve markers#########
 tar<-"integrated_snn_res.0.3"
 Idents(object = Tcell_PW_Adult_kit_Human.integrated) <- tar
 tar1<-unique(as.numeric(Tcell_PW_Adult_kit_Human.integrated@meta.data$"integrated_snn_res.0.3"))
 tar1<-sort(tar1-1)
 DEG<-data.frame()
 for (i in 1:length(tar1))
 {
  markers <- FindConservedMarkers(Tcell_PW_Adult_kit_Human.integrated, ident.1 = tar1[i], grouping.var = "orig.ident", verbose = T,only.pos = T,logfc.threshold = 0.25,min.cells.group=1)
  markers[,"cluster"]<-tar1[i]
  markers[,'gene']<-row.names(markers)
  tar2<-grep('log2FC',colnames(markers))
  markers[,'average_Log2FC']<-apply(markers[,colnames(markers)[tar2]],1,mean)
  tag<-markers[,'max_pval']<=0.05
  markers<-markers[tag,]
  markers<-markers[order(markers[,'average_Log2FC'],decreasing=TRUE),]
  DEG<-rbind(DEG,markers)
 }





myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/Seruat_Tcell_PW_Adult_kit_Human.integrated__finished.RDS"
Tcell_PW_Adult_kit_Human.integrated<-readRDS(myinf)

 DefaultAssay(Tcell_PW_Adult_kit_Human.integrated) <- "RNA"

data_RNA<-Tcell_PW_Adult_kit_Human.integrated[['RNA']]@data
data_integrated<-Tcell_PW_Adult_kit_Human.integrated[['integrated']]@data
data_RNA<-data_RNA[Immune_general,]
data_integrated<-data_integrated[Immune_general,]
cells<-sample(colnames(data_RNA),20)


 DefaultAssay(Tcell_PW_Adult_kit_Human.integrated) <- "RNA"
Immune_general<-c('CD3G','CD3E','CD8A','CD4','ITGAM','KLRB1','CD69','ITGAE','IFNG','IL4','IL13','SELL','S1PR1','CX3CR1','IL17A')

p<-FeaturePlot(object = Tcell_PW_Adult_kit_Human.integrated, features = Immune_general,cols= c('#deebf7', "#bd0026"),
               reduction = "umap",slot = "data",max.cutoff=4,pt.size=0.3,combine=F)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + theme(plot.title = element_text(size=40),
                          axis.title = element_text(family="Helvetica",face = "bold", color = "black",size=20), 
                           axis.text=element_text(family="Helvetica",face = "bold", color = "black",size=16),
                           axis.line = element_line(colour = 'black', size = 1),
                           axis.ticks = element_line(colour = "black", size = 1),
                           axis.ticks.length=unit(.3, "cm"))
}
 
myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/Immune_general_marker.pdf"
pdf(myoutf,width=50,height=30)
print(cowplot::plot_grid(plotlist = p,ncol=5))
dev.off()






tar1 <- c("integrated_snn_res.0.3","integrated_snn_res.0.4")

 for(i in 1 : length(tar1))
 {
  cat("\r",i)
  Idents(object = Tcell_PW_Adult_kit_Human.integrated) <- tar1[i]
  tmp <- Tcell_PW_Adult_kit_Human.integrated
  tmp.markers <- FindAllMarkers(object = tmp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/cluster_markers_",tar1[i],"_ScRNA_RNAassay.xls")
  write.table(tmp.markers,myoutf,sep="\t",quote=F)
 }


myeloid_general<-c('SELL','IRF4','F5','CSF1R','VCAN','MERTK','MRC1','CD14','ITGAX','ITGAM','GATA6','CD1C','LYVE1','CCR2','TIMD4','SELP','CX3CR1','CD163',
'FLT3','MARCO','IL3RA','AXL','SIGLEC6','RETN','CD226',"CLEC4C","CLEC9A",'CD209','VSIG4')

flow_panel<-c('PTPRC','CD14','ITGAM','ITGAX','LYVE1','TIMD4','CX3CR1','SELP','CCR2','HLA-DRA','CCR7','CD209','MARCO','MRC1','FOLR2',
'LILRB4','CD1C','CLEC9A','SIGLEC6','AXL','IL3RA','CD274','CLEC4C','SIRPA','CD163')

important<-c('GATA6','SELP','TIMD4','LYVE1','MARCO','CCR2','CD1C','CD14','PDCD1LG2','IL3RA','MRC1','CD163')

p<-FeaturePlot(object = Tcell_PW_Adult_kit_Human.integrated, features = myeloid_general,cols= c('#deebf7', "#bd0026"),
               reduction = "umap",slot = "data",max.cutoff=4,pt.size=0.3,combine=F,order = TRUE)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + theme(plot.title = element_text(size=40),
                          axis.title = element_text(family="Helvetica",face = "bold", color = "black",size=20), 
                           axis.text=element_text(family="Helvetica",face = "bold", color = "black",size=16),
                           axis.line = element_line(colour = 'black', size = 1),
                           axis.ticks = element_line(colour = "black", size = 1),
                           axis.ticks.length=unit(.3, "cm"))
}
 
myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/Myeloid_marker.pdf"
pdf(myoutf,width=50,height=50)
print(cowplot::plot_grid(plotlist = p,ncol=5))
dev.off()




p<-FeaturePlot(object = Tcell_PW_Adult_kit_Human.integrated, features = flow_panel,cols= c('#deebf7', "#bd0026"),
               reduction = "umap",slot = "data",max.cutoff=4,pt.size=0.3,combine=F,order = TRUE)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + theme(plot.title = element_text(size=40),
                          axis.title = element_text(family="Helvetica",face = "bold", color = "black",size=20), 
                           axis.text=element_text(family="Helvetica",face = "bold", color = "black",size=16),
                           axis.line = element_line(colour = 'black', size = 1),
                           axis.ticks = element_line(colour = "black", size = 1),
                           axis.ticks.length=unit(.3, "cm"))
}
 
myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/Flow_panel_marker.pdf"
pdf(myoutf,width=50,height=50)
print(cowplot::plot_grid(plotlist = p,ncol=5))
dev.off()




p<-FeaturePlot(object = Tcell_PW_Adult_kit_Human.integrated, features = important,cols= c('#deebf7', "#bd0026"),
               reduction = "umap",slot = "data",max.cutoff=4,pt.size=0.3,combine=F,order = TRUE)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + theme(plot.title = element_text(size=50),
                          axis.title = element_text(family="Helvetica",face = "bold", color = "black",size=30), 
                           axis.text=element_text(family="Helvetica",face = "bold", color = "black",size=20),
                           axis.line = element_line(colour = 'black', size = 1),
                           axis.ticks = element_line(colour = "black", size = 1),
                           axis.ticks.length=unit(.3, "cm"))
}
 
myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/important_marker.pdf"
pdf(myoutf,width=65,height=50)
print(cowplot::plot_grid(plotlist = p,ncol=4))
dev.off()


####Plot adult and children myeloid cell markers seperately seperating sex

tag1<-as.numeric(Tcell_PW_Adult_kit_Human.integrated@meta.data$age)>18&Tcell_PW_Adult_kit_Human.integrated@meta.data$gender=='M'
tag2<-as.numeric(Tcell_PW_Adult_kit_Human.integrated@meta.data$age)>18&Tcell_PW_Adult_kit_Human.integrated@meta.data$gender=='F'
tag3<-as.numeric(Tcell_PW_Adult_kit_Human.integrated@meta.data$age)<=18&Tcell_PW_Adult_kit_Human.integrated@meta.data$gender=='M'
tag4<-as.numeric(Tcell_PW_Adult_kit_Human.integrated@meta.data$age)<=18&Tcell_PW_Adult_kit_Human.integrated@meta.data$gender=='F'


p<-FeaturePlot(object = Tcell_PW_Adult_kit_Human.integrated, cells =row.names(Tcell_PW_Adult_kit_Human.integrated@meta.data)[tag1], features = myeloid_general,cols= c('#deebf7', "#bd0026"),
               reduction = "umap",slot = "data",max.cutoff=4,pt.size=0.3,combine=F,order = TRUE)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + theme(plot.title = element_text(size=40),
                          axis.title = element_text(family="Helvetica",face = "bold", color = "black",size=20), 
                           axis.text=element_text(family="Helvetica",face = "bold", color = "black",size=16),
                           axis.line = element_line(colour = 'black', size = 1),
                           axis.ticks = element_line(colour = "black", size = 1),
                           axis.ticks.length=unit(.3, "cm"))
}
 
myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/Myeloid_marker_adult_male.pdf"
pdf(myoutf,width=50,height=50)
print(cowplot::plot_grid(plotlist = p,ncol=5))
dev.off()



p<-FeaturePlot(object = Tcell_PW_Adult_kit_Human.integrated, cells =row.names(Tcell_PW_Adult_kit_Human.integrated@meta.data)[tag2], features = myeloid_general,cols= c('#deebf7', "#bd0026"),
               reduction = "umap",slot = "data",max.cutoff=4,pt.size=0.3,combine=F,order = TRUE)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + theme(plot.title = element_text(size=40),
                          axis.title = element_text(family="Helvetica",face = "bold", color = "black",size=20), 
                           axis.text=element_text(family="Helvetica",face = "bold", color = "black",size=16),
                           axis.line = element_line(colour = 'black', size = 1),
                           axis.ticks = element_line(colour = "black", size = 1),
                           axis.ticks.length=unit(.3, "cm"))
}
 
myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/Myeloid_marker_adult_female.pdf"
pdf(myoutf,width=50,height=50)
print(cowplot::plot_grid(plotlist = p,ncol=5))
dev.off()




p<-FeaturePlot(object = Tcell_PW_Adult_kit_Human.integrated, cells =row.names(Tcell_PW_Adult_kit_Human.integrated@meta.data)[tag3], features = myeloid_general,cols= c('#deebf7', "#bd0026"),
               reduction = "umap",slot = "data",max.cutoff=4,pt.size=0.3,combine=F,order = TRUE)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + theme(plot.title = element_text(size=40),
                          axis.title = element_text(family="Helvetica",face = "bold", color = "black",size=20), 
                           axis.text=element_text(family="Helvetica",face = "bold", color = "black",size=16),
                           axis.line = element_line(colour = 'black', size = 1),
                           axis.ticks = element_line(colour = "black", size = 1),
                           axis.ticks.length=unit(.3, "cm"))
}
 
myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/Myeloid_marker_kid_male.pdf"
pdf(myoutf,width=50,height=50)
print(cowplot::plot_grid(plotlist = p,ncol=5))
dev.off()



p<-FeaturePlot(object = Tcell_PW_Adult_kit_Human.integrated, cells =row.names(Tcell_PW_Adult_kit_Human.integrated@meta.data)[tag4], features = myeloid_general,cols= c('#deebf7', "#bd0026"),
               reduction = "umap",slot = "data",max.cutoff=4,pt.size=0.3,combine=F,order = TRUE)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + theme(plot.title = element_text(size=40),
                          axis.title = element_text(family="Helvetica",face = "bold", color = "black",size=20), 
                           axis.text=element_text(family="Helvetica",face = "bold", color = "black",size=16),
                           axis.line = element_line(colour = 'black', size = 1),
                           axis.ticks = element_line(colour = "black", size = 1),
                           axis.ticks.length=unit(.3, "cm"))
}
 
myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/Myeloid_marker_kid_female.pdf"
pdf(myoutf,width=50,height=50)
print(cowplot::plot_grid(plotlist = p,ncol=5))
dev.off()






#############Cluster contribution by each different patient##########
   
   Idents(Tcell_PW_Adult_kit_Human.integrated)<-"integrated_snn_res.0.5"
 
 xx<-table(Tcell_PW_Adult_kit_Human.integrated@meta.data$integrated_snn_res.0.3,Tcell_PW_Adult_kit_Human.integrated@meta.data$orig.ident)
 yy<-t(apply(xx,1,function(x){x/sum(x)})   )
 myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.2/cluster_patient_distribution/patient_distribution_to_each_cluster.xls"
 write.table(yy,myoutf,sep="\t",quote=F)
 
 data<-matrix(0,ncol(yy)*nrow(yy),3)
 colnames(data)<-c("cluster","patient","value")

 data<-as.data.frame(data)
 data$cluster<-rep(row.names(yy),ncol(yy))

 data$patient<-as.character(sapply(colnames(yy),function(x){rep(x,nrow(yy))}))

 for (i in 1:nrow(data))
 {
  number<-yy[data$cluster[i],data$patient[i]]
  data$value[i]<-number
 }


 data$cluster<-factor(data$cluster,levels=unique(data$cluster))
 data$patient<-factor(data$patient,levels=c("PW_normal","PW_CRC1","PW_CRC2",'Kid1','Kid2','Kid4','Kid6','Kid7','Kid8','Kid9','Kid10'))

  col<-c("#2171b5","#cb181d","#74c476")

 myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.3/cluster_patient_distribution/stackplot_patient_contribution_each_cluster_proportion.pdf"
 pdf(myoutf,width=8,height=8)
 ggplot() + geom_bar(aes(y = value, x = cluster, fill = patient), data = data,
                    stat="identity")+theme(axis.text.x = element_text(angle = 60, hjust = 1,size=10),
                                           axis.title.x = element_blank())

 dev.off()






myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/Seruat_Tcell_PW_Adult_kit_Human.integrated__finished.RDS"
Tcell_PW_Adult_kit_Human.integrated<-readRDS(myinf)

tar<-unique(Tcell_PW_Adult_kit_Human.integrated@meta.data$integrated_snn_res.0.5)
tar<-tar[order(tar)]
chi<-as.matrix(table(Tcell_PW_Adult_kit_Human.integrated@meta.data$orig.ident,Tcell_PW_Adult_kit_Human.integrated@meta.data$integrated_snn_res.0.5))

chi_square<-chisq.test(chi)

chi_observed<-chi_square$observed
chi_expected<-chi_square$expected
ROE<-chi_observed/chi_expected

myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.5/cluster_patient_distribution/ROE.xls"
write.table(ROE,myoutf,sep="\t",quote=F)



each_cluster_all<-apply(chi,2,sum)
chi_enrichment<-chi
for (i in 1:nrow(chi))
{
    for (k in 1:ncol(chi))
    {
      chi_enrichment[i,k]<-chi[i,k]*sum(chi)/sum(chi[i,])/sum(chi[,k])
    }
}

myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.5/cluster_patient_distribution/hyper_enrichment.xls"
write.table(chi_enrichment,myoutf,sep="\t",quote=F)











myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/Seruat_Tcell_PW_Adult_kit_Human.integrated__finished.RDS"
Tcell_PW_Adult_kit_Human.integrated<-readRDS(myinf)

tar<-unique(Tcell_PW_Adult_kit_Human.integrated@meta.data$integrated_snn_res.0.5)
tar<-tar[order(tar)]
chi<-as.matrix(table(Tcell_PW_Adult_kit_Human.integrated@meta.data$orig.ident,Tcell_PW_Adult_kit_Human.integrated@meta.data$integrated_snn_res.0.5))

chi_square<-chisq.test(chi)

chi_observed<-chi_square$observed
chi_expected<-chi_square$expected
ROE<-chi_observed/chi_expected

myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.5/cluster_patient_distribution/ROE.xls"
write.table(ROE,myoutf,sep="\t",quote=F)



each_cluster_all<-apply(chi,2,sum)
chi_enrichment<-chi
for (i in 1:nrow(chi))
{
    for (k in 1:ncol(chi))
    {
      chi_enrichment[i,k]<-chi[i,k]*sum(chi)/sum(chi[i,])/sum(chi[,k])
    }
}

myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.7/cluster_patient_distribution/hyper_enrichment.xls"
write.table(chi_enrichment,myoutf,sep="\t",quote=F)




#################Output cell of each orig.ident and adult versus kid#############
myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/Seruat_Tcell_PW_Adult_kit_Human.integrated__finished.RDS"
Tcell_PW_Adult_kit_Human.integrated<-readRDS(myinf)

Idents(Tcell_PW_Adult_kit_Human.integrated)<-"integrated_snn_res.0.5"
 tar<-unique(Tcell_PW_Adult_kit_Human.integrated@meta.data$orig.ident)
 for (i in 1:length(tar))
 {
    tag<-Tcell_PW_Adult_kit_Human.integrated@meta.data$orig.ident==tar[i]
    cell_ID<-row.names(Tcell_PW_Adult_kit_Human.integrated@meta.data)[tag]
    myoutf1 <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.5/each_sample/UMAP_SB6_",tar[i],"_undefined.pdf")
  pdf(myoutf1,7,5)
  print(DimPlot(Tcell_PW_Adult_kit_Human.integrated, reduction = "umap",pt.size=0.1,cells =cell_ID))
  dev.off()
 
 }
  
  



myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/Seruat_Tcell_PW_Adult_kit_Human.integrated__finished.RDS"
Tcell_PW_Adult_kit_Human.integrated<-readRDS(myinf)

Idents(Tcell_PW_Adult_kit_Human.integrated)<-"integrated_snn_res.0.5"
cell_ID_adult<-row.names(Tcell_PW_Adult_kit_Human.integrated@meta.data)[as.numeric(Tcell_PW_Adult_kit_Human.integrated@meta.data$age)>18]
cell_ID_kid<-row.names(Tcell_PW_Adult_kit_Human.integrated@meta.data)[as.numeric(Tcell_PW_Adult_kit_Human.integrated@meta.data$age)<18]


    myoutf1 <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.5/each_sample/All_Adult.pdf")
  pdf(myoutf1,7,5)
  print(DimPlot(Tcell_PW_Adult_kit_Human.integrated, reduction = "umap",pt.size=0.1,cells =cell_ID_adult))
  dev.off()
 

myoutf1 <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.5/each_sample/All_Kid.pdf")
  pdf(myoutf1,7,5)
  print(DimPlot(Tcell_PW_Adult_kit_Human.integrated, reduction = "umap",pt.size=0.1,cells =cell_ID_kid))
  dev.off()
 
  
  







################Find differentially expressed genes between adult, kid, male and female even in the same cluster###############################

tar <- unique(Tcell_PW_Adult_kit_Human.integrated@meta.data$group)
tar1<-combn(tar,2)
tar2<-tar1[,c(2,3,4,5)]

Idents(object = Tcell_PW_Adult_kit_Human.integrated) <- 'integrated_snn_res.0.5'
 DefaultAssay(Tcell_PW_Adult_kit_Human.integrated) <- "integrated"




 for(i in 1 : ncol(tar2))
 {
  cat("\r",i)
  tmp <- Tcell_PW_Adult_kit_Human.integrated
  tmp.markers <- FindMarkers(object = tmp, ident.1=tar2[,i][1],ident.2 =tar2[,i][2],min.pct = 0.1, logfc.threshold = 0.1)
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/DEG_",tar2[,i][1],'_',tar2[,i][2],".xls")
  write.table(tmp.markers,myoutf,sep="\t",quote=F)
 }

tar3<-as.character(unique(Tcell_PW_Adult_kit_Human.integrated@meta.data$integrated_snn_res.0.5))


for(i in 1 : ncol(tar2))
 {
  cat("\r",i)
  markers_all<-data.frame()
  for (j in 1:length(tar3))
  {
  tmp <- subset(Tcell_PW_Adult_kit_Human.integrated,idents=tar3[j])
  Idents(object = tmp) <- 'group'
  tmp.markers <- FindMarkers(object = tmp, ident.1=tar2[,i][1],ident.2 =tar2[,i][2],min.pct = 0.1, logfc.threshold = 0.1)
  tmp.markers[,'cluster']<-rep(tar3[j],nrow(tmp.markers))
  markers_all<-rbind(markers_all,tmp.markers)
  }
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/DEG_",tar2[,i][1],'_',tar2[,i][2],".xls")
  write.table(markers_all,myoutf,sep="\t",quote=F)
 }











Tcell_PW_Adult_kit_Human.integrated@meta.data[,'age_group_only']<-rep('unknown',nrow(Tcell_PW_Adult_kit_Human.integrated@meta.data))
tag1<-as.numeric(Tcell_PW_Adult_kit_Human.integrated@meta.data$age)>=18
tag2<-as.numeric(Tcell_PW_Adult_kit_Human.integrated@meta.data$age)<18
Tcell_PW_Adult_kit_Human.integrated@meta.data$age_group_only[tag1]<-'Adult'
Tcell_PW_Adult_kit_Human.integrated@meta.data$age_group_only[tag2]<-'Kid'


  markers_all<-data.frame()
  for (j in 1:length(tar3))
  {
  tmp <- subset(Tcell_PW_Adult_kit_Human.integrated,idents=tar3[j])
  Idents(object = tmp) <- 'age_group_only'
  tmp.markers <- FindMarkers(object = tmp, ident.1='Adult',ident.2 ='Kid',min.pct = 0.1, logfc.threshold = 0.1)
  tmp.markers[,'cluster']<-rep(tar3[j],nrow(tmp.markers))
  markers_all<-rbind(markers_all,tmp.markers)
  }
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/DEG_Adult_vs_Kid.xls")
  write.table(markers_all,myoutf,sep="\t",quote=F)







################Cell cycle scoring###################
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/Seruat_Tcell_PW_Adult_kit_Human.integrated__finished.RDS"
Tcell_PW_Adult_kit_Human.integrated<-readRDS(myinf)

DefaultAssay(Tcell_PW_Adult_kit_Human.integrated) <- "RNA"
Tcell_PW_Adult_kit_Human.integrated<-CellCycleScoring(Tcell_PW_Adult_kit_Human.integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
cell_cycle_score<-t(Tcell_PW_Adult_kit_Human.integrated@meta.data[,c("S.Score",'G2M.Score')])
adt_assay <- CreateAssayObject(counts = cell_cycle_score)
Tcell_PW_Adult_kit_Human.integrated[['cell_cycle']]<-adt_assay
DefaultAssay(Tcell_PW_Adult_kit_Human.integrated) <- "cell_cycle"
p<-FeaturePlot(object = Tcell_PW_Adult_kit_Human.integrated, features = c('S.Score','G2M.Score'),cols= c('#deebf7', "#bd0026"),
               reduction = "umap",slot = "data",max.cutoff=1.5,pt.size=0.3,combine=F)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + theme(plot.title = element_text(size=40),
                          axis.title = element_text(family="Helvetica",face = "bold", color = "black",size=20), 
                           axis.text=element_text(family="Helvetica",face = "bold", color = "black",size=16),
                           axis.line = element_line(colour = 'black', size = 1),
                           axis.ticks = element_line(colour = "black", size = 1),
                           axis.ticks.length=unit(.3, "cm"))
}
 
myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.5/cell_cycle/Cell_Cycle_Score.pdf"
pdf(myoutf,width=25,height=10)
print(cowplot::plot_grid(plotlist = p,ncol=2))
dev.off()

myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.5/cell_cycle/By_cell_cycle_Tcell_PW_Adult_kit_Human.integrated.pdf"
  pdf(myoutf,width=8,height=5)
 DimPlot(object = Tcell_PW_Adult_kit_Human.integrated,group.by='Phase',reduction="umap",pt.size=0.3)
 dev.off()

myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/Seruat_Tcell_PW_Adult_kit_Human.integrated__finished.RDS"
 saveRDS(Tcell_PW_Adult_kit_Human.integrated, file = myoutf)

Idents(Tcell_PW_Adult_kit_Human.integrated)<-"integrated_snn_res.0.5"
 
myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.5/cell_cycle/Cell_cycle_violin.pdf"
  pdf(myoutf,width=13,height=5)
 VlnPlot(object = Tcell_PW_Adult_kit_Human.integrated,assay='cell_cycle',features = c('S.Score','G2M.Score'),pt.size=0.3)
 dev.off()




Myeloid_marker<-c('XCR1','IL3RA','CLEC4C','SIGLEC6','ITGAM','ITGAX','MRC1','CD14','CD68','GATA6','CD1C','LYVE1','CCR2','TIMD4','CLEC9A')


p<-FeaturePlot(object = Tcell_PW_Adult_kit_Human.integrated, features = Myeloid_marker,cols= c('#deebf7', "#bd0026"),
               reduction = "umap",slot = "data",max.cutoff=4,pt.size=0.3,combine=F)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + theme(plot.title = element_text(size=40),
                          axis.title = element_text(family="Helvetica",face = "bold", color = "black",size=20), 
                           axis.text=element_text(family="Helvetica",face = "bold", color = "black",size=16),
                           axis.line = element_line(colour = 'black', size = 1),
                           axis.ticks = element_line(colour = "black", size = 1),
                           axis.ticks.length=unit(.3, "cm"))
}
 
myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/Myeloid_marker.pdf"
pdf(myoutf,width=50,height=30)
print(cowplot::plot_grid(plotlist = p,ncol=5))
dev.off()
















##########First get average gene expression of each cluster and then do z transformation to do heatmap for selected genes###########

myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/Seruat_Tcell_PW_Adult_kit_Human.integrated__finished.RDS"
exp<-readRDS(myinf)
DefaultAssay(exp) <- "RNA"

myoutf1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.5/expression_data/CD8_normalized_data.xls"
myoutf2 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.5/expression_data/CD8_normalized_cell_info.xls"


info <- exp@meta.data
res <- GetAssayData(object = exp,slot='data')
res <- as.matrix(res)
row.names(res) <- toupper(row.names(res))
avg <- apply(res,1, mean)
tag <- avg >0
res <- res[tag,]

res <- unique(res)

tag <- duplicated(row.names(res))

dup_gene <- row.names(res)[tag]

tag <- row.names(res) == dup_gene
#tag==logical(0) Thus ignore the res1,res STBep
#res1 <- res[tag==0,]

#res2 <- res[tag==1,]
#avg <- apply(res2,1, mean)

#res <- rbind(res1, res2[2,])
res <- as.data.frame(res)

write.table(res,myoutf1,sep="\t",quote=F)
write.table(info,myoutf2,sep="\t",quote=F)






myinf1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.5/expression_data/CD8_normalized_data.xls"
myinf2 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.5/expression_data/CD8_normalized_cell_info.xls"

data<-read.table(myinf1,sep="\t",quote=NULL)
info<-read.table(myinf2,sep="\t",quote=NULL)

DC_markers<-c('TCF4','IL3RA','CLEC4C','LILRB4','FCER1A','FCER1G','LILRA4',
            'IRF8','IRF4','ID2','BATF3','CLEC9A','XCR1','BTLA',
            'CD163','S100A9','S100A8','VCAN','FCN1','CD14','BST1','CD36',
            'FCGR2B','CD1C','HLA-DQA2','HLA-DQA1','HLA-DQB1','HLA-DQP1',
            'ZEB2','KLF4','ITGAX','ITGAM','CD2','SIRPA','CLEC4A','CLEC10A',
            'AXL','SIGLEC6','CX3CR1','SIGLEC1')

PM_markers<-c('GARNL3','CAV1','HGF','LAPTM5','GATA6','CD9','TGFB2','APOC1','F5','HDC','LYZ1',
               'FSCN1','ID1','PLXDC2','SELP','CXCL13','FABP5','ADGRE1','GAS6','ID2',
               'MRC1','FOLR2','LYVE1','MARCO','CD163','CD209','CD74','FPR1',
               'ANXA5','CBR2','CLEC4A1','MAF','ITGB1','CXCL13',
               'FABP7','CCR2','CSF1R','SELL','ITGAX','HLA-DQA2','HLA-DQA1','HLA-DQB1','HLA-DQP1','ITGAM',
               'RCN3','FCGRT','CCL24','IRF7','EMP1','LGALS1','HAL','RPS18')



#cluster_markers<-c('S100A9','S100A8','VCAN','FCN1','CD14','MAFB','ITGAM','SPP1',
#                    'FCER1A','CD1C','CD1E','HLA-DPB1','HLA-DQB1','HLA-DRA','CLEC10A','FCGR2B','CLEC4A','TCF4','IL3RA','ITGAX','ZEB2',
#                    'LILRB4','CSF1R','GSN','GPR183','RGS1','IL1B','FOS','NR4A1',
#                    'C1QC','C1QA','C1QB','FN1','LYVE1','FOLR2','MARCO','CD163','FCGR3A','TIMD4','APOE','CD68',
#                    'IL32','CCL2','GZMA','CCL5','NKG7','GNLY','ID3','STAT1','CD36','CCL20','ITGA4','CD3G','CD8A',
#                    'CCL4','CXCL10','CCL3','CXCL9','CXCL11',
#                    'CXCR6','GATA3','RORA','TIGIT','CD69',
#                    'MKI67','TK1','CCNA2','TUBA1B',
#                    'IRF8','ID2','BATF3','CLEC9A','XCR1','BTLA',
#                     'NEAT1','CIITA','NFKBIZ','FOXP1','TLR2','IRF4'
#                    )

data_DC<-data[DC_markers,]
data_DC<-na.omit(data_DC)

tar <- unique(info$integrated_snn_res.0.5)
#tar <- as.numeric(tar)
tar <- tar[order(tar)]

for(i in 1 : length(tar))
{
  res <- matrix(0, nrow(data_DC),4)
  row.names(res) <- row.names(data_DC)
  colnames(res) <- c("C1_avg_expression","C1_median_expression","C2_avg_expression","C2_median_expression")
  
  tag <- info$integrated_snn_res.0.5 == tar[i]
  
  sam1 <- row.names(info)[tag==1]
  sam2 <- row.names(info)[tag==0]
  
  for(k in 1 : nrow(data_DC))
  {
    cat("\r",i,"----->",k)
    xx1 <- as.numeric(data_DC[k, sam1])
    xx2 <- as.numeric(data_DC[k, sam2])
    
    res[k,"C1_avg_expression"] <- mean(xx1)
    res[k,"C2_avg_expression"] <- mean(xx2)
    res[k,"C1_median_expression"]<-median(xx1)
    res[k,"C2_median_expression"] <- median(xx2)
  }
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.5/cluster_identity/pheatmap/Human_DC_gene_",tar[i],"_vs_others_normalized_expression.xls")
  write.table(res,myoutf,sep="\t",quote=F)
}




res<-matrix(0,nrow(data_DC),length(tar))
 row.names(res)<-row.names(data_DC)
 colnames(res)<-tar[order(tar)]



for (i in 1:length(tar))
{
  myinf<-paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.5/cluster_identity/pheatmap/Human_DC_gene_",tar[i],"_vs_others_normalized_expression.xls")
  xx<-read.table(myinf,sep='\t')
  res[,i]<-xx$C1_avg_expression
}

res <- apply(res,1, function(arg) (arg-mean(arg))/sd(arg))
res <- t(res)


#res<-res[,c("2",'5','4','8','9','0','1','7','3','6')]

#ts<-res[sort(row.names(res)),]
breaklist_1<-seq(-2.5,0,0.2)
breaklist_2<-seq(0,2.6,0.2)
breaklist<-unique(c(breaklist_1,breaklist_2))

myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.5/cluster_identity/pheatmap/heatmap/Human_DC_gene_each_cluster.pdf")
pdf(myoutf,width=25,height=15)
pheatmap(mat=res,
         color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(length(breaklist)),
         breaks = breaklist,
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         border_color      = NA,
         show_colnames     = T,
         show_rownames     = T,
         drop_levels       = TRUE,
         scale='none',
         
         fontsize=20)
dev.off()












data_PM<-data[PM_markers,]
data_PM<-na.omit(data_PM)

tar <- unique(info$integrated_snn_res.0.5)
#tar <- as.numeric(tar)
tar <- tar[order(tar)]

for(i in 1 : length(tar))
{
  res <- matrix(0, nrow(data_PM),4)
  row.names(res) <- row.names(data_PM)
  colnames(res) <- c("C1_avg_expression","C1_median_expression","C2_avg_expression","C2_median_expression")
  
  tag <- info$integrated_snn_res.0.5 == tar[i]
  
  sam1 <- row.names(info)[tag==1]
  sam2 <- row.names(info)[tag==0]
  
  for(k in 1 : nrow(data_PM))
  {
    cat("\r",i,"----->",k)
    xx1 <- as.numeric(data_PM[k, sam1])
    xx2 <- as.numeric(data_PM[k, sam2])
    
    res[k,"C1_avg_expression"] <- mean(xx1)
    res[k,"C2_avg_expression"] <- mean(xx2)
    res[k,"C1_median_expression"]<-median(xx1)
    res[k,"C2_median_expression"] <- median(xx2)
  }
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.5/cluster_identity/pheatmap/Mouse_Peritoneal_MAC_gene_",tar[i],"_vs_others_normalized_expression.xls")
  write.table(res,myoutf,sep="\t",quote=F)
}




res<-matrix(0,nrow(data_PM),length(tar))
 row.names(res)<-row.names(data_PM)
 colnames(res)<-tar[order(tar)]



for (i in 1:length(tar))
{
  myinf<-paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.5/cluster_identity/pheatmap/Mouse_Peritoneal_MAC_gene_",tar[i],"_vs_others_normalized_expression.xls")
  xx<-read.table(myinf,sep='\t')
  res[,i]<-xx$C1_avg_expression
}

res <- apply(res,1, function(arg) (arg-mean(arg))/sd(arg))
res <- t(res)


#res<-res[,c("2",'5','4','8','9','0','1','7','3','6')]

#ts<-res[sort(row.names(res)),]
breaklist_1<-seq(-2.5,0,0.2)
breaklist_2<-seq(0,2.6,0.2)
breaklist<-unique(c(breaklist_1,breaklist_2))

myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.5/cluster_identity/pheatmap/heatmap/Mouse_Peritoneal_MAC_gene_each_cluster.pdf")
pdf(myoutf,width=15,height=15)
pheatmap(mat=res,
         color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(length(breaklist)),
         breaks = breaklist,
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         border_color      = NA,
         show_colnames     = T,
         show_rownames     = T,
         drop_levels       = TRUE,
         scale='none',
         
         fontsize=20)
dev.off()










data_each_cluster<-data[cluster_markers,]
data_each_cluster<-na.omit(data_each_cluster)

tar <- unique(info$integrated_snn_res.0.3)
#tar <- as.numeric(tar)
tar <- tar[order(tar)]

for(i in 1 : length(tar))
{
  res <- matrix(0, nrow(data_each_cluster),4)
  row.names(res) <- row.names(data_each_cluster)
  colnames(res) <- c("C1_avg_expression","C1_median_expression","C2_avg_expression","C2_median_expression")
  
  tag <- info$integrated_snn_res.0.3 == tar[i]
  
  sam1 <- row.names(info)[tag==1]
  sam2 <- row.names(info)[tag==0]
  
  for(k in 1 : nrow(data_each_cluster))
  {
    cat("\r",i,"----->",k)
    xx1 <- as.numeric(data_each_cluster[k, sam1])
    xx2 <- as.numeric(data_each_cluster[k, sam2])
    
    res[k,"C1_avg_expression"] <- mean(xx1)
    res[k,"C2_avg_expression"] <- mean(xx2)
    res[k,"C1_median_expression"]<-median(xx1)
    res[k,"C2_median_expression"] <- median(xx2)
  }
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.3/cluster_identity/pheatmap/marker_each_cluster_",tar[i],"_vs_others_normalized_expression.xls")
  write.table(res,myoutf,sep="\t",quote=F)
}




res<-matrix(0,nrow(data_each_cluster),length(tar))
 row.names(res)<-row.names(data_each_cluster)
 colnames(res)<-tar[order(tar)]



for (i in 1:length(tar))
{
  myinf<-paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.3/cluster_identity/pheatmap/marker_each_cluster_",tar[i],"_vs_others_normalized_expression.xls")
  xx<-read.table(myinf,sep='\t')
  res[,i]<-xx$C1_avg_expression
}

res <- apply(res,1, function(arg) (arg-mean(arg))/sd(arg))
res <- t(res)


#res<-res[,c("2",'5','4','8','9','0','1','7','3','6')]

#ts<-res[sort(row.names(res)),]
breaklist_1<-seq(-2.5,0,0.2)
breaklist_2<-seq(0,2.6,0.2)
breaklist<-unique(c(breaklist_1,breaklist_2))

myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.3/cluster_identity/pheatmap/heatmap/marker_each_cluster.pdf")
pdf(myoutf,width=8,height=17)
pheatmap(mat=res,
         color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(length(breaklist)),
         breaks = breaklist,
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         border_color      = NA,
         show_colnames     = T,
         show_rownames     = T,
         drop_levels       = TRUE,
         scale='none',
         
         fontsize=14)
dev.off()








data_headtmap_paper<-data[Paper_markers,]
data_headtmap_paper<-na.omit(data_headtmap_paper)

tar <- unique(info$integrated_snn_res.0.5)
#tar <- as.numeric(tar)
tar <- tar[order(tar)]

for(i in 1 : length(tar))
{
  res <- matrix(0, nrow(data_PM),4)
  row.names(res) <- row.names(data_PM)
  colnames(res) <- c("C1_avg_expression","C1_median_expression","C2_avg_expression","C2_median_expression")
  
  tag <- info$integrated_snn_res.0.5 == tar[i]
  
  sam1 <- row.names(info)[tag==1]
  sam2 <- row.names(info)[tag==0]
  
  for(k in 1 : nrow(data_PM))
  {
    cat("\r",i,"----->",k)
    xx1 <- as.numeric(data_PM[k, sam1])
    xx2 <- as.numeric(data_PM[k, sam2])
    
    res[k,"C1_avg_expression"] <- mean(xx1)
    res[k,"C2_avg_expression"] <- mean(xx2)
    res[k,"C1_median_expression"]<-median(xx1)
    res[k,"C2_median_expression"] <- median(xx2)
  }
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.5/cluster_identity/pheatmap/Mouse_Peritoneal_MAC_gene_",tar[i],"_vs_others_normalized_expression.xls")
  write.table(res,myoutf,sep="\t",quote=F)
}




res<-matrix(0,nrow(data_PM),length(tar))
 row.names(res)<-row.names(data_PM)
 colnames(res)<-tar[order(tar)]



for (i in 1:length(tar))
{
  myinf<-paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.5/cluster_identity/pheatmap/Mouse_Peritoneal_MAC_gene_",tar[i],"_vs_others_normalized_expression.xls")
  xx<-read.table(myinf,sep='\t')
  res[,i]<-xx$C1_avg_expression
}

res <- apply(res,1, function(arg) (arg-mean(arg))/sd(arg))
res <- t(res)


#res<-res[,c("2",'5','4','8','9','0','1','7','3','6')]

#ts<-res[sort(row.names(res)),]
breaklist_1<-seq(-2.5,0,0.2)
breaklist_2<-seq(0,2.6,0.2)
breaklist<-unique(c(breaklist_1,breaklist_2))

myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.5/cluster_identity/pheatmap/heatmap/Mouse_Peritoneal_MAC_gene_each_cluster.pdf")
pdf(myoutf,width=15,height=15)
pheatmap(mat=res,
         color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(length(breaklist)),
         breaks = breaklist,
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         border_color      = NA,
         show_colnames     = T,
         show_rownames     = T,
         drop_levels       = TRUE,
         scale='none',
         
         fontsize=20)
dev.off()


























PM_markers<-c('GARNL3','CAV1','HGF','LAPTM5','GATA6','CD9','TGFB2','APOC1','F5','HDC','LYZ1',
               'FSCN1','ID1','PLXDC2','SELP','CXCL13','FABP5','ADGRE1','GAS6','ID2',
               'MRC1','FOLR2','LYVE1','MARCO','CD163','CD209','CD74','FPR1',
               'ANXA5','CBR2','CLEC4A1','MAF','ITGB1','PF4',
               'FABP7','CCR2','CSF1R','SELL','ITGAX','HLA-DQA2','HLA-DQA1','HLA-DQB1','HLA-DQP1','ITGAM',
               'RCN3','FCGRT','CCL24','IRF7','EMP1','LGALS1','HAL','RPS18')



PM_markers<-c('GARNL3','CAV1','HGF','LAPTM5','GATA6','CD9','TGFB2','APOC1','F5','HDC','LYZ1',
               'FSCN1','ID1','PLXDC2','SELP','CXCL13','FABP5','ADGRE1','GAS6','ID2',
               'MRC1','FOLR2','LYVE1','MARCO','CD163','CD209','CD74','FPR1',
               'ANXA5','CBR2','CLEC4A1','MAF','ITGB1','CXCL13',
               'FABP7','CCR2','CSF1R','SELL','ITGAX','HLA-DQA2','HLA-DQA1','HLA-DQB1','HLA-DQP1','ITGAM',
               'RCN3','FCGRT','CCL24','IRF7','EMP1','LGALS1','HAL','RPS18')



myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.5/cluster_identity/Dotplot/Dot_plot_MAC_heterogeneity.pdf"
pdf(myoutf,width=15,height=6)
print(DotPlot(Tcell_PW_Adult_kit_Human.integrated, cols=c("#08519c","#ef3b2c"),group.by="integrated_snn_res.0.5",features = PM_markers)+
        theme(axis.text.x = element_text(angle = 45,size=10,hjust = 1)))

dev.off()


DC_markers<-c('TCF4','IL3RA','CLEC4C','LILRB4','FCER1A','FCER1G','LILRA4',
            'IRF8','IRF4','ID2','BATF3','CLEC9A','XCR1','BTLA',
            'CD163','S100A9','S100A8','VCAN','FCN1','CD14','BST1','CD36',
            'FCGR2B','CD1C','HLA-DQA2','HLA-DQA1','HLA-DQB1','HLA-DQP1',
            'ZEB2','KLF4','ITGAX','ITGAM','CD2','SIRPA','CLEC4A','CLEC10A',
            'AXL','SIGLEC6','CX3CR1','SIGLEC1')


myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.3/cluster_identity/Dotplot/Dot_plot_DC_heterogeneity.pdf"
pdf(myoutf,width=15,height=6)
print(DotPlot(Tcell_PW_Adult_kit_Human.integrated, cols=c("#08519c","#ef3b2c"),group.by="integrated_snn_res.0.5",features = DC_markers)+
        theme(axis.text.x = element_text(angle = 45,size=10,hjust = 1)))

dev.off()





Idents(object = Tcell_PW_Adult_kit_Human.integrated) <- "integrated_snn_res.0.5"
  tmp <- Tcell_PW_Adult_kit_Human.integrated
  tmp.markers <- FindAllMarkers(object = tmp, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.3,assay='integrated')


tar<-sort(unique(Tcell_PW_Adult_kit_Human.integrated@meta.data$integrated_snn_res.0.5))
all_cell<-character()
for(i in 1:(length(tar)-2))
{
  cat("\r",i)
  len<-sum(Tcell_PW_Adult_kit_Human.integrated@meta.data$integrated_snn_res.0.5==tar[i])
  tag<-sample(1:len,245)
  info<-Tcell_PW_Adult_kit_Human.integrated@meta.data[Tcell_PW_Adult_kit_Human.integrated$integrated_snn_res.0.5==tar[i],]
  barcode_use<-row.names(info)[tag]
  all_cell<-c(all_cell,barcode_use)
}
cluster15_16_barcode<-row.names(Tcell_PW_Adult_kit_Human.integrated@meta.data)[Tcell_PW_Adult_kit_Human.integrated@meta.data$integrated_snn_res.0.5%in%c('15','16')]

all_cell<-c(all_cell,cluster15_16_barcode)

myoutf1 <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.5/cluster_identity/heatmap_marker_genes_ordered_245cells_each_cluster.pdf")

tmp.markers <- tmp.markers[tmp.markers$p_val_adj<0.05,]

ts<-tmp.markers[order(factor(tmp.markers$cluster, levels = tar)),]

top20 <- ts %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

Tcell_PW_Adult_kit_Human.integrated@meta.data[,"ordered_ident"]<-factor(Tcell_PW_Adult_kit_Human.integrated@meta.data$integrated_snn_res.0.5,
                                                                    levels=seq(0,16,1))

pdf(myoutf1,12,10)
print(DoHeatmap(Tcell_PW_Adult_kit_Human.integrated, cells=all_cell,features = top20$gene,group.by="integrated_snn_res.0.5",assay='integrated') + 
        scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')), 
                              mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')), 
                              midpoint = 0, guide = "colourbar", aesthetics = "fill")+NoLegend())
dev.off()






##########Check GATA6################
GATA6_exp<-as.matrix(GetAssayData(Tcell_PW_Adult_kit_Human.integrated,assay='RNA',slot='data')[c('GATA6','SELP'),])
apply(GATA6_exp,1,summary)

#              GATA6        SELP
#Min.    0.000000000 0.000000000
#1st Qu. 0.000000000 0.000000000
#Median  0.000000000 0.000000000
#Mean    0.009995213 0.007311954
#3rd Qu. 0.000000000 0.000000000
#Max.    2.931832980 2.733619776

tag1<-GATA6_exp[1,]>mean(GATA6_exp[1,])
tag2<-GATA6_exp[2,]>mean(GATA6_exp[2,])

GATA6_pos<-colnames(GATA6_exp)[tag1]
SELP_pos<-colnames(GATA6_exp)[tag2]

Tcell_PW_Adult_kit_Human.integrated@meta.data[,'GATA6']<-rep("Neg",nrow(Tcell_PW_Adult_kit_Human.integrated@meta.data))
Tcell_PW_Adult_kit_Human.integrated@meta.data[,'SELP']<-rep("Neg",nrow(Tcell_PW_Adult_kit_Human.integrated@meta.data))
Tcell_PW_Adult_kit_Human.integrated@meta.data[GATA6_pos,'GATA6']<-'Pos'
Tcell_PW_Adult_kit_Human.integrated@meta.data[SELP_pos,'SELP']<-'Pos'

GATA6_dis<-table(Tcell_PW_Adult_kit_Human.integrated@meta.data[,'GATA6'],Tcell_PW_Adult_kit_Human.integrated@meta.data[,'orig.ident'])
SELP_dis<-table(Tcell_PW_Adult_kit_Human.integrated@meta.data[,'SELP'],Tcell_PW_Adult_kit_Human.integrated@meta.data[,'orig.ident'])



chi_square<-chisq.test(GATA6_dis)

chi_observed<-chi_square$observed
chi_expected<-chi_square$expected
ROE_GATA6<-chi_observed/chi_expected


chi_square<-chisq.test(SELP_dis)

chi_observed<-chi_square$observed
chi_expected<-chi_square$expected
ROE_SELP<-chi_observed/chi_expected



GATA6_enrichment<-GATA6_dis
for (i in 1:nrow(GATA6_dis))
{
    for (k in 1:ncol(GATA6_dis))
    {
      GATA6_enrichment[i,k]<-GATA6_dis[i,k]*sum(GATA6_dis)/sum(GATA6_dis[i,])/sum(GATA6_dis[,k])
    }
}



myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.7/cluster_patient_distribution/GATA6_enrichment.xls"
write.table(ROE_GATA6,myoutf,sep="\t",quote=F)

myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.7/cluster_patient_distribution/SELP_enrichment.xls"
write.table(ROE_SELP,myoutf,sep="\t",quote=F)





##############Find doublets in the dataset#############
 
 myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/Seruat_Tcell_PW_Adult_kit_Human.integrated__finished.RDS"
 Tcell_PW_Adult_kit_Human.integrated<-readRDS(myinf)
 cell_ID<-row.names(Tcell_PW_Adult_kit_Human.integrated@meta.data)

 myinf1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/raw/Nan_PW_Normal"
 PW_normal.data <- Read10X(data.dir=myinf1)
 PW_normal <- CreateSeuratObject(counts = PW_normal.data, project = "PW_normal",min.cells=3,min.features=200)
 PW_normal<-RenameCells(object=PW_normal,add.cell.id=unique(PW_normal@meta.data$orig.ident))
 PW_normal<-RenameCells(object=PW_normal,new.names=gsub('-1','',Cells(PW_normal)))


  myinf2 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/raw/Nan_PW_CRC1/sample_feature_bc_matrix"
 PW_CRC1.data <- Read10X(data.dir=myinf2)
 PW_CRC1 <- CreateSeuratObject(counts = PW_CRC1.data, project = "PW_CRC1",min.cells=3,min.features=200)
 PW_CRC1<-RenameCells(object=PW_CRC1,add.cell.id=unique(PW_CRC1@meta.data$orig.ident))
 PW_CRC1<-RenameCells(object=PW_CRC1,new.names=gsub('-1','',Cells(PW_CRC1)))


 myinf3 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/raw/Nan_PW_CRC2/sample_feature_bc_matrix"
 PW_CRC2.data <- Read10X(data.dir=myinf3)
 PW_CRC2 <- CreateSeuratObject(counts = PW_CRC2.data, project = "PW_CRC2",min.cells=3,min.features=200)
 PW_CRC2<-RenameCells(object=PW_CRC2,add.cell.id=unique(PW_CRC2@meta.data$orig.ident))
 PW_CRC2<-RenameCells(object=PW_CRC2,new.names=gsub('-1','',Cells(PW_CRC2)))

 xx<- merge(PW_normal, y = c(PW_CRC1, PW_CRC2), project = "PW_Adult_kit")


 Tcell_PW_Adult_kit_Human.doublet<-subset(x=xx,cells=cell_ID)
 
 Tcell_PW_Adult_kit_Human.doublet.list <- SplitObject(Tcell_PW_Adult_kit_Human.doublet, split.by = "orig.ident")

 Tcell_PW_Adult_kit_Human.doublet.list <- lapply(X = Tcell_PW_Adult_kit_Human.doublet.list, FUN = function(x) {
    x <- NormalizeData(x,assay='RNA')
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000,assay='RNA')
    x <- ScaleData(x, verbose = T)
    x <- RunPCA(x,verbose = T)
    x<-FindNeighbors(x, dims = 1:15)
   x <- FindClusters(x,resolution = 0.3)
 
 })









 Tcell_PW_Adult_kit_Human.doublet.list <- lapply(X = Tcell_PW_Adult_kit_Human.doublet.list, FUN = function(x) {
    yy <- paramSweep_v3(x, PCs = 1:15, sct = FALSE)
    zz<-summarizeSweep(yy, GT = FALSE)
    tt<-find.pK(zz)
    pK<-tt %>% filter(BCmetric == max(BCmetric)) %>% select (pK)
    pK<- as.numeric (as.character(pK[[1]]))
    homotypic.prop <- modelHomotypic(x@meta.data$"integrated_snn_res.0.3")
    nExp_poi <- round(0.055*nrow(x@meta.data)) 
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    x <- doubletFinder_v3(x, PCs = 1:15, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
 })


 Tcell_PW_Adult_kit_Human.doublet<-merge(
  Tcell_PW_Adult_kit_Human.doublet.list[[1]],
  y =c(Tcell_PW_Adult_kit_Human.doublet.list[[2]], Tcell_PW_Adult_kit_Human.doublet.list[[3]]),
  add.cell.ids = NULL,
  project = "Tcell_PW_Adult_kit_Human.doublet"
)

tag<-Tcell_PW_Adult_kit_Human.doublet@meta.data[,"DF.classifications_0.25_0.17_448"]=='Doublet'
doublet_ID<-row.names(Tcell_PW_Adult_kit_Human.doublet@meta.data[tag,])

Tcell_PW_Adult_kit_Human.integrated@meta.data[,"doublet"]<-rep('Singlet',nrow(Tcell_PW_Adult_kit_Human.integrated@meta.data))
Tcell_PW_Adult_kit_Human.integrated@meta.data[doublet_ID,"doublet"]<-"Doublet"

myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/doublet_finder.pdf"
  pdf(myoutf,width=8,height=5)
 DimPlot(object = Tcell_PW_Adult_kit_Human.integrated,group.by='doublet',reduction="umap",pt.size=0.3,cells.highlight=doublet_ID,cols.highlight = "#e6550d",
  sizes.highlight = 0.3)
 dev.off()














####################Perform GSEA########

screen -r 182535


myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/Seruat_Tcell_PW_Adult_kit_Human.integrated__finished.RDS"
Tcell_PW_Adult_kit_Human.integrated<-readRDS(myinf)
DefaultAssay(Tcell_PW_Adult_kit_Human.integrated) <- "RNA"
tar<-"integrated_snn_res.0.5"
 Idents(object = Tcell_PW_Adult_kit_Human.integrated) <- tar
 tar1<-unique(as.numeric(Tcell_PW_Adult_kit_Human.integrated@meta.data$"integrated_snn_res.0.5"))
 tar1<-sort(tar1-1)
 for (i in 1:length(tar1))
 {
  markers <- FindConservedMarkers(Tcell_PW_Adult_kit_Human.integrated, ident.1 = tar1[i], grouping.var = "orig.ident", verbose = T,logfc.threshold = 0.001,min.cells.group=1,min.pct=0.001)
  markers[,"cluster"]<-tar1[i]
  markers[,'gene']<-row.names(markers)
  tar2<-grep('log2FC',colnames(markers))
  markers[,'average_Log2FC']<-apply(markers[,colnames(markers)[tar2]],1,mean)
  markers<-markers[order(markers[,'average_Log2FC'],decreasing=TRUE),]
  myoutf<-paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.5/GSEA/marker/Cluster_",tar1[i],"differentially_expressed_markers.xls")
  write.table(markers,myoutf,quote=F,sep="\t")
 }


dir <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/PM_KA_manuscript/Human_bulkRNA/analysis/marker/"
files <- list.files(dir)
file_nam <- gsub(".xls","",files)
for(i in 1:length(files))
{
  cat("\r",i)
  tmpinf <- paste0(dir,files[i])
  res <- read.table(tmpinf,sep="\t",quote=NULL)
  res <- res[order(res$Log2FC,decreasing=T),]
  
  xx <- res$Log2FC
  names(xx) <- row.names(res)
  names(xx) <- toupper(names(xx))
  
  tag = duplicated(names(xx))
  res = res[tag==0,]
  row.names(res) <- toupper(row.names(res))
  res[,"Name"] <- row.names(res)
  res[,"metric"] <- res$Log2FC
  res <- res[,c("Name","metric")]
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/PM_KA_manuscript/Human_bulkRNA/analysis/GSEA/",file_nam[i],".rnk")
  write.table(res,file=myoutf,quote=F,sep="\t",row.names=F)
}





dir <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.5/GSEA/marker/"
files <- list.files(dir)
file_nam <- gsub(".xls","",files)
for(i in 1:length(files))
{
  cat("\r",i)
  tmpinf <- paste0(dir,files[i])
  res <- read.table(tmpinf,sep="\t",quote=NULL)
  #tag<-res$average_Log2FC>0
  #res<-res[tag,]
  res <- res[order(res$average_Log2FC,decreasing=T),]
  
  xx <- res$average_Log2FC
  names(xx) <- row.names(res)
  names(xx) <- toupper(names(xx))
  
  tag = duplicated(names(xx))
  res = res[tag==0,]
  row.names(res) <- toupper(row.names(res))
  res[,"Name"] <- row.names(res)
  res[,"metric"] <- res$average_Log2FC
  res <- res[,c("Name","metric")]
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.5/GSEA/preranked_GSEA/",file_nam[i],".rnk")
  write.table(res,file=myoutf,quote=F,sep="\t",row.names=F)
}























##########First get average gene expression of each cluster and then do z transformation to do heatmap for selected genes###########

myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/Seruat_Tcell_PW_Adult_kit_Human.integrated__finished.RDS"
exp <- readRDS(myinf)
DefaultAssay(exp) <- "integrated"

myoutf1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.3/expression_data/CD8_normalized_data.xls"
myoutf2 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.3/expression_data/CD8_normalized_cell_info.xls"


info <- exp@meta.data
res <- GetAssayData(object = exp)
res <- as.matrix(res)
row.names(res) <- toupper(row.names(res))
avg <- apply(res,1, mean)
tag <- avg >0
res <- res[tag,]

res <- unique(res)

tag <- duplicated(row.names(res))

dup_gene <- row.names(res)[tag]

tag <- row.names(res) == dup_gene
#tag==logical(0) Thus ignore the res1,res STBep
#res1 <- res[tag==0,]

#res2 <- res[tag==1,]
#avg <- apply(res2,1, mean)

#res <- rbind(res1, res2[2,])
res <- as.data.frame(res)

write.table(res,myoutf1,sep="\t",quote=F)
write.table(info,myoutf2,sep="\t",quote=F)

myinf1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.3/expression_data/CD8_normalized_data.xls"
myinf2 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.3/expression_data/CD8_normalized_cell_info.xls"

data<-read.table(myinf1,sep="\t",quote=NULL)
info<-read.table(myinf2,sep="\t",quote=NULL)


dir <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/PM_KA_manuscript/Human_bulkRNA/analysis/marker/"
files <- list.files(dir)
file_nam <- gsub(".txt","",files)
for(i in 1:length(files))
{
  cat("\r",i)
  tmpinf <- paste0(dir,files[i])
  res1 <- read.table(tmpinf,sep="\t",quote=NULL)
  #tag<-res$average_Log2FC>0
  #res<-res[tag,]
  res1 <- res1[order(res1$Log2FC,decreasing=T),]
  FC2_sig<-row.names(res1[res1$Log2FC>2,])
  FC1_sig<-row.names(res1[res1$Log2FC>1,])
  top100_sig<-row.names(res1[1:100,])
  
  data_FC2<-data[row.names(data)%in%FC2_sig,]
  data_FC1<-data[row.names(data)%in%FC1_sig,]
  data_top100<-data[row.names(data)%in%top100_sig,]
  tar <- unique(info$integrated_snn_res.0.3)
 #tar <- as.numeric(tar)
 tar <- tar[order(tar)]

 for(m in 1 : length(tar))
 {
  res <- matrix(0, nrow(data_FC2),4)
  row.names(res) <- row.names(data_FC2)
  colnames(res) <- c("C1_avg_expression","C1_median_expression","C2_avg_expression","C2_median_expression")
  
  tag <- info$integrated_snn_res.0.3 == tar[m]
  
  sam1 <- row.names(info)[tag==1]
  sam2 <- row.names(info)[tag==0]
  
  for(k in 1 : nrow(data_FC2))
  {
    cat("\r",m,"----->",k)
    xx1 <- as.numeric(data_FC2[k, sam1])
    xx2 <- as.numeric(data_FC2[k, sam2])
    
    res[k,"C1_avg_expression"] <- mean(xx1)
    res[k,"C2_avg_expression"] <- mean(xx2)
    res[k,"C1_median_expression"]<-median(xx1)
    res[k,"C2_median_expression"] <- median(xx2)
  }
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.3/scRNA_bulk/",file_nam[i],"_FC2_genes_cluster_",tar[m],"_vs_others_normalized_expression.xls")
  write.table(res,myoutf,sep="\t",quote=F)
 }

 for(m in 1 : length(tar))
 {
  res <- matrix(0, nrow(data_FC1),4)
  row.names(res) <- row.names(data_FC1)
  colnames(res) <- c("C1_avg_expression","C1_median_expression","C2_avg_expression","C2_median_expression")
  
  tag <- info$integrated_snn_res.0.3 == tar[m]
  
  sam1 <- row.names(info)[tag==1]
  sam2 <- row.names(info)[tag==0]
  
  for(k in 1 : nrow(data_FC1))
  {
    cat("\r",m,"----->",k)
    xx1 <- as.numeric(data_FC1[k, sam1])
    xx2 <- as.numeric(data_FC1[k, sam2])
    
    res[k,"C1_avg_expression"] <- mean(xx1)
    res[k,"C2_avg_expression"] <- mean(xx2)
    res[k,"C1_median_expression"]<-median(xx1)
    res[k,"C2_median_expression"] <- median(xx2)
  }
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.3/scRNA_bulk/",file_nam[i],"_FC1_genes_cluster_",tar[m],"_vs_others_normalized_expression.xls")
  write.table(res,myoutf,sep="\t",quote=F)
}

for(m in 1 : length(tar))
 {
  res <- matrix(0, nrow(data_top100),4)
  row.names(res) <- row.names(data_top100)
  colnames(res) <- c("C1_avg_expression","C1_median_expression","C2_avg_expression","C2_median_expression")
  
  tag <- info$integrated_snn_res.0.3 == tar[m]
  
  sam1 <- row.names(info)[tag==1]
  sam2 <- row.names(info)[tag==0]
  
  for(k in 1 : nrow(data_top100))
  {
    cat("\r",m,"----->",k)
    xx1 <- as.numeric(data_top100[k, sam1])
    xx2 <- as.numeric(data_top100[k, sam2])
    
    res[k,"C1_avg_expression"] <- mean(xx1)
    res[k,"C2_avg_expression"] <- mean(xx2)
    res[k,"C1_median_expression"]<-median(xx1)
    res[k,"C2_median_expression"] <- median(xx2)
  }
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.3/scRNA_bulk/",file_nam[i],"_top100_genes_cluster_",tar[m],"_vs_others_normalized_expression.xls")
  write.table(res,myoutf,sep="\t",quote=F)
}

}



dir <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/PM_KA_manuscript/Human_bulkRNA/analysis/marker/"
files <- list.files(dir)
file_nam <- gsub(".txt","",files)
for (k in 1:length(file_nam))
{
  myinf<-paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.3/scRNA_bulk/",file_nam[k],"_FC2_genes_cluster_0_vs_others_normalized_expression.xls")
  xx<-read.table(myinf,sep='\t')
  
  res<-matrix(0,nrow(xx),length(tar))
 row.names(res)<-row.names(xx)
 colnames(res)<-tar[order(tar)]
  for (i in 1:length(tar))
{
  myinf<-paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.3/scRNA_bulk/",file_nam[k],"_FC2_genes_cluster_",tar[i],"_vs_others_normalized_expression.xls")
  xx<-read.table(myinf,sep='\t')
  res[,i]<-xx$C1_avg_expression
}
res <- apply(res,1, function(arg) (arg-mean(arg))/sd(arg))
res <- t(res)


#res<-res[,c("2",'5','4','8','9','0','1','7','3','6')]

#ts<-res[sort(row.names(res)),]
breaklist_1<-seq(-2.5,0,0.2)
breaklist_2<-seq(0,2.6,0.2)
breaklist<-unique(c(breaklist_1,breaklist_2))

myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.3/scRNA_bulk/heatmaps/",file_nam[k],"_FC2_heatmap_everycluster.pdf")
pdf(myoutf,width=15,height=30)
pheatmap(mat=res,
         color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(length(breaklist)),
         breaks = breaklist,
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         border_color      = NA,
         show_colnames     = T,
         show_rownames     = T,
         drop_levels       = TRUE,
         scale='none',
         
         fontsize=20)
dev.off()


}






myinf1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.3/expression_data/CD8_normalized_data.xls"
myinf2 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.3/expression_data/CD8_normalized_cell_info.xls"

data<-read.table(myinf1,sep="\t",quote=NULL)
info<-read.table(myinf2,sep="\t",quote=NULL)

myinf3<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/PM_KA_manuscript/Mouse_PM/Meta_data/GATA6KO_vs_WT_DEG.txt"
myinf4<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/PM_KA_manuscript/Mouse_PM/Meta_data/GATA6KO_vs_WT_DEG.txt"



data3<-read.table(myinf3, sep="\t", header=T)
data4<-read.table(myinf4, sep="\t", header=T)

GATA6KO_sig<-toupper(data3[data3$Log2FC>1,'X'])
LPMWT_sig<-toupper(data3[data3$Log2FC<=-1,'X'])


##Mouse_to_human <- function(x){

   #human <- useEnsembl(biomart = "ensembl", 
   #                dataset = "hsapiens_gene_ensembl", 
   #                mirror = "uswest")
   #mouse <- useEnsembl(biomart = "ensembl", 
   #                dataset = "mmusculus_gene_ensembl", 
   #                mirror = "uswest")


#human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
#mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

#   genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
#
 #  humanx <- unique(genesV2[, 2])

 #  return(humanx)
#}

data_GATA6KO<-data[GATA6KO_sig,]
data_GATA6KO<-na.omit(data_GATA6KO)

tar <- unique(info$integrated_snn_res.0.3)
#tar <- as.numeric(tar)
tar <- tar[order(tar)]

for(i in 1 : length(tar))
{
  res <- matrix(0, nrow(data_GATA6KO),4)
  row.names(res) <- row.names(data_GATA6KO)
  colnames(res) <- c("C1_avg_expression","C1_median_expression","C2_avg_expression","C2_median_expression")
  
  tag <- info$integrated_snn_res.0.3 == tar[i]
  
  sam1 <- row.names(info)[tag==1]
  sam2 <- row.names(info)[tag==0]
  
  for(k in 1 : nrow(data_GATA6KO))
  {
    cat("\r",i,"----->",k)
    xx1 <- as.numeric(data_GATA6KO[k, sam1])
    xx2 <- as.numeric(data_GATA6KO[k, sam2])
    
    res[k,"C1_avg_expression"] <- mean(xx1)
    res[k,"C2_avg_expression"] <- mean(xx2)
    res[k,"C1_median_expression"]<-median(xx1)
    res[k,"C2_median_expression"] <- median(xx2)
  }
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.3/scRNA_bulk/GATA6KO_sig_",tar[i],"_vs_others_normalized_expression.xls")
  write.table(res,myoutf,sep="\t",quote=F)
}




data_LPMWT<-data[LPMWT_sig,]
data_LPMWT<-na.omit(data_LPMWT)

tar <- unique(info$integrated_snn_res.0.3)
#tar <- as.numeric(tar)
tar <- tar[order(tar)]

for(i in 1 : length(tar))
{
  res <- matrix(0, nrow(data_LPMWT),4)
  row.names(res) <- row.names(data_LPMWT)
  colnames(res) <- c("C1_avg_expression","C1_median_expression","C2_avg_expression","C2_median_expression")
  
  tag <- info$integrated_snn_res.0.3 == tar[i]
  
  sam1 <- row.names(info)[tag==1]
  sam2 <- row.names(info)[tag==0]
  
  for(k in 1 : nrow(data_LPMWT))
  {
    cat("\r",i,"----->",k)
    xx1 <- as.numeric(data_LPMWT[k, sam1])
    xx2 <- as.numeric(data_LPMWT[k, sam2])
    
    res[k,"C1_avg_expression"] <- mean(xx1)
    res[k,"C2_avg_expression"] <- mean(xx2)
    res[k,"C1_median_expression"]<-median(xx1)
    res[k,"C2_median_expression"] <- median(xx2)
  }
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.3/scRNA_bulk/data_LPMWT_sig_",tar[i],"_vs_others_normalized_expression.xls")
  write.table(res,myoutf,sep="\t",quote=F)
}






for (i in 1:length(tar))
{
  myinf<-paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.3/scRNA_bulk/GATA6KO_sig_",tar[i],"_vs_others_normalized_expression.xls")
  xx<-read.table(myinf,sep='\t')
  res[,i]<-xx$C1_avg_expression
}

res <- apply(res,1, function(arg) (arg-mean(arg))/sd(arg))
res <- t(res)


#res<-res[,c("2",'5','4','8','9','0','1','7','3','6')]

#ts<-res[sort(row.names(res)),]
breaklist_1<-seq(-2.5,0,0.2)
breaklist_2<-seq(0,2.6,0.2)
breaklist<-unique(c(breaklist_1,breaklist_2))

myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.3/scRNA_bulk/heatmaps/GATA6KO_signature.pdf")
pdf(myoutf,width=15,height=30)
pheatmap(mat=res,
         color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(length(breaklist)),
         breaks = breaklist,
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         border_color      = NA,
         show_colnames     = T,
         show_rownames     = T,
         drop_levels       = TRUE,
         scale='none',
         
         fontsize=20)
dev.off()



for (i in 1:length(tar))
{
  myinf<-paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.3/scRNA_bulk/data_LPMWT_sig_",tar[i],"_vs_others_normalized_expression.xls")
  xx<-read.table(myinf,sep='\t')
  res[,i]<-xx$C1_avg_expression
}

res <- apply(res,1, function(arg) (arg-mean(arg))/sd(arg))
res <- t(res)


#res<-res[,c("2",'5','4','8','9','0','1','7','3','6')]

#ts<-res[sort(row.names(res)),]
breaklist_1<-seq(-2.5,0,0.2)
breaklist_2<-seq(0,2.6,0.2)
breaklist<-unique(c(breaklist_1,breaklist_2))

myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.3/scRNA_bulk/heatmaps/LPMWT_signature.pdf")
pdf(myoutf,width=15,height=30)
pheatmap(mat=res,
         color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(length(breaklist)),
         breaks = breaklist,
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         border_color      = NA,
         show_colnames     = T,
         show_rownames     = T,
         drop_levels       = TRUE,
         scale='none',
         
         fontsize=20)
dev.off()












################heatmap single cell Seurat.
myoutf1 <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.3/cluster_identity/heatmap_marker_genes_different_phenotype_Tcell.pdf")


Idents(Tcell_PW_Adult_kit_Human.integrated)<-"integrated_snn_res.0.3"
 
tmp.markers <- FindAllMarkers(object = Tcell_PW_Adult_kit_Human.integrated, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.25,return.thresh=0.05)


top50 <- tmp.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)


pdf(myoutf1,35,35)
print(DoHeatmap(Tcell_PW_Adult_kit_Human.integrated, features = top50$gene) + 
        scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')), 
                              mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')), 
                              midpoint = 0, guide = "colourbar", aesthetics = "fill"))
dev.off()

tar<-c('0','1','3','6','7',"2",'5','4','8','9')
all_cell<-character()
for(i in 1:length(tar))
{
  cat("\r",i)
  len<-sum(Tcell_PW_Adult_kit_Human.integrated@meta.data$integrated_snn_res.0.3==tar[i])
  tag<-sample(1:len,242)
  info<-Tcell_PW_Adult_kit_Human.integrated@meta.data[Tcell_PW_Adult_kit_Human.integrated@meta.data$integrated_snn_res.0.3==tar[i],]
  barcode_use<-row.names(info)[tag]
  all_cell<-combine(all_cell,barcode_use)
}



myoutf1 <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.3/cluster_identity/heatmap_marker_genes_ordered_242cells_each_cluster.pdf")


Idents(Tcell_PW_Adult_kit_Human.integrated)<-"integrated_snn_res.0.3"

tmp.markers <- FindAllMarkers(object = Tcell_PW_Adult_kit_Human.integrated, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.15,return.thresh=0.05)

ts<-tmp.markers[order(factor(tmp.markers$cluster, levels = tar)),]

top30 <- ts %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)

Tcell_PW_Adult_kit_Human.integrated@meta.data[,"ordered_ident"]<-factor(Tcell_PW_Adult_kit_Human.integrated@meta.data$integrated_snn_res.0.3,levels=c('0','1','3','6','7',"2",'5','4','8','9'))

pdf(myoutf1,12,10)
print(DoHeatmap(Tcell_PW_Adult_kit_Human.integrated, cells=all_cell,features = top30$gene,group.by="ordered_ident") + 
        scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')), 
                              mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')), 
                              midpoint = 0, guide = "colourbar", aesthetics = "fill"))
dev.off()

myoutf2 <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.3/cluster_identity/heatmap_marker_top30.xls")

write.table(top30,myoutf2,sep='\t')







tmp.markers <- FindAllMarkers(object = Tcell_PW_Adult_kit_Human.integrated,only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.3/cluster_identity/markers_each_cluster.xls")
write.table(tmp.markers,myoutf,sep="\t",quote=F)

tmp.markers <- FindAllMarkers(object = Tcell_PW_Adult_kit_Human.integrated,only.pos = T, logfc.threshold = 0.01)
myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.3/cluster_identity/markers_each_cluster_0.01.xls")
write.table(tmp.markers,myoutf,sep="\t",quote=F)


myinf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.3/cluster_identity/markers_each_cluster_0.01.xls")
markers_all<-read.table(myinf,sep="\t")
marker0<-subset(markers_all,cluster=="0")$gene
marker1<-subset(markers_all,cluster=="1")$gene

###Do venn diagram show overlap of differentially expressed genes between two Trm clusters

myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Tcell_PW_Adult_kit_Human_Human/PCx/res.0.3/cluster_identity/Venn_Two_Trm_Overlap.pdf")
pdf(myoutf,width=10,height=10)
draw.pairwise.venn(area1 = 109, area2 = 81,  cross.area = 45,category = c("C1-FOS","C2-IFNG"),
                 lwd = rep(2, 2), lty = rep("solid", 2), col =
                   c("#54278f","#fc4e2a"), fill = c("#54278f","#fc4e2a"), alpha = rep(0.5, 2))

dev.off()



myoutf <- paste0("/dartfs-hpc/rc/lab/T/TurkLab/Jichang/bulkTCR/Fred/tumor635.pdf")
pdf(myoutf,width=10,height=10)
draw.pairwise.venn(area1 = 2240, area2 = 2901,  cross.area = 1738,category = c("RNAseq","TCRseq"),
                   lwd = rep(2, 2), lty = rep("solid", 2), col =
                     c("#54278f","#fc4e2a"), fill = c("#54278f","#fc4e2a"), alpha = rep(0.5, 2))

dev.off()




















































































#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##################Subset macrophages for downstream analysis#########################
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

screen -r 128617

export LSF_DOCKER_VOLUMES='/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/:/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/'
bsub -Is -G compute-gjrandolph -q general-interactive -R 'rusage[mem=128GB]' -a 'docker(jichanghan/scrna_seurat4_monocle3_seuratwrapper_1003:latest)' /bin/bash

library(Seurat)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)   
library(sctransform)
library(harmony)
library(metap)
library(multtest)
library(biomaRt)
library(monocle3)

rm(list=ls()) 



myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Myeloid_PW_Adult_kit_Human/Seruat_Myeloid_PW_Adult_kit.integrated__finished.RDS"
Myeloid_PW_Adult_kit.integrated<-readRDS(myinf)
tar<-"integrated_snn_res.0.5"
 Idents(object = Myeloid_PW_Adult_kit.integrated) <- tar
 subcluster<-c('0','2','3','5','6','7','15','13')
 tag<-Myeloid_PW_Adult_kit.integrated@meta.data[,tar]%in%subcluster
 tag1<-GetAssayData(Myeloid_PW_Adult_kit.integrated,assay='RNA', slot='data')['CD14',]>0.1|GetAssayData(Myeloid_PW_Adult_kit.integrated,assay='RNA', slot='data')['CD1C',]<0.5
 cell_ID<-row.names(Myeloid_PW_Adult_kit.integrated@meta.data)[tag&tag1]
 
 myinf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/Seruat_PW_PW_Adult_kit_Human__raw.RDS"
 xx<-readRDS(myinf)

 
 Macrophage_PW_Adult_kit.integrated<-subset(x=xx,cells=cell_ID)
table(Macrophage_PW_Adult_kit.integrated@meta.data$orig.ident)
  
   #Kid1   Kid10   Kid11    Kid2    Kid4    Kid6    Kid7    Kid8    Kid9 PW_Adu1 
   #1128    1195    1556     744     578     802    2136     155    1491    2580 
#PW_Adu2    PW_C    PW_Q    PW_S    PW_U    PW_W 
#   2132    4249    1687    2939     753    2809 

 Macrophage_PW_Adult_kit.integrated.list <- SplitObject(Macrophage_PW_Adult_kit.integrated, split.by = "orig.ident")


 Macrophage_PW_Adult_kit.integrated.list <- lapply(X = Macrophage_PW_Adult_kit.integrated.list, FUN = function(x) {
    x <- NormalizeData(x,assay='RNA')
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 5000,assay='RNA')
 })


############################ Using RPCA to anchor############ 
 features <- SelectIntegrationFeatures(object.list = Macrophage_PW_Adult_kit.integrated.list,nfeatures = 4000)
 Macrophage_PW_Adult_kit.integrated.list <- lapply(X = Macrophage_PW_Adult_kit.integrated.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = T)
    x <- RunPCA(x, features = features, verbose = T)
 })

 immune.anchors <- FindIntegrationAnchors(object.list = Macrophage_PW_Adult_kit.integrated.list, anchor.features = features, reduction = "rpca")
 Macrophage_PW_Adult_kit.integrated <- IntegrateData(anchorset = immune.anchors)
 
 DefaultAssay(Macrophage_PW_Adult_kit.integrated) <- "integrated"



 Macrophage_PW_Adult_kit.integrated <- ScaleData(Macrophage_PW_Adult_kit.integrated, verbose = T)
  Macrophage_PW_Adult_kit.integrated <- RunPCA(Macrophage_PW_Adult_kit.integrated, npcs = 30, verbose = T)


 x<-c(1:25)


 eigs <- Macrophage_PW_Adult_kit.integrated[["pca"]]@stdev^2 #PCA的STBandard dev的平方#
  print(paste0("Variance captured by 35 PCs: ",sum(eigs[x] / sum(eigs)))) #0.954#




# Macrophage_PW_Adult_kit.integrated do clustering and TSNE plot No regresSTBissue ---------------------------------------------


 #[5]Do clustering
 Macrophage_PW_Adult_kit.integrated<-FindNeighbors(Macrophage_PW_Adult_kit.integrated, dims = x)
  Macrophage_PW_Adult_kit.integrated <- FindClusters(Macrophage_PW_Adult_kit.integrated,resolution = seq(0.1,1.0,0.1))

PrintFindClustersParams(object = Macrophage_PW_Adult_kit.integrated)

   clusters_resolution<-sapply(grep("^integrated_snn_res",colnames(Macrophage_PW_Adult_kit.integrated@meta.data),value = TRUE),
                            function(x) length(unique(Macrophage_PW_Adult_kit.integrated@meta.data[,x])))
 clusters_resolution


 #PC=x Macrophage_PW_Adult_kit.integrated
#integrated_snn_res.0.1 integrated_snn_res.0.2 integrated_snn_res.0.3 
#                     4                      7                      9 
#integrated_snn_res.0.4 integrated_snn_res.0.5 integrated_snn_res.0.6 
#                    12                     14                     14 
#integrated_snn_res.0.7 integrated_snn_res.0.8 integrated_snn_res.0.9 
#                    14                     14                     16 
#  integrated_snn_res.1 
#                    19 


 #Set up a loop to look for markers

 tar <- c("integrated_snn_res.0.1","integrated_snn_res.0.2","integrated_snn_res.0.3","integrated_snn_res.0.4",
         "integrated_snn_res.0.5","integrated_snn_res.0.6","integrated_snn_res.0.7","integrated_snn_res.0.8",
         "integrated_snn_res.0.9","integrated_snn_res.1")




 Macrophage_PW_Adult_kit.integrated <- RunUMAP(Macrophage_PW_Adult_kit.integrated, dims= x)

 
 for(i in 1:length(tar))
 {
  Idents(object = Macrophage_PW_Adult_kit.integrated) <- tar[i]
  myoutf1 <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/UMAP_SB6_",tar[i],"_undefined.pdf")
  pdf(myoutf1,7,5)
  print(DimPlot(Macrophage_PW_Adult_kit.integrated, reduction = "umap",pt.size=0.2))
  dev.off()
 }



 myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/TSNE_MFMo_PW_Adult_kit.integrated_groupby_orig.ident.pdf"
  pdf(myoutf,width=8,height=5)
 DimPlot(object = Macrophage_PW_Adult_kit.integrated,group.by='orig.ident',reduction="umap",pt.size=0.2)
 dev.off()



 myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/Seruat_MFMo_PW_Adult_kit.integrated__finished.RDS"
 saveRDS(Macrophage_PW_Adult_kit.integrated, file = myoutf)



DefaultAssay(Macrophage_PW_Adult_kit.integrated) <- "integrated"
tar1 <- c("integrated_snn_res.0.3","integrated_snn_res.0.4","integrated_snn_res.0.5","integrated_snn_res.0.2","integrated_snn_res.0.6","integrated_snn_res.0.1")

 for(i in 1 : length(tar1))
 {
  cat("\r",i)
  Idents(object = Macrophage_PW_Adult_kit.integrated) <- tar1[i]
  tmp <- Macrophage_PW_Adult_kit.integrated
  tmp.markers <- FindAllMarkers(object = tmp, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.2,assay='integrated')
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/cluster_markers_",tar1[i],"_ScRNA.xls")
  write.table(tmp.markers,myoutf,sep="\t",quote=F)
  
  cell.list <- WhichCells(tmp, downsample = 200)
  myoutf1 <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/heatmap_markers_",tar1[i],"_.pdf")
  pdf(myoutf1,30,35)
  top20 <- tmp.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
  print(DoHeatmap(tmp, cells=cell.list,features = top20$gene,assay='integrated') + NoLegend())
  dev.off()
 }



#########################Find DEG of each cluster using findconserve markers#########
 tar<-"integrated_snn_res.0.3"
 Idents(object = Macrophage_PW_Adult_kit.integrated) <- tar
 tar1<-unique(as.numeric(Macrophage_PW_Adult_kit.integrated@meta.data$"integrated_snn_res.0.3"))
 tar1<-sort(tar1-1)
 DEG<-data.frame()
 for (i in 1:length(tar1))
 {
  markers <- FindConservedMarkers(Macrophage_PW_Adult_kit.integrated, ident.1 = tar1[i], grouping.var = "orig.ident", verbose = T,only.pos = T,logfc.threshold = 0.25,min.cells.group=1)
  markers[,"cluster"]<-tar1[i]
  markers[,'gene']<-row.names(markers)
  tar2<-grep('log2FC',colnames(markers))
  markers[,'average_Log2FC']<-apply(markers[,colnames(markers)[tar2]],1,mean)
  tag<-markers[,'max_pval']<=0.05
  markers<-markers[tag,]
  markers<-markers[order(markers[,'average_Log2FC'],decreasing=TRUE),]
  DEG<-rbind(DEG,markers)
 }





screen -r 127707

myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/Seruat_MFMo_PW_Adult_kit.integrated__finished.RDS"
Macrophage_PW_Adult_kit.integrated<-readRDS(myinf)

 DefaultAssay(Macrophage_PW_Adult_kit.integrated) <- "RNA"

data_RNA<-Macrophage_PW_Adult_kit.integrated[['RNA']]@data
data_integrated<-Macrophage_PW_Adult_kit.integrated[['integrated']]@data
data_RNA<-data_RNA[Immune_general,]
data_integrated<-data_integrated[Immune_general,]
cells<-sample(colnames(data_RNA),20)


 DefaultAssay(Macrophage_PW_Adult_kit.integrated) <- "RNA"

Immune_general<-c('SELL','IRF4','F5','CSF1R','VCAN','MERTK','MRC1','CD14','S100A8','GATA6','CD1C','LYVE1','CCR2','TIMD4','SELP','MARCO')


p<-FeaturePlot(object = Macrophage_PW_Adult_kit.integrated, features = Immune_general,cols= c('#deebf7', "#bd0026"),
               reduction = "umap",slot = "data",max.cutoff=4,pt.size=0.3,combine=F)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + theme(plot.title = element_text(size=40),
                          axis.title = element_text(family="Helvetica",face = "bold", color = "black",size=20), 
                           axis.text=element_text(family="Helvetica",face = "bold", color = "black",size=16),
                           axis.line = element_line(colour = 'black', size = 1),
                           axis.ticks = element_line(colour = "black", size = 1),
                           axis.ticks.length=unit(.3, "cm"))
}
 
myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/Myeloid_marker.pdf"
pdf(myoutf,width=50,height=50)
print(cowplot::plot_grid(plotlist = p,ncol=4))
dev.off()



#############Cluster contribution by each different patient##########
   
   Idents(Macrophage_PW_Adult_kit.integrated)<-"integrated_snn_res.0.7"
 
 xx<-table(Macrophage_PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.3,Macrophage_PW_Adult_kit.integrated@meta.data$orig.ident)
 yy<-t(apply(xx,1,function(x){x/sum(x)})   )
 myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.2/cluster_patient_distribution/patient_distribution_to_each_cluster.xls"
 write.table(yy,myoutf,sep="\t",quote=F)
 
 data<-matrix(0,ncol(yy)*nrow(yy),3)
 colnames(data)<-c("cluster","patient","value")

 data<-as.data.frame(data)
 data$cluster<-rep(row.names(yy),ncol(yy))

 data$patient<-as.character(sapply(colnames(yy),function(x){rep(x,nrow(yy))}))

 for (i in 1:nrow(data))
 {
  number<-yy[data$cluster[i],data$patient[i]]
  data$value[i]<-number
 }


 data$cluster<-factor(data$cluster,levels=unique(data$cluster))
 data$patient<-factor(data$patient,levels=c("PW_normal","PW_CRC1","PW_CRC2",'Kid1','Kid2','Kid4','Kid6','Kid7','Kid8','Kid9','Kid10'))

  col<-c("#2171b5","#cb181d","#74c476")

 myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.3/cluster_patient_distribution/stackplot_patient_contribution_each_cluster_proportion.pdf"
 pdf(myoutf,width=8,height=8)
 ggplot() + geom_bar(aes(y = value, x = cluster, fill = patient), data = data,
                    stat="identity")+theme(axis.text.x = element_text(angle = 60, hjust = 1,size=10),
                                           axis.title.x = element_blank())

 dev.off()






myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/Seruat_MFMo_PW_Adult_kit.integrated__finished.RDS"
Macrophage_PW_Adult_kit.integrated<-readRDS(myinf)

tar<-unique(Macrophage_PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.7)
tar<-tar[order(tar)]
chi<-as.matrix(table(Macrophage_PW_Adult_kit.integrated@meta.data$orig.ident,Macrophage_PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.7))

chi_square<-chisq.test(chi)

chi_observed<-chi_square$observed
chi_expected<-chi_square$expected
ROE<-chi_observed/chi_expected

myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.7/cluster_patient_distribution/ROE.xls"
write.table(ROE,myoutf,sep="\t",quote=F)


each_cluster_all<-apply(chi,2,sum)
chi_enrichment<-chi
for (i in 1:nrow(chi))
{
    for (k in 1:ncol(chi))
    {
      chi_enrichment[i,k]<-chi[i,k]*sum(chi)/sum(chi[i,])/sum(chi[,k])
    }
}

myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.7/cluster_patient_distribution/hyper_enrichment.xls"
write.table(chi_enrichment,myoutf,sep="\t",quote=F)




#################Output cell of each orig.ident#############
myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/Seruat_MFMo_PW_Adult_kit.integrated__finished.RDS"
Macrophage_PW_Adult_kit.integrated<-readRDS(myinf)

Idents(Macrophage_PW_Adult_kit.integrated)<-"integrated_snn_res.0.7"
 tar<-unique(Macrophage_PW_Adult_kit.integrated@meta.data$orig.ident)
 for (i in 1:length(tar))
 {
    tag<-Macrophage_PW_Adult_kit.integrated@meta.data$orig.ident==tar[i]
    cell_ID<-row.names(Macrophage_PW_Adult_kit.integrated@meta.data)[tag]
    myoutf1 <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.7/each_sample/UMAP_SB6_",tar[i],"_undefined.pdf")
  pdf(myoutf1,7,5)
  print(DimPlot(Macrophage_PW_Adult_kit.integrated, reduction = "umap",pt.size=0.1,cells =cell_ID))
  dev.off()
 

 }
  
  










################Cell cycle scoring###################
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/Seruat_MFMo_PW_Adult_kit.integrated__finished.RDS"
Macrophage_PW_Adult_kit.integrated<-readRDS(myinf)

DefaultAssay(Macrophage_PW_Adult_kit.integrated) <- "RNA"
Macrophage_PW_Adult_kit.integrated<-CellCycleScoring(Macrophage_PW_Adult_kit.integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
cell_cycle_score<-t(Macrophage_PW_Adult_kit.integrated@meta.data[,c("S.Score",'G2M.Score')])
adt_assay <- CreateAssayObject(counts = cell_cycle_score)
Macrophage_PW_Adult_kit.integrated[['cell_cycle']]<-adt_assay
DefaultAssay(Macrophage_PW_Adult_kit.integrated) <- "cell_cycle"
p<-FeaturePlot(object = Macrophage_PW_Adult_kit.integrated, features = c('S.Score','G2M.Score'),cols= c('#deebf7', "#bd0026"),
               reduction = "umap",slot = "data",max.cutoff=1.5,pt.size=0.3,combine=F)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + theme(plot.title = element_text(size=40),
                          axis.title = element_text(family="Helvetica",face = "bold", color = "black",size=20), 
                           axis.text=element_text(family="Helvetica",face = "bold", color = "black",size=16),
                           axis.line = element_line(colour = 'black', size = 1),
                           axis.ticks = element_line(colour = "black", size = 1),
                           axis.ticks.length=unit(.3, "cm"))
}
 
myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.3/cell_cycle/Cell_Cycle_Score.pdf"
pdf(myoutf,width=25,height=10)
print(cowplot::plot_grid(plotlist = p,ncol=2))
dev.off()

myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.3/cell_cycle/By_cell_cycle_MFMo_PW_Adult_kit.integrated.pdf"
  pdf(myoutf,width=8,height=5)
 DimPlot(object = Macrophage_PW_Adult_kit.integrated,group.by='Phase',reduction="umap",pt.size=0.3)
 dev.off()

myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/Seruat_MFMo_PW_Adult_kit.integrated__finished.RDS"
 saveRDS(Macrophage_PW_Adult_kit.integrated, file = myoutf)

Idents(Macrophage_PW_Adult_kit.integrated)<-"integrated_snn_res.0.3"
 
myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.3/cell_cycle/Cell_cycle_violin.pdf"
  pdf(myoutf,width=13,height=5)
 VlnPlot(object = Macrophage_PW_Adult_kit.integrated,assay='cell_cycle',features = c('S.Score','G2M.Score'),pt.size=0.3)
 dev.off()




Myeloid_marker<-c('XCR1','IL3RA','CLEC4C','SIGLEC6','ITGAM','ITGAX','MRC1','CD14','CD68','GATA6','CD1C','LYVE1','CCR2','TIMD4','CLEC9A')


p<-FeaturePlot(object = Macrophage_PW_Adult_kit.integrated, features = Myeloid_marker,cols= c('#deebf7', "#bd0026"),
               reduction = "umap",slot = "data",max.cutoff=4,pt.size=0.3,combine=F)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + theme(plot.title = element_text(size=40),
                          axis.title = element_text(family="Helvetica",face = "bold", color = "black",size=20), 
                           axis.text=element_text(family="Helvetica",face = "bold", color = "black",size=16),
                           axis.line = element_line(colour = 'black', size = 1),
                           axis.ticks = element_line(colour = "black", size = 1),
                           axis.ticks.length=unit(.3, "cm"))
}
 
myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/Myeloid_marker.pdf"
pdf(myoutf,width=50,height=30)
print(cowplot::plot_grid(plotlist = p,ncol=5))
dev.off()
















##########First get average gene expression of each cluster and then do z transformation to do heatmap for selected genes###########

myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/Seruat_MFMo_PW_Adult_kit.integrated__finished.RDS"
exp<-readRDS(myinf)
DefaultAssay(exp) <- "integrated"

myoutf1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.7/expression_data/CD8_normalized_data.xls"
myoutf2 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.7/expression_data/CD8_normalized_cell_info.xls"


info <- exp@meta.data
res <- GetAssayData(object = exp)
res <- as.matrix(res)
row.names(res) <- toupper(row.names(res))
avg <- apply(res,1, mean)
tag <- avg >0
res <- res[tag,]

res <- unique(res)

tag <- duplicated(row.names(res))

dup_gene <- row.names(res)[tag]

tag <- row.names(res) == dup_gene
#tag==logical(0) Thus ignore the res1,res STBep
#res1 <- res[tag==0,]

#res2 <- res[tag==1,]
#avg <- apply(res2,1, mean)

#res <- rbind(res1, res2[2,])
res <- as.data.frame(res)

write.table(res,myoutf1,sep="\t",quote=F)
write.table(info,myoutf2,sep="\t",quote=F)






myinf1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.7/expression_data/CD8_normalized_data.xls"
myinf2 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.7/expression_data/CD8_normalized_cell_info.xls"

data<-read.table(myinf1,sep="\t",quote=NULL)
info<-read.table(myinf2,sep="\t",quote=NULL)

DC_markers<-c('TCF4','IL3RA','CLEC4C','LILRB4','FCER1A','FCER1G','LILRA4',
            'IRF8','IRF4','ID2','BATF3','CLEC9A','XCR1','BTLA',
            'CD163','S100A9','S100A8','VCAN','FCN1','CD14','BST1','CD36',
            'FCGR2B','CD1C','HLA-DQA2','HLA-DQA1','HLA-DQB1','HLA-DQP1',
            'ZEB2','KLF4','ITGAX','ITGAM','CD2','SIRPA','CLEC4A','CLEC10A',
            'AXL','SIGLEC6','CX3CR1','SIGLEC1')

PM_markers<-c('GARNL3','CAV1','HGF','LAPTM5','GATA6','CD9','TGFB2','APOC1','F5','HDC','LYZ1',
               'FSCN1','ID1','PLXDC2','SELP','CXCL13','FABP5','ADGRE1','GAS6','ID2',
               'MRC1','FOLR2','LYVE1','MARCO','CD163','CD209','CD74','FPR1',
               'ANXA5','CBR2','CLEC4A1','MAF','ITGB1','CXCL13',
               'FABP7','CCR2','CSF1R','SELL','ITGAX','HLA-DQA2','HLA-DQA1','HLA-DQB1','HLA-DQP1','ITGAM',
               'RCN3','FCGRT','CCL24','IRF7','EMP1','LGALS1','HAL','RPS18')


#cluster_markers<-c('S100A9','S100A8','VCAN','FCN1','CD14','MAFB','ITGAM','SPP1',
#                    'FCER1A','CD1C','CD1E','HLA-DPB1','HLA-DQB1','HLA-DRA','CLEC10A','FCGR2B','CLEC4A','TCF4','IL3RA','ITGAX','ZEB2',
#                    'LILRB4','CSF1R','GSN','GPR183','RGS1','IL1B','FOS','NR4A1',
#                    'C1QC','C1QA','C1QB','FN1','LYVE1','FOLR2','MARCO','CD163','FCGR3A','TIMD4','APOE','CD68',
#                    'IL32','CCL2','GZMA','CCL5','NKG7','GNLY','ID3','STAT1','CD36','CCL20','ITGA4','CD3G','CD8A',
#                    'CCL4','CXCL10','CCL3','CXCL9','CXCL11',
#                    'CXCR6','GATA3','RORA','TIGIT','CD69',
#                    'MKI67','TK1','CCNA2','TUBA1B',
#                    'IRF8','ID2','BATF3','CLEC9A','XCR1','BTLA',
#                     'NEAT1','CIITA','NFKBIZ','FOXP1','TLR2','IRF4'
#                    )

data_DC<-data[DC_markers,]
data_DC<-na.omit(data_DC)

tar <- unique(info$integrated_snn_res.0.7)
#tar <- as.numeric(tar)
tar <- tar[order(tar)]

for(i in 1 : length(tar))
{
  res <- matrix(0, nrow(data_DC),4)
  row.names(res) <- row.names(data_DC)
  colnames(res) <- c("C1_avg_expression","C1_median_expression","C2_avg_expression","C2_median_expression")
  
  tag <- info$integrated_snn_res.0.7 == tar[i]
  
  sam1 <- row.names(info)[tag==1]
  sam2 <- row.names(info)[tag==0]
  
  for(k in 1 : nrow(data_DC))
  {
    cat("\r",i,"----->",k)
    xx1 <- as.numeric(data_DC[k, sam1])
    xx2 <- as.numeric(data_DC[k, sam2])
    
    res[k,"C1_avg_expression"] <- mean(xx1)
    res[k,"C2_avg_expression"] <- mean(xx2)
    res[k,"C1_median_expression"]<-median(xx1)
    res[k,"C2_median_expression"] <- median(xx2)
  }
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.7/cluster_identity/pheatmap/Human_DC_gene_",tar[i],"_vs_others_normalized_expression.xls")
  write.table(res,myoutf,sep="\t",quote=F)
}




res<-matrix(0,nrow(data_DC),length(tar))
 row.names(res)<-row.names(data_DC)
 colnames(res)<-tar[order(tar)]



for (i in 1:length(tar))
{
  myinf<-paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.7/cluster_identity/pheatmap/Human_DC_gene_",tar[i],"_vs_others_normalized_expression.xls")
  xx<-read.table(myinf,sep='\t')
  res[,i]<-xx$C1_avg_expression
}

res <- apply(res,1, function(arg) (arg-mean(arg))/sd(arg))
res <- t(res)


#res<-res[,c("2",'5','4','8','9','0','1','7','3','6')]

#ts<-res[sort(row.names(res)),]
breaklist_1<-seq(-2.5,0,0.2)
breaklist_2<-seq(0,2.6,0.2)
breaklist<-unique(c(breaklist_1,breaklist_2))

myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.7/cluster_identity/pheatmap/heatmap/Human_DC_gene_each_cluster.pdf")
pdf(myoutf,width=25,height=15)
pheatmap(mat=res,
         color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(length(breaklist)),
         breaks = breaklist,
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         border_color      = NA,
         show_colnames     = T,
         show_rownames     = T,
         drop_levels       = TRUE,
         scale='none',
         
         fontsize=20)
dev.off()












data_PM<-data[PM_markers,]
data_PM<-na.omit(data_PM)

tar <- unique(info$integrated_snn_res.0.7)
#tar <- as.numeric(tar)
tar <- tar[order(tar)]

for(i in 1 : length(tar))
{
  res <- matrix(0, nrow(data_PM),4)
  row.names(res) <- row.names(data_PM)
  colnames(res) <- c("C1_avg_expression","C1_median_expression","C2_avg_expression","C2_median_expression")
  
  tag <- info$integrated_snn_res.0.7 == tar[i]
  
  sam1 <- row.names(info)[tag==1]
  sam2 <- row.names(info)[tag==0]
  
  for(k in 1 : nrow(data_PM))
  {
    cat("\r",i,"----->",k)
    xx1 <- as.numeric(data_PM[k, sam1])
    xx2 <- as.numeric(data_PM[k, sam2])
    
    res[k,"C1_avg_expression"] <- mean(xx1)
    res[k,"C2_avg_expression"] <- mean(xx2)
    res[k,"C1_median_expression"]<-median(xx1)
    res[k,"C2_median_expression"] <- median(xx2)
  }
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.7/cluster_identity/pheatmap/Mouse_Peritoneal_MAC_gene_",tar[i],"_vs_others_normalized_expression.xls")
  write.table(res,myoutf,sep="\t",quote=F)
}




res<-matrix(0,nrow(data_PM),length(tar))
 row.names(res)<-row.names(data_PM)
 colnames(res)<-tar[order(tar)]



for (i in 1:length(tar))
{
  myinf<-paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.7/cluster_identity/pheatmap/Mouse_Peritoneal_MAC_gene_",tar[i],"_vs_others_normalized_expression.xls")
  xx<-read.table(myinf,sep='\t')
  res[,i]<-xx$C1_avg_expression
}

res <- apply(res,1, function(arg) (arg-mean(arg))/sd(arg))
res <- t(res)


#res<-res[,c("2",'5','4','8','9','0','1','7','3','6')]

#ts<-res[sort(row.names(res)),]
breaklist_1<-seq(-2.5,0,0.2)
breaklist_2<-seq(0,2.6,0.2)
breaklist<-unique(c(breaklist_1,breaklist_2))

myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.7/cluster_identity/pheatmap/heatmap/Mouse_Peritoneal_MAC_gene_each_cluster.pdf")
pdf(myoutf,width=15,height=15)
pheatmap(mat=res,
         color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(length(breaklist)),
         breaks = breaklist,
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         border_color      = NA,
         show_colnames     = T,
         show_rownames     = T,
         drop_levels       = TRUE,
         scale='none',
         
         fontsize=20)
dev.off()










data_each_cluster<-data[cluster_markers,]
data_each_cluster<-na.omit(data_each_cluster)

tar <- unique(info$integrated_snn_res.0.3)
#tar <- as.numeric(tar)
tar <- tar[order(tar)]

for(i in 1 : length(tar))
{
  res <- matrix(0, nrow(data_each_cluster),4)
  row.names(res) <- row.names(data_each_cluster)
  colnames(res) <- c("C1_avg_expression","C1_median_expression","C2_avg_expression","C2_median_expression")
  
  tag <- info$integrated_snn_res.0.3 == tar[i]
  
  sam1 <- row.names(info)[tag==1]
  sam2 <- row.names(info)[tag==0]
  
  for(k in 1 : nrow(data_each_cluster))
  {
    cat("\r",i,"----->",k)
    xx1 <- as.numeric(data_each_cluster[k, sam1])
    xx2 <- as.numeric(data_each_cluster[k, sam2])
    
    res[k,"C1_avg_expression"] <- mean(xx1)
    res[k,"C2_avg_expression"] <- mean(xx2)
    res[k,"C1_median_expression"]<-median(xx1)
    res[k,"C2_median_expression"] <- median(xx2)
  }
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.3/cluster_identity/pheatmap/marker_each_cluster_",tar[i],"_vs_others_normalized_expression.xls")
  write.table(res,myoutf,sep="\t",quote=F)
}




res<-matrix(0,nrow(data_each_cluster),length(tar))
 row.names(res)<-row.names(data_each_cluster)
 colnames(res)<-tar[order(tar)]



for (i in 1:length(tar))
{
  myinf<-paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.3/cluster_identity/pheatmap/marker_each_cluster_",tar[i],"_vs_others_normalized_expression.xls")
  xx<-read.table(myinf,sep='\t')
  res[,i]<-xx$C1_avg_expression
}

res <- apply(res,1, function(arg) (arg-mean(arg))/sd(arg))
res <- t(res)


#res<-res[,c("2",'5','4','8','9','0','1','7','3','6')]

#ts<-res[sort(row.names(res)),]
breaklist_1<-seq(-2.5,0,0.2)
breaklist_2<-seq(0,2.6,0.2)
breaklist<-unique(c(breaklist_1,breaklist_2))

myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.3/cluster_identity/pheatmap/heatmap/marker_each_cluster.pdf")
pdf(myoutf,width=8,height=17)
pheatmap(mat=res,
         color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(length(breaklist)),
         breaks = breaklist,
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         border_color      = NA,
         show_colnames     = T,
         show_rownames     = T,
         drop_levels       = TRUE,
         scale='none',
         
         fontsize=14)
dev.off()
















PM_markers<-c('GARNL3','CAV1','HGF','LAPTM5','GATA6','CD9','TGFB2','APOC1','F5','HDC','LYZ1',
               'FSCN1','ID1','PLXDC2','SELP','CXCL13','FABP5','ADGRE1','GAS6','ID2',
               'MRC1','FOLR2','LYVE1','MARCO','CD163','CD209','CD74','FPR1',
               'ANXA5','CBR2','CLEC4A1','MAF','ITGB1','PF4',
               'FABP7','CCR2','CSF1R','SELL','ITGAX','HLA-DQA2','HLA-DQA1','HLA-DQB1','HLA-DQP1','ITGAM',
               'RCN3','FCGRT','CCL24','IRF7','EMP1','LGALS1','HAL','RPS18')



PM_markers<-c('GARNL3','CAV1','HGF','LAPTM5','GATA6','CD9','TGFB2','APOC1','F5','HDC','LYZ1',
               'FSCN1','ID1','PLXDC2','SELP','CXCL13','FABP5','ADGRE1','GAS6','ID2',
               'MRC1','FOLR2','LYVE1','MARCO','CD163','CD209','CD74','FPR1',
               'ANXA5','CBR2','CLEC4A1','MAF','ITGB1','CXCL13',
               'FABP7','CCR2','CSF1R','SELL','ITGAX','HLA-DQA2','HLA-DQA1','HLA-DQB1','HLA-DQP1','ITGAM',
               'RCN3','FCGRT','CCL24','IRF7','EMP1','LGALS1','HAL','RPS18')



myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.3/cluster_identity/Dotplot/Dot_plot_MAC_heterogeneity.pdf"
pdf(myoutf,width=15,height=6)
print(DotPlot(Macrophage_PW_Adult_kit.integrated, cols=c("#08519c","#ef3b2c"),group.by="integrated_snn_res.0.3",features = PM_markers)+
        theme(axis.text.x = element_text(angle = 45,size=10,hjust = 1)))

dev.off()


DC_markers<-c('TCF4','IL3RA','CLEC4C','LILRB4','FCER1A','FCER1G','LILRA4',
            'IRF8','IRF4','ID2','BATF3','CLEC9A','XCR1','BTLA',
            'CD163','S100A9','S100A8','VCAN','FCN1','CD14','BST1','CD36',
            'FCGR2B','CD1C','HLA-DQA2','HLA-DQA1','HLA-DQB1','HLA-DQP1',
            'ZEB2','KLF4','ITGAX','ITGAM','CD2','SIRPA','CLEC4A','CLEC10A',
            'AXL','SIGLEC6','CX3CR1','SIGLEC1')


myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.3/cluster_identity/Dotplot/Dot_plot_DC_heterogeneity.pdf"
pdf(myoutf,width=15,height=6)
print(DotPlot(Macrophage_PW_Adult_kit.integrated, cols=c("#08519c","#ef3b2c"),group.by="integrated_snn_res.0.3",features = DC_markers)+
        theme(axis.text.x = element_text(angle = 45,size=10,hjust = 1)))

dev.off()





##########Check GATA6################
GATA6_exp<-as.matrix(GetAssayData(Macrophage_PW_Adult_kit.integrated,assay='RNA',slot='data')[c('GATA6','SELP'),])
apply(GATA6_exp,1,summary)

#              GATA6        SELP
#Min.    0.000000000 0.000000000
#1st Qu. 0.000000000 0.000000000
#Median  0.000000000 0.000000000
#Mean    0.009995213 0.007311954
#3rd Qu. 0.000000000 0.000000000
#Max.    2.931832980 2.733619776

tag1<-GATA6_exp[1,]>mean(GATA6_exp[1,])
tag2<-GATA6_exp[2,]>mean(GATA6_exp[2,])

GATA6_pos<-colnames(GATA6_exp)[tag1]
SELP_pos<-colnames(GATA6_exp)[tag2]

Macrophage_PW_Adult_kit.integrated@meta.data[,'GATA6']<-rep("Neg",nrow(Macrophage_PW_Adult_kit.integrated@meta.data))
Macrophage_PW_Adult_kit.integrated@meta.data[,'SELP']<-rep("Neg",nrow(Macrophage_PW_Adult_kit.integrated@meta.data))
Macrophage_PW_Adult_kit.integrated@meta.data[GATA6_pos,'GATA6']<-'Pos'
Macrophage_PW_Adult_kit.integrated@meta.data[SELP_pos,'SELP']<-'Pos'

GATA6_dis<-table(Macrophage_PW_Adult_kit.integrated@meta.data[,'GATA6'],Macrophage_PW_Adult_kit.integrated@meta.data[,'orig.ident'])
SELP_dis<-table(Macrophage_PW_Adult_kit.integrated@meta.data[,'SELP'],Macrophage_PW_Adult_kit.integrated@meta.data[,'orig.ident'])



chi_square<-chisq.test(GATA6_dis)

chi_observed<-chi_square$observed
chi_expected<-chi_square$expected
ROE_GATA6<-chi_observed/chi_expected


chi_square<-chisq.test(SELP_dis)

chi_observed<-chi_square$observed
chi_expected<-chi_square$expected
ROE_SELP<-chi_observed/chi_expected



GATA6_enrichment<-GATA6_dis
for (i in 1:nrow(GATA6_dis))
{
    for (k in 1:ncol(GATA6_dis))
    {
      GATA6_enrichment[i,k]<-GATA6_dis[i,k]*sum(GATA6_dis)/sum(GATA6_dis[i,])/sum(GATA6_dis[,k])
    }
}



myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.7/cluster_patient_distribution/GATA6_enrichment.xls"
write.table(ROE_GATA6,myoutf,sep="\t",quote=F)

myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.7/cluster_patient_distribution/SELP_enrichment.xls"
write.table(ROE_SELP,myoutf,sep="\t",quote=F)





##############Find doublets in the dataset#############
 
 myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/Seruat_MFMo_PW_Adult_kit.integrated__finished.RDS"
 Macrophage_PW_Adult_kit.integrated<-readRDS(myinf)
 cell_ID<-row.names(Macrophage_PW_Adult_kit.integrated@meta.data)

 myinf1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/raw/Nan_PW_Normal"
 PW_normal.data <- Read10X(data.dir=myinf1)
 PW_normal <- CreateSeuratObject(counts = PW_normal.data, project = "PW_normal",min.cells=3,min.features=200)
 PW_normal<-RenameCells(object=PW_normal,add.cell.id=unique(PW_normal@meta.data$orig.ident))
 PW_normal<-RenameCells(object=PW_normal,new.names=gsub('-1','',Cells(PW_normal)))


  myinf2 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/raw/Nan_PW_CRC1/sample_feature_bc_matrix"
 PW_CRC1.data <- Read10X(data.dir=myinf2)
 PW_CRC1 <- CreateSeuratObject(counts = PW_CRC1.data, project = "PW_CRC1",min.cells=3,min.features=200)
 PW_CRC1<-RenameCells(object=PW_CRC1,add.cell.id=unique(PW_CRC1@meta.data$orig.ident))
 PW_CRC1<-RenameCells(object=PW_CRC1,new.names=gsub('-1','',Cells(PW_CRC1)))


 myinf3 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/raw/Nan_PW_CRC2/sample_feature_bc_matrix"
 PW_CRC2.data <- Read10X(data.dir=myinf3)
 PW_CRC2 <- CreateSeuratObject(counts = PW_CRC2.data, project = "PW_CRC2",min.cells=3,min.features=200)
 PW_CRC2<-RenameCells(object=PW_CRC2,add.cell.id=unique(PW_CRC2@meta.data$orig.ident))
 PW_CRC2<-RenameCells(object=PW_CRC2,new.names=gsub('-1','',Cells(PW_CRC2)))

 xx<- merge(PW_normal, y = c(PW_CRC1, PW_CRC2), project = "PW_Adult_kit")


 Myeloid_PW_Adult_kit.doublet<-subset(x=xx,cells=cell_ID)
 
 Myeloid_PW_Adult_kit.doublet.list <- SplitObject(Myeloid_PW_Adult_kit.doublet, split.by = "orig.ident")

 Myeloid_PW_Adult_kit.doublet.list <- lapply(X = Myeloid_PW_Adult_kit.doublet.list, FUN = function(x) {
    x <- NormalizeData(x,assay='RNA')
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000,assay='RNA')
    x <- ScaleData(x, verbose = T)
    x <- RunPCA(x,verbose = T)
    x<-FindNeighbors(x, dims = 1:15)
   x <- FindClusters(x,resolution = 0.3)
 
 })









 Myeloid_PW_Adult_kit.doublet.list <- lapply(X = Myeloid_PW_Adult_kit.doublet.list, FUN = function(x) {
    yy <- paramSweep_v3(x, PCs = 1:15, sct = FALSE)
    zz<-summarizeSweep(yy, GT = FALSE)
    tt<-find.pK(zz)
    pK<-tt %>% filter(BCmetric == max(BCmetric)) %>% select (pK)
    pK<- as.numeric (as.character(pK[[1]]))
    homotypic.prop <- modelHomotypic(x@meta.data$"integrated_snn_res.0.3")
    nExp_poi <- round(0.055*nrow(x@meta.data)) 
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    x <- doubletFinder_v3(x, PCs = 1:15, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
 })


 Myeloid_PW_Adult_kit.doublet<-merge(
  Myeloid_PW_Adult_kit.doublet.list[[1]],
  y =c(Myeloid_PW_Adult_kit.doublet.list[[2]], Myeloid_PW_Adult_kit.doublet.list[[3]]),
  add.cell.ids = NULL,
  project = "Myeloid_PW_Adult_kit.doublet"
)

tag<-Myeloid_PW_Adult_kit.doublet@meta.data[,"DF.classifications_0.25_0.17_448"]=='Doublet'
doublet_ID<-row.names(Myeloid_PW_Adult_kit.doublet@meta.data[tag,])

Macrophage_PW_Adult_kit.integrated@meta.data[,"doublet"]<-rep('Singlet',nrow(Macrophage_PW_Adult_kit.integrated@meta.data))
Macrophage_PW_Adult_kit.integrated@meta.data[doublet_ID,"doublet"]<-"Doublet"

myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/doublet_finder.pdf"
  pdf(myoutf,width=8,height=5)
 DimPlot(object = Macrophage_PW_Adult_kit.integrated,group.by='doublet',reduction="umap",pt.size=0.3,cells.highlight=doublet_ID,cols.highlight = "#e6550d",
  sizes.highlight = 0.3)
 dev.off()














####################Perform GSEA########

myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/Seruat_MFMo_PW_Adult_kit.integrated__finished.RDS"
Macrophage_PW_Adult_kit.integrated<-readRDS(myinf)
DefaultAssay(Macrophage_PW_Adult_kit.integrated) <- "RNA"
tar<-"integrated_snn_res.0.3"
 Idents(object = Macrophage_PW_Adult_kit.integrated) <- tar
 tar1<-unique(as.numeric(Macrophage_PW_Adult_kit.integrated@meta.data$"integrated_snn_res.0.3"))
 tar1<-sort(tar1-1)
 for (i in 1:length(tar1))
 {
  markers <- FindConservedMarkers(Macrophage_PW_Adult_kit.integrated, ident.1 = tar1[i], grouping.var = "orig.ident", verbose = T,logfc.threshold = 0.001,min.cells.group=1,min.pct=0.001)
  markers[,"cluster"]<-tar1[i]
  markers[,'gene']<-row.names(markers)
  tar2<-grep('log2FC',colnames(markers))
  markers[,'average_Log2FC']<-apply(markers[,colnames(markers)[tar2]],1,mean)
  markers<-markers[order(markers[,'average_Log2FC'],decreasing=TRUE),]
  myoutf<-paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.3/GSEA/marker/Cluster_",tar1[i],"differentially_expressed_markers.xls")
  write.table(markers,myoutf,quote=F,sep="\t")
 }


dir <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/PM_KA_manuscript/Human_bulkRNA/analysis/marker/"
files <- list.files(dir)
file_nam <- gsub(".xls","",files)
for(i in 1:length(files))
{
  cat("\r",i)
  tmpinf <- paste0(dir,files[i])
  res <- read.table(tmpinf,sep="\t",quote=NULL)
  res <- res[order(res$Log2FC,decreasing=T),]
  
  xx <- res$Log2FC
  names(xx) <- row.names(res)
  names(xx) <- toupper(names(xx))
  
  tag = duplicated(names(xx))
  res = res[tag==0,]
  row.names(res) <- toupper(row.names(res))
  res[,"Name"] <- row.names(res)
  res[,"metric"] <- res$Log2FC
  res <- res[,c("Name","metric")]
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/PM_KA_manuscript/Human_bulkRNA/analysis/GSEA/",file_nam[i],".rnk")
  write.table(res,file=myoutf,quote=F,sep="\t",row.names=F)
}





dir <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.3/GSEA/marker/"
files <- list.files(dir)
file_nam <- gsub(".xls","",files)
for(i in 1:length(files))
{
  cat("\r",i)
  tmpinf <- paste0(dir,files[i])
  res <- read.table(tmpinf,sep="\t",quote=NULL)
  #tag<-res$average_Log2FC>0
  #res<-res[tag,]
  res <- res[order(res$average_Log2FC,decreasing=T),]
  
  xx <- res$average_Log2FC
  names(xx) <- row.names(res)
  names(xx) <- toupper(names(xx))
  
  tag = duplicated(names(xx))
  res = res[tag==0,]
  row.names(res) <- toupper(row.names(res))
  res[,"Name"] <- row.names(res)
  res[,"metric"] <- res$average_Log2FC
  res <- res[,c("Name","metric")]
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.3/GSEA/preranked_GSEA/",file_nam[i],".rnk")
  write.table(res,file=myoutf,quote=F,sep="\t",row.names=F)
}























##########First get average gene expression of each cluster and then do z transformation to do heatmap for selected genes###########

myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/Seruat_MFMo_PW_Adult_kit.integrated__finished.RDS"
exp <- readRDS(myinf)
DefaultAssay(exp) <- "integrated"

myoutf1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.3/expression_data/CD8_normalized_data.xls"
myoutf2 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.3/expression_data/CD8_normalized_cell_info.xls"


info <- exp@meta.data
res <- GetAssayData(object = exp)
res <- as.matrix(res)
row.names(res) <- toupper(row.names(res))
avg <- apply(res,1, mean)
tag <- avg >0
res <- res[tag,]

res <- unique(res)

tag <- duplicated(row.names(res))

dup_gene <- row.names(res)[tag]

tag <- row.names(res) == dup_gene
#tag==logical(0) Thus ignore the res1,res STBep
#res1 <- res[tag==0,]

#res2 <- res[tag==1,]
#avg <- apply(res2,1, mean)

#res <- rbind(res1, res2[2,])
res <- as.data.frame(res)

write.table(res,myoutf1,sep="\t",quote=F)
write.table(info,myoutf2,sep="\t",quote=F)

myinf1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.3/expression_data/CD8_normalized_data.xls"
myinf2 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.3/expression_data/CD8_normalized_cell_info.xls"

data<-read.table(myinf1,sep="\t",quote=NULL)
info<-read.table(myinf2,sep="\t",quote=NULL)


dir <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/PM_KA_manuscript/Human_bulkRNA/analysis/marker/"
files <- list.files(dir)
file_nam <- gsub(".txt","",files)
for(i in 1:length(files))
{
  cat("\r",i)
  tmpinf <- paste0(dir,files[i])
  res1 <- read.table(tmpinf,sep="\t",quote=NULL)
  #tag<-res$average_Log2FC>0
  #res<-res[tag,]
  res1 <- res1[order(res1$Log2FC,decreasing=T),]
  FC2_sig<-row.names(res1[res1$Log2FC>2,])
  FC1_sig<-row.names(res1[res1$Log2FC>1,])
  top100_sig<-row.names(res1[1:100,])
  
  data_FC2<-data[row.names(data)%in%FC2_sig,]
  data_FC1<-data[row.names(data)%in%FC1_sig,]
  data_top100<-data[row.names(data)%in%top100_sig,]
  tar <- unique(info$integrated_snn_res.0.3)
 #tar <- as.numeric(tar)
 tar <- tar[order(tar)]

 for(m in 1 : length(tar))
 {
  res <- matrix(0, nrow(data_FC2),4)
  row.names(res) <- row.names(data_FC2)
  colnames(res) <- c("C1_avg_expression","C1_median_expression","C2_avg_expression","C2_median_expression")
  
  tag <- info$integrated_snn_res.0.3 == tar[m]
  
  sam1 <- row.names(info)[tag==1]
  sam2 <- row.names(info)[tag==0]
  
  for(k in 1 : nrow(data_FC2))
  {
    cat("\r",m,"----->",k)
    xx1 <- as.numeric(data_FC2[k, sam1])
    xx2 <- as.numeric(data_FC2[k, sam2])
    
    res[k,"C1_avg_expression"] <- mean(xx1)
    res[k,"C2_avg_expression"] <- mean(xx2)
    res[k,"C1_median_expression"]<-median(xx1)
    res[k,"C2_median_expression"] <- median(xx2)
  }
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.3/scRNA_bulk/",file_nam[i],"_FC2_genes_cluster_",tar[m],"_vs_others_normalized_expression.xls")
  write.table(res,myoutf,sep="\t",quote=F)
 }

 for(m in 1 : length(tar))
 {
  res <- matrix(0, nrow(data_FC1),4)
  row.names(res) <- row.names(data_FC1)
  colnames(res) <- c("C1_avg_expression","C1_median_expression","C2_avg_expression","C2_median_expression")
  
  tag <- info$integrated_snn_res.0.3 == tar[m]
  
  sam1 <- row.names(info)[tag==1]
  sam2 <- row.names(info)[tag==0]
  
  for(k in 1 : nrow(data_FC1))
  {
    cat("\r",m,"----->",k)
    xx1 <- as.numeric(data_FC1[k, sam1])
    xx2 <- as.numeric(data_FC1[k, sam2])
    
    res[k,"C1_avg_expression"] <- mean(xx1)
    res[k,"C2_avg_expression"] <- mean(xx2)
    res[k,"C1_median_expression"]<-median(xx1)
    res[k,"C2_median_expression"] <- median(xx2)
  }
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.3/scRNA_bulk/",file_nam[i],"_FC1_genes_cluster_",tar[m],"_vs_others_normalized_expression.xls")
  write.table(res,myoutf,sep="\t",quote=F)
}

for(m in 1 : length(tar))
 {
  res <- matrix(0, nrow(data_top100),4)
  row.names(res) <- row.names(data_top100)
  colnames(res) <- c("C1_avg_expression","C1_median_expression","C2_avg_expression","C2_median_expression")
  
  tag <- info$integrated_snn_res.0.3 == tar[m]
  
  sam1 <- row.names(info)[tag==1]
  sam2 <- row.names(info)[tag==0]
  
  for(k in 1 : nrow(data_top100))
  {
    cat("\r",m,"----->",k)
    xx1 <- as.numeric(data_top100[k, sam1])
    xx2 <- as.numeric(data_top100[k, sam2])
    
    res[k,"C1_avg_expression"] <- mean(xx1)
    res[k,"C2_avg_expression"] <- mean(xx2)
    res[k,"C1_median_expression"]<-median(xx1)
    res[k,"C2_median_expression"] <- median(xx2)
  }
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.3/scRNA_bulk/",file_nam[i],"_top100_genes_cluster_",tar[m],"_vs_others_normalized_expression.xls")
  write.table(res,myoutf,sep="\t",quote=F)
}

}



dir <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/PM_KA_manuscript/Human_bulkRNA/analysis/marker/"
files <- list.files(dir)
file_nam <- gsub(".txt","",files)
for (k in 1:length(file_nam))
{
  myinf<-paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.3/scRNA_bulk/",file_nam[k],"_FC2_genes_cluster_0_vs_others_normalized_expression.xls")
  xx<-read.table(myinf,sep='\t')
  
  res<-matrix(0,nrow(xx),length(tar))
 row.names(res)<-row.names(xx)
 colnames(res)<-tar[order(tar)]
  for (i in 1:length(tar))
{
  myinf<-paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.3/scRNA_bulk/",file_nam[k],"_FC2_genes_cluster_",tar[i],"_vs_others_normalized_expression.xls")
  xx<-read.table(myinf,sep='\t')
  res[,i]<-xx$C1_avg_expression
}
res <- apply(res,1, function(arg) (arg-mean(arg))/sd(arg))
res <- t(res)


#res<-res[,c("2",'5','4','8','9','0','1','7','3','6')]

#ts<-res[sort(row.names(res)),]
breaklist_1<-seq(-2.5,0,0.2)
breaklist_2<-seq(0,2.6,0.2)
breaklist<-unique(c(breaklist_1,breaklist_2))

myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.3/scRNA_bulk/heatmaps/",file_nam[k],"_FC2_heatmap_everycluster.pdf")
pdf(myoutf,width=15,height=30)
pheatmap(mat=res,
         color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(length(breaklist)),
         breaks = breaklist,
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         border_color      = NA,
         show_colnames     = T,
         show_rownames     = T,
         drop_levels       = TRUE,
         scale='none',
         
         fontsize=20)
dev.off()


}






myinf1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.3/expression_data/CD8_normalized_data.xls"
myinf2 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.3/expression_data/CD8_normalized_cell_info.xls"

data<-read.table(myinf1,sep="\t",quote=NULL)
info<-read.table(myinf2,sep="\t",quote=NULL)

myinf3<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/PM_KA_manuscript/Mouse_PM/Meta_data/GATA6KO_vs_WT_DEG.txt"
myinf4<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/PM_KA_manuscript/Mouse_PM/Meta_data/GATA6KO_vs_WT_DEG.txt"



data3<-read.table(myinf3, sep="\t", header=T)
data4<-read.table(myinf4, sep="\t", header=T)

GATA6KO_sig<-toupper(data3[data3$Log2FC>1,'X'])
LPMWT_sig<-toupper(data3[data3$Log2FC<=-1,'X'])


##Mouse_to_human <- function(x){

   #human <- useEnsembl(biomart = "ensembl", 
   #                dataset = "hsapiens_gene_ensembl", 
   #                mirror = "uswest")
   #mouse <- useEnsembl(biomart = "ensembl", 
   #                dataset = "mmusculus_gene_ensembl", 
   #                mirror = "uswest")


#human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
#mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

#   genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
#
 #  humanx <- unique(genesV2[, 2])

 #  return(humanx)
#}

data_GATA6KO<-data[GATA6KO_sig,]
data_GATA6KO<-na.omit(data_GATA6KO)

tar <- unique(info$integrated_snn_res.0.3)
#tar <- as.numeric(tar)
tar <- tar[order(tar)]

for(i in 1 : length(tar))
{
  res <- matrix(0, nrow(data_GATA6KO),4)
  row.names(res) <- row.names(data_GATA6KO)
  colnames(res) <- c("C1_avg_expression","C1_median_expression","C2_avg_expression","C2_median_expression")
  
  tag <- info$integrated_snn_res.0.3 == tar[i]
  
  sam1 <- row.names(info)[tag==1]
  sam2 <- row.names(info)[tag==0]
  
  for(k in 1 : nrow(data_GATA6KO))
  {
    cat("\r",i,"----->",k)
    xx1 <- as.numeric(data_GATA6KO[k, sam1])
    xx2 <- as.numeric(data_GATA6KO[k, sam2])
    
    res[k,"C1_avg_expression"] <- mean(xx1)
    res[k,"C2_avg_expression"] <- mean(xx2)
    res[k,"C1_median_expression"]<-median(xx1)
    res[k,"C2_median_expression"] <- median(xx2)
  }
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.3/scRNA_bulk/GATA6KO_sig_",tar[i],"_vs_others_normalized_expression.xls")
  write.table(res,myoutf,sep="\t",quote=F)
}




data_LPMWT<-data[LPMWT_sig,]
data_LPMWT<-na.omit(data_LPMWT)

tar <- unique(info$integrated_snn_res.0.3)
#tar <- as.numeric(tar)
tar <- tar[order(tar)]

for(i in 1 : length(tar))
{
  res <- matrix(0, nrow(data_LPMWT),4)
  row.names(res) <- row.names(data_LPMWT)
  colnames(res) <- c("C1_avg_expression","C1_median_expression","C2_avg_expression","C2_median_expression")
  
  tag <- info$integrated_snn_res.0.3 == tar[i]
  
  sam1 <- row.names(info)[tag==1]
  sam2 <- row.names(info)[tag==0]
  
  for(k in 1 : nrow(data_LPMWT))
  {
    cat("\r",i,"----->",k)
    xx1 <- as.numeric(data_LPMWT[k, sam1])
    xx2 <- as.numeric(data_LPMWT[k, sam2])
    
    res[k,"C1_avg_expression"] <- mean(xx1)
    res[k,"C2_avg_expression"] <- mean(xx2)
    res[k,"C1_median_expression"]<-median(xx1)
    res[k,"C2_median_expression"] <- median(xx2)
  }
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.3/scRNA_bulk/data_LPMWT_sig_",tar[i],"_vs_others_normalized_expression.xls")
  write.table(res,myoutf,sep="\t",quote=F)
}






for (i in 1:length(tar))
{
  myinf<-paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.3/scRNA_bulk/GATA6KO_sig_",tar[i],"_vs_others_normalized_expression.xls")
  xx<-read.table(myinf,sep='\t')
  res[,i]<-xx$C1_avg_expression
}

res <- apply(res,1, function(arg) (arg-mean(arg))/sd(arg))
res <- t(res)


#res<-res[,c("2",'5','4','8','9','0','1','7','3','6')]

#ts<-res[sort(row.names(res)),]
breaklist_1<-seq(-2.5,0,0.2)
breaklist_2<-seq(0,2.6,0.2)
breaklist<-unique(c(breaklist_1,breaklist_2))

myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.3/scRNA_bulk/heatmaps/GATA6KO_signature.pdf")
pdf(myoutf,width=15,height=30)
pheatmap(mat=res,
         color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(length(breaklist)),
         breaks = breaklist,
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         border_color      = NA,
         show_colnames     = T,
         show_rownames     = T,
         drop_levels       = TRUE,
         scale='none',
         
         fontsize=20)
dev.off()



for (i in 1:length(tar))
{
  myinf<-paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.3/scRNA_bulk/data_LPMWT_sig_",tar[i],"_vs_others_normalized_expression.xls")
  xx<-read.table(myinf,sep='\t')
  res[,i]<-xx$C1_avg_expression
}

res <- apply(res,1, function(arg) (arg-mean(arg))/sd(arg))
res <- t(res)


#res<-res[,c("2",'5','4','8','9','0','1','7','3','6')]

#ts<-res[sort(row.names(res)),]
breaklist_1<-seq(-2.5,0,0.2)
breaklist_2<-seq(0,2.6,0.2)
breaklist<-unique(c(breaklist_1,breaklist_2))

myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.3/scRNA_bulk/heatmaps/LPMWT_signature.pdf")
pdf(myoutf,width=15,height=30)
pheatmap(mat=res,
         color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(length(breaklist)),
         breaks = breaklist,
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         border_color      = NA,
         show_colnames     = T,
         show_rownames     = T,
         drop_levels       = TRUE,
         scale='none',
         
         fontsize=20)
dev.off()












################heatmap single cell Seurat.
myoutf1 <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.3/cluster_identity/heatmap_marker_genes_different_phenotype_Tcell.pdf")


Idents(Macrophage_PW_Adult_kit.integrated)<-"integrated_snn_res.0.3"
 
tmp.markers <- FindAllMarkers(object = Macrophage_PW_Adult_kit.integrated, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.25,return.thresh=0.05)


top50 <- tmp.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)


pdf(myoutf1,35,35)
print(DoHeatmap(Macrophage_PW_Adult_kit.integrated, features = top50$gene) + 
        scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')), 
                              mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')), 
                              midpoint = 0, guide = "colourbar", aesthetics = "fill"))
dev.off()

tar<-c('0','1','3','6','7',"2",'5','4','8','9')
all_cell<-character()
for(i in 1:length(tar))
{
  cat("\r",i)
  len<-sum(Macrophage_PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.3==tar[i])
  tag<-sample(1:len,242)
  info<-Macrophage_PW_Adult_kit.integrated@meta.data[Macrophage_PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.3==tar[i],]
  barcode_use<-row.names(info)[tag]
  all_cell<-combine(all_cell,barcode_use)
}



myoutf1 <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.3/cluster_identity/heatmap_marker_genes_ordered_242cells_each_cluster.pdf")


Idents(Macrophage_PW_Adult_kit.integrated)<-"integrated_snn_res.0.3"

tmp.markers <- FindAllMarkers(object = Macrophage_PW_Adult_kit.integrated, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.15,return.thresh=0.05)

ts<-tmp.markers[order(factor(tmp.markers$cluster, levels = tar)),]

top30 <- ts %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)

Macrophage_PW_Adult_kit.integrated@meta.data[,"ordered_ident"]<-factor(Macrophage_PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.3,levels=c('0','1','3','6','7',"2",'5','4','8','9'))

pdf(myoutf1,12,10)
print(DoHeatmap(Macrophage_PW_Adult_kit.integrated, cells=all_cell,features = top30$gene,group.by="ordered_ident") + 
        scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')), 
                              mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')), 
                              midpoint = 0, guide = "colourbar", aesthetics = "fill"))
dev.off()

myoutf2 <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.3/cluster_identity/heatmap_marker_top30.xls")

write.table(top30,myoutf2,sep='\t')







tmp.markers <- FindAllMarkers(object = Macrophage_PW_Adult_kit.integrated,only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.3/cluster_identity/markers_each_cluster.xls")
write.table(tmp.markers,myoutf,sep="\t",quote=F)

tmp.markers <- FindAllMarkers(object = Macrophage_PW_Adult_kit.integrated,only.pos = T, logfc.threshold = 0.01)
myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.3/cluster_identity/markers_each_cluster_0.01.xls")
write.table(tmp.markers,myoutf,sep="\t",quote=F)


myinf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.3/cluster_identity/markers_each_cluster_0.01.xls")
markers_all<-read.table(myinf,sep="\t")
marker0<-subset(markers_all,cluster=="0")$gene
marker1<-subset(markers_all,cluster=="1")$gene

###Do venn diagram show overlap of differentially expressed genes between two Trm clusters

myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/Macrophage_PW_Adult_kit_Human/PCx/res.0.3/cluster_identity/Venn_Two_Trm_Overlap.pdf")
pdf(myoutf,width=10,height=10)
draw.pairwise.venn(area1 = 109, area2 = 81,  cross.area = 45,category = c("C1-FOS","C2-IFNG"),
                 lwd = rep(2, 2), lty = rep("solid", 2), col =
                   c("#54278f","#fc4e2a"), fill = c("#54278f","#fc4e2a"), alpha = rep(0.5, 2))

dev.off()



myoutf <- paste0("/dartfs-hpc/rc/lab/T/TurkLab/Jichang/bulkTCR/Fred/tumor635.pdf")
pdf(myoutf,width=10,height=10)
draw.pairwise.venn(area1 = 2240, area2 = 2901,  cross.area = 1738,category = c("RNAseq","TCRseq"),
                   lwd = rep(2, 2), lty = rep("solid", 2), col =
                     c("#54278f","#fc4e2a"), fill = c("#54278f","#fc4e2a"), alpha = rep(0.5, 2))

dev.off()


































######Tree cluster##############


Idents(object = PW_Adult_kit.integrated) <- 'integrated_snn_res.0.3'


PW_Adult_kit.integrated_tree <- BuildClusterTree(
  PW_Adult_kit.integrated,
  reorder = F,
  reorder.numeric = F,
  dims=1:15)

myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/Cluster_tree_validation.pdf"
pdf(myoutf,width=5,height=5)
print(PlotClusterTree(PW_Adult_kit.integrated_tree))
dev.off()




VG<-VariableFeatures(object = PW_Adult_kit.integrated)
viriablegene<-GetAssayData(object = PW_Adult_kit.integrated)[VG,]



# Preranked GSEA -----------------------------------------------------------


myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/Seruat_PW_Adult_kit.integrated__finished.RDS"
PW_Adult_kit.integrated<-readRDS(myinf)
Idents(object = PW_Adult_kit.integrated) <- "integrated_snn_res.0.3"
tar<-unique(PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.3)
tar<-tar[order(tar)]
for(i in 1:length(tar))
{
  cluster.markers <- FindMarkers(PW_Adult_kit.integrated, ident.1 =tar[i],logfc.threshold = 0.01,min.pct = 0.01)
  myoutf<-paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/GSEA/marker/Cluster_",tar[i],"differentially_expressed_markers.xls")
  write.table(cluster.markers,myoutf,quote=F,sep="\t")
}



dir <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/GSEA/marker/"
files <- list.files(dir)
file_nam <- gsub(".xls","",files)
for(i in 1:length(files))
{
  cat("\r",i)
  tmpinf <- paste0(dir,files[i])
  res <- read.table(tmpinf,sep="\t",quote=NULL)
  res <- res[order(res$avg_logFC,decreasing=T),]
  
  xx <- res$avg_logFC
  names(xx) <- row.names(res)
  names(xx) <- toupper(names(xx))
  
  tag = duplicated(names(xx))
  res = res[tag==0,]
  row.names(res) <- toupper(row.names(res))
  res[,"Name"] <- row.names(res)
  res[,"metric"] <- res$avg_logFC
  res <- res[,c("Name","metric")]
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/GSEA/preranked_GSEA/",file_nam[i],".rnk")
  write.table(res,file=myoutf,quote=F,sep="\t",row.names=F)
}










myinf<-"/dartfs-hpc/rc/lab/T/TurkLab/Jichang/yanding/prerank_geneset.xls"

res <- read.table(myinf,sep="\t",quote=NULL)

res <- res[order(res$T.score,decreasing=T),]

xx <- res$T.score
names(xx) <- row.names(res)
names(xx) <- toupper(names(xx))

tag = duplicated(names(xx))
res = res[tag==0,]
row.names(res) <- toupper(row.names(res))
res[,"Name"] <- row.names(res)
res[,"metric"] <- res$T.score
res <- res[,c("Name","metric")]

myoutf <- paste0("/dartfs-hpc/rc/lab/T/TurkLab/Jichang/yanding/preranked_GSEA_Yanding.rnk")
write.table(res,file=myoutf,quote=F,sep="\t",row.names=F)










############Dot plot show Trm heterogeneity

Hetero_Trm<-c("IFNG","TNF","CCL4","CCL3","CCL5","GZMA","GZMB","CD69","CTLA4","PDCD1","TOX","HAVCR2","TIGIT","LAG3")
Trm<-c("CD69","NR4A1","RGS1","CXCR6")
IFNG<-c("CCL4","IFNG","TNF","CCL5")
JUN<-c("FOS","JUNB","FOSB","JUN")
Tcir<-c("KLF2","SELL","TCF7","S1PR1")



Idents(PW_Adult_kit.integrated)<-"integrated_snn_res.0.3"

PW_Adult_kit.integrated@meta.data[,"Dotplot"]<-rep("Other",nrow(PW_Adult_kit.integrated@meta.data))
tag1<-PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.3==0
tag2<-PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.3==1
tag3<-PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.3==7

PW_Adult_kit.integrated@meta.data$Dotplot[tag1]<-"C0-AP1-Trm"
PW_Adult_kit.integrated@meta.data$Dotplot[tag2]<-"C1-IFNG-Trm"
PW_Adult_kit.integrated@meta.data$Dotplot[tag3]<-"C7-CXCL13-Tex"


myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/Trm_signature/Dot plot Trm heterogeneity.pdf"
pdf(myoutf,width=8,height=6)
print(DotPlot(PW_Adult_kit.integrated, cols=c("#08519c","#ef3b2c"),group.by="Dotplot",features = Hetero_Trm)+
        theme(axis.text.x = element_text(angle = 45,size=10,hjust = 1)))

dev.off()

myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/Trm_signature/Dot plot Trm heterogeneity_Trm.pdf"
pdf(myoutf,width=6,height=3)
print(DotPlot(PW_Adult_kit.integrated, cols=c("#08519c","#ef3b2c"),col.min=-1,col.max=1,group.by="Dotplot",features = Trm)+
        theme(axis.text.x = element_text(angle = 45,size=10,hjust = 1)))

dev.off()


myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/Trm_signature/Dot plot Trm heterogeneity_IFNG.pdf"
pdf(myoutf,width=6,height=3)
print(DotPlot(PW_Adult_kit.integrated, cols=c("#08519c","#ef3b2c"),col.min=-1,col.max=1,group.by="Dotplot",features = IFNG)+
        theme(axis.text.x = element_text(angle = 45,size=10,hjust = 1)))

dev.off()



myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/Trm_signature/Dot plot Trm heterogeneity_Tcir.pdf"
pdf(myoutf,width=6,height=3)
print(DotPlot(PW_Adult_kit.integrated, cols=c("#08519c","#ef3b2c"),col.min=-1,col.max=1,group.by="Dotplot",features = Tcir)+
        theme(axis.text.x = element_text(angle = 45,size=10,hjust = 1)))

dev.off()








# Get normalized data and z transformed data ------------------------------

###Get normalized data and average expression for each cluster###
myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/Seruat_PW_Adult_kit.integrated__finished.RDS"
exp <- readRDS(myinf)

myoutf1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/expression_data/CD8_normalized_data.xls"
myoutf2 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/expression_data/CD8_normalized_cell_info.xls"


info <- exp@meta.data
res <- GetAssayData(object = exp)
res <- as.matrix(res)
row.names(res) <- toupper(row.names(res))
avg <- apply(res,1, mean)
tag <- avg >0
res <- res[tag,]

res <- unique(res)

tag <- duplicated(row.names(res))

dup_gene <- row.names(res)[tag]

tag <- row.names(res) == dup_gene
#tag==logical(0) Thus ignore the res1,res STBep
#res1 <- res[tag==0,]

#res2 <- res[tag==1,]
#avg <- apply(res2,1, mean)

#res <- rbind(res1, res2[2,])
res <- as.data.frame(res)

write.table(res,myoutf1,sep="\t",quote=F)
write.table(info,myoutf2,sep="\t",quote=F)















############Calculate the tissue distribution of each cluster each patient#############
#[1]chi-square

myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/Seruat_PW_Adult_kit.integrated__finished.RDS"
PW_Adult_kit.integrated<-readRDS(myinf)

tar<-unique(PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.3)
tar<-tar[order(tar)]
chi<-matrix(0,12,length(tar))
number<-sort(rep((unique(PW_Adult_kit.integrated@meta.data$orig.ident)),3))
tissue<-rep(c('skin','tumor','blood'),4)
row.names(chi)<-paste0(number,tissue)
colnames(chi)<-tar

tar2<-unique(PW_Adult_kit.integrated@meta.data$orig.ident)
tar1<-unique(PW_Adult_kit.integrated@meta.data$label)

for(j in 1:length(tar2))
{
  yy<-PW_Adult_kit.integrated@meta.data[PW_Adult_kit.integrated@meta.data$orig.ident==tar2[j],]
  
  for (i in 1:length(tar1))
  {
    tag<-yy$label==tar1[i]
    xx<-yy[tag,]
    for (k in 1:length(tar))
    {
      chi[paste0(tar2[j],tar1[i]),tar[k]]<-sum(xx$integrated_snn_res.0.3==tar[k])
    }
    
  }
  
  
}


chi_square<-chisq.test(chi)

chi_observed<-chi_square$observed
chi_expected<-chi_square$expected
ROE<-chi_observed/chi_expected

myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/cluster_tissue_distribution/ROE.xls"
write.table(ROE,myoutf,sep="\t",quote=F)




































################TCR matched clonotypes gene expression profile################


tar<-c(1,3,5,6)
barcode_info_all<-data.frame()

for(i in 1:length(tar))
{
  myinf1<-paste0("/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scTCRseq/Run_",tar[i],"/TCR/data/skin/skin_filtered_contig_annotations.csv")
  myinf2<-paste0("/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scTCRseq/Run_",tar[i],"/TCR/data/blood/blood_filtered_contig_annotations.csv")
  myinf3<-paste0("/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scTCRseq/Run_",tar[i],"/TCR/data/tumor/tumor_filtered_contig_annotations.csv")
  data_skin<-read.csv(myinf1)
  data_blood<-read.csv(myinf2)
  data_tumor<-read.csv(myinf3)
  data_skin[,"label"]<-rep("skin",nrow(data_skin))
  data_blood[,"label"]<-rep("blood",nrow(data_blood))
  data_tumor[,"label"]<-rep("tumor",nrow(data_tumor))
  
  data_skin[,"labeled_clonotype_id"]<-paste0(data_skin$raw_clonotype_id,"_",data_skin$label,tar[i])
  data_blood[,"labeled_clonotype_id"]<-paste0(data_blood$raw_clonotype_id,"_",data_blood$label,tar[i])
  data_tumor[,"labeled_clonotype_id"]<-paste0(data_tumor$raw_clonotype_id,"_",data_tumor$label,tar[i])
  
  barcode_info<-rbind(data_skin,data_blood,data_tumor)
  
  
  tag<-barcode_info$is_cell=="True"
  barcode_info<-barcode_info[tag,]
  tag<-barcode_info$high_confidence=="True"
  barcode_info<-barcode_info[tag,]
  tag<-barcode_info$productive=="True"
  barcode_info<-barcode_info[tag,]
  
  barcode_info$barcode<-gsub("-1","",barcode_info$barcode)
  barcode_info$barcode<-paste0(barcode_info$label,"_",barcode_info$barcode)
  barcode_info$barcode<-gsub("_",paste0("_",tar[i],"_"),barcode_info$barcode)
  barcode_info_all<-rbind(barcode_info,barcode_info_all)
}

barcode_info_all_brief<-unique(barcode_info_all[,c("barcode","label","labeled_clonotype_id")])







###way1
tar<-c(1,3,5,6)
clonotype_STB_skin<-character()
clonotype_STB_tumor<-character()
clonotype_ST_skin<-character()
clonotype_ST_tumor<-character()




clonotype_BT_blood_all<-character()

all_info_ST_skin<-data.frame()
all_info_blood<-data.frame()

clonotype_label_STB_all<-character()
clonotype_label_ST_all<-character()


for(i in 1:length(tar))
{
  cat('\r',i)
  myinf1<-paste0("/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scTCRseq/Run_",tar[i],"/TCR/data/skin/skin_clonotypes.csv")
  myinf2<-paste0("/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scTCRseq/Run_",tar[i],"/TCR/data/blood/blood_clonotypes.csv")
  myinf3<-paste0("/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scTCRseq/Run_",tar[i],"/TCR/data/tumor/tumor_clonotypes.csv")
  data_skin<-read.csv(myinf1)
  data_blood<-read.csv(myinf2)
  data_tumor<-read.csv(myinf3)
  data_skin[,"label"]<-rep("skin",nrow(data_skin))
  data_blood[,"label"]<-rep("blood",nrow(data_blood))
  data_tumor[,"label"]<-rep("tumor",nrow(data_tumor))
  
  data_skin[,"labeled_clonotype_id"]<-paste0(data_skin$clonotype_id,"_",data_skin$label,tar[i])
  data_blood[,"labeled_clonotype_id"]<-paste0(data_blood$clonotype_id,"_",data_blood$label,tar[i])
  data_tumor[,"labeled_clonotype_id"]<-paste0(data_tumor$clonotype_id,"_",data_tumor$label,tar[i])
  
  skin_tumor_match<-merge(data_skin,data_tumor,by.x="cdr3s_nt",by.y="cdr3s_nt")
  
  tag<-skin_tumor_match$cdr3s_nt%in%data_blood$cdr3s_nt
  STB_match<-skin_tumor_match[tag,]
  ST_specific_match<-skin_tumor_match[!tag,]
  
  clonotype_STB_skin<-combine(clonotype_STB_skin,STB_match$labeled_clonotype_id.x)
  clonotype_STB_tumor<-combine(clonotype_STB_tumor, STB_match$labeled_clonotype_id.y)
  
  
  
  clonotype_ST_skin<-combine(clonotype_ST_skin,ST_specific_match$labeled_clonotype_id.x)
  clonotype_ST_tumor<-combine(clonotype_ST_tumor, ST_specific_match$labeled_clonotype_id.y)
  
  blood_tumor_match<-merge(data_blood,data_tumor,by.x="cdr3s_nt",by.y="cdr3s_nt")
  clonotypeSB_blood
  
  all_info_STB_skin<-rbind(all_info_ST_skin,skin_tumor_match)
  
  all_info_blood<-rbind(all_info_blood,data_blood)
  
  blood_tumor_match<-merge(data_blood,data_tumor,by.x="cdr3s_nt",by.y="cdr3s_nt")
  
  clonotype_BT_blood_all<-combine(clonotype_BT_blood_all,blood_tumor_match$labeled_clonotype_id.x)
  clonotype_tumor_all<-combine( clonotype_tumor_all,blood_tumor_match$labeled_clonotype_id.y)
  
  
  
  
  
  #tag<- skin_tumor_match$cdr3s_aa%in%data_blood$cdr3s_aa
  #skin_tumor_match[,"match"]<-rep("skin tumor specific matched",nrow(skin_tumor_match))
  #skin_tumor_match$match[tag]<-"skin tumor blood matched"
  #STB<-merge(skin_tumor_match,data_blood,by="cdr3s_aa")
  #tag1<-skin_tumor_match$match=="skin tumor specific matched"
  #ST<-skin_tumor_match[tag1,]
  
  
  #clonotype_STB_all<-rbind(clonotype_STB_all,STB)
  #clonotype_ST_all<-rbind(clonotype_ST_all,ST)
  
  #clonotype_STB_label<-combine(STB$labeled_clonotype_id.x,STB$labeled_clonotype_id.y,STB$labeled_clonotype_id)
  #clonotype_ST_label<-combine(ST$labeled_clonotype_id.x,ST$labeled_clonotype_id.y)
  
  #clonotype_label_STB_all<-combine(clonotype_label_STB_all,clonotype_STB_label)
  #clonotype_label_ST_all<-combine(clonotype_label_ST_all,clonotype_ST_label)
  
}

clonotype_tumor_all<-unique(clonotype_tumor_all)

tag1<-barcode_info_all_brief$labeled_clonotype_id%in%clonotype_ST_skin_all
tag2<-barcode_info_all_brief$labeled_clonotype_id%in%clonotype_BT_blood_all
tag3<-barcode_info_all_brief$labeled_clonotype_id%in%clonotype_tumor_all

barcode_skin<-unique(barcode_info_all_brief[tag1,"barcode"])
barcode_blood<-unique(barcode_info_all_brief[tag2,"barcode"])
barcode_tumor<-unique(barcode_info_all_brief[tag3,"barcode"])


myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/Seruat_PW_Adult_kit.integrated__finished.RDS"
PW_Adult_kit.integrated<-readRDS(myinf)


tag1<-row.names(PW_Adult_kit.integrated@meta.data)%in%barcode_skin
tag2<-row.names(PW_Adult_kit.integrated@meta.data)%in%barcode_blood
tag3<-row.names(PW_Adult_kit.integrated@meta.data)%in%barcode_tumor
tag4<-PW_Adult_kit.integrated@meta.data$label=="blood"


PW_Adult_kit.integrated@meta.data[,"match"]<-rep("Not Matched",nrow(PW_Adult_kit.integrated@meta.data))

PW_Adult_kit.integrated@meta.data$match[tag1]<-"Skin match with tumor"
PW_Adult_kit.integrated@meta.data$match[tag2]<-"Blood match with tumor"
PW_Adult_kit.integrated@meta.data$match[tag3]<-"Tumor match with skin"
PW_Adult_kit.integrated@meta.data$match[tag4]<-"Blood Not Show"


tag<-PW_Adult_kit.integrated@meta.data$match=="Skin match with tumor"
xx<-PW_Adult_kit.integrated@meta.data[tag,]

#PT_1 PT_3 PT_5 PT_6 
#56   25  410  369 

myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/matched_TCR/UMAP_matched_clonotypes_show_only_skin_notshow_blood.pdf"
pdf(myoutf,width=7.5,height=5)
DimPlot(object = PW_Adult_kit.integrated,group.by='match',shape.by="label",reduction="umap",
        cols = c('white','#f0f0f0','#de2d26',"#f0f0f0"),pt.size = 0.3)
dev.off()



###Way2

tar<-c(1,3,5,6)
clonotype_STB_all_skin<-character()
clonotype_STB_all_tumor<-character()
clonotype_STB_all_blood<-character()

clonotype_ST_all_tumor<-character()
clonotype_ST_all_skin<-character()

clonotype_BT_all_tumor<-character()
clonotype_BT_all_blood<-character()



all_info_STB<-data.frame()

all_info_ST<-data.frame()

all_info_BT<-data.frame()


clonotype_label_STB_all<-character()
clonotype_label_ST_all<-character()


for(i in 1:length(tar))
{
   cat('\r',i)
   myinf1<-paste0("/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scTCRseq/Run_",tar[i],"/TCR/data/skin/skin_clonotypes.csv")
   myinf2<-paste0("/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scTCRseq/Run_",tar[i],"/TCR/data/blood/blood_clonotypes.csv")
   myinf3<-paste0("/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scTCRseq/Run_",tar[i],"/TCR/data/tumor/tumor_clonotypes.csv")
   data_skin<-read.csv(myinf1)
   data_blood<-read.csv(myinf2)
   data_tumor<-read.csv(myinf3)
   data_skin[,"label"]<-rep("skin",nrow(data_skin))
   data_blood[,"label"]<-rep("blood",nrow(data_blood))
   data_tumor[,"label"]<-rep("tumor",nrow(data_tumor))
   
   data_skin[,"labeled_clonotype_id"]<-paste0(data_skin$clonotype_id,"_",data_skin$label,tar[i])
   data_blood[,"labeled_clonotype_id"]<-paste0(data_blood$clonotype_id,"_",data_blood$label,tar[i])
   data_tumor[,"labeled_clonotype_id"]<-paste0(data_tumor$clonotype_id,"_",data_tumor$label,tar[i])
   
   ST_match<-merge(data_skin,data_tumor,by.x="cdr3s_nt",by.y="cdr3s_nt")
   STB_match<-merge( ST_match,data_blood,by="cdr3s_nt")
   clonotype_STB_all_skin<-combine(clonotype_STB_all_skin,STB_match$labeled_clonotype_id.x)
   clonotype_STB_all_tumor<-combine(clonotype_STB_all_tumor,STB_match$labeled_clonotype_id.y)
   clonotype_STB_all_blood<-combine(clonotype_STB_all_blood,STB_match$labeled_clonotype_id)
   all_info_STB<-rbind(all_info_STB,STB_match)
  
   
   tag<-ST_match$cdr3s_nt%in%data_blood$cdr3s_nt
   ST_spec_match<-ST_match[!tag,]
   clonotype_ST_all_tumor<-combine(clonotype_ST_all_tumor,ST_spec_match$labeled_clonotype_id.y)
   clonotype_ST_all_skin<-combine(clonotype_ST_all_skin,ST_spec_match$labeled_clonotype_id.x)
   all_info_ST<-rbind( all_info_ST,ST_spec_match)
   
   
   BT_match<-merge(data_blood,data_tumor,by.x="cdr3s_nt",by.y="cdr3s_nt")
   tag<-BT_match$cdr3s_nt%in%STB_match$cdr3s_nt
   BT_spec_match<-BT_match[!tag,]
   clonotype_BT_all_tumor<-combine(clonotype_BT_all_tumor,BT_spec_match$labeled_clonotype_id.y)
   clonotype_BT_all_blood<-combine(clonotype_BT_all_blood,BT_spec_match$labeled_clonotype_id.x)
   all_info_BT<-rbind(all_info_BT, BT_spec_match)

   
   #tag<- skin_tumor_match$cdr3s_aa%in%data_blood$cdr3s_aa
   #skin_tumor_match[,"match"]<-rep("skin tumor specific matched",nrow(skin_tumor_match))
   #skin_tumor_match$match[tag]<-"skin tumor blood matched"
   #STB<-merge(skin_tumor_match,data_blood,by="cdr3s_aa")
   #tag1<-skin_tumor_match$match=="skin tumor specific matched"
   #ST<-skin_tumor_match[tag1,]
   
   
   #clonotype_STB_all<-rbind(clonotype_STB_all,STB)
   #clonotype_ST_all<-rbind(clonotype_ST_all,ST)
   
   #clonotype_STB_label<-combine(STB$labeled_clonotype_id.x,STB$labeled_clonotype_id.y,STB$labeled_clonotype_id)
   #clonotype_ST_label<-combine(ST$labeled_clonotype_id.x,ST$labeled_clonotype_id.y)
   
   #clonotype_label_STB_all<-combine(clonotype_label_STB_all,clonotype_STB_label)
   #clonotype_label_ST_all<-combine(clonotype_label_ST_all,clonotype_ST_label)
   
}





tag1<-barcode_info_all_brief$labeled_clonotype_id%in%clonotype_STB_all_skin
tag2<-barcode_info_all_brief$labeled_clonotype_id%in%clonotype_STB_all_tumor
tag3<-barcode_info_all_brief$labeled_clonotype_id%in%clonotype_STB_all_blood

barcode_STB_skin<-unique(barcode_info_all_brief[tag1,"barcode"])
barcode_STB_tumor<-unique(barcode_info_all_brief[tag2,"barcode"])
barcode_STB_blood<-unique(barcode_info_all_brief[tag3,"barcode"])


tag1<-barcode_info_all_brief$labeled_clonotype_id%in%clonotype_ST_all_skin
tag2<-barcode_info_all_brief$labeled_clonotype_id%in%clonotype_ST_all_tumor

barcode_ST_skin<-unique(barcode_info_all_brief[tag1,"barcode"])
barcode_ST_tumor<-unique(barcode_info_all_brief[tag2,"barcode"])

tag1<-barcode_info_all_brief$labeled_clonotype_id%in%clonotype_BT_all_blood
tag2<-barcode_info_all_brief$labeled_clonotype_id%in%clonotype_BT_all_tumor

barcode_BT_blood<-unique(barcode_info_all_brief[tag1,"barcode"])
barcode_BT_tumor<-unique(barcode_info_all_brief[tag2,"barcode"])



myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/Seruat_PW_Adult_kit.integrated__finished.RDS"
PW_Adult_kit.integrated<-readRDS(myinf)


tag1<-row.names(PW_Adult_kit.integrated@meta.data)%in%barcode_STB_skin
tag2<-row.names(PW_Adult_kit.integrated@meta.data)%in%barcode_STB_tumor
tag3<-row.names(PW_Adult_kit.integrated@meta.data)%in%barcode_STB_blood

PW_Adult_kit.integrated@meta.data[,"STB_match"]<-rep("Not Matched",nrow(PW_Adult_kit.integrated@meta.data))

PW_Adult_kit.integrated@meta.data$STB_match[tag1]<-"STB_skin"
PW_Adult_kit.integrated@meta.data$STB_match[tag2]<-"STB_tumor"
PW_Adult_kit.integrated@meta.data$STB_match[tag3]<-"STB_blood"

table(PW_Adult_kit.integrated@meta.data$STB_match,PW_Adult_kit.integrated@meta.data$orig.ident)
#              PT_1 PT_3 PT_5 PT_6
#Not Matched 1714 2572 2585 1730
#STB_blood     83   16   57  768
#STB_skin      25    7  185  159
#STB_tumor     41    9  428  279

tag1<-row.names(PW_Adult_kit.integrated@meta.data)%in%barcode_ST_skin
tag2<-row.names(PW_Adult_kit.integrated@meta.data)%in%barcode_ST_tumor

PW_Adult_kit.integrated@meta.data[,"ST_match"]<-rep("Not Matched",nrow(PW_Adult_kit.integrated@meta.data))

PW_Adult_kit.integrated@meta.data$ST_match[tag1]<-"ST_skin"
PW_Adult_kit.integrated@meta.data$ST_match[tag2]<-"ST_tumor"

table(PW_Adult_kit.integrated@meta.data$ST_match,PW_Adult_kit.integrated@meta.data$orig.ident)
#             PT_1 PT_3 PT_5 PT_6
#Not Matched 1767 2560 2837 2575
#ST_skin       31   18  225  210
#ST_tumor      65   26  193  151



tag1<-row.names(PW_Adult_kit.integrated@meta.data)%in%barcode_BT_blood
tag2<-row.names(PW_Adult_kit.integrated@meta.data)%in%barcode_BT_tumor

PW_Adult_kit.integrated@meta.data[,"BT_match"]<-rep("Not Matched",nrow(PW_Adult_kit.integrated@meta.data))

PW_Adult_kit.integrated@meta.data$BT_match[tag1]<-"BT_blood"
PW_Adult_kit.integrated@meta.data$BT_match[tag2]<-"BT_tumor"

table(PW_Adult_kit.integrated@meta.data$BT_match,PW_Adult_kit.integrated@meta.data$orig.ident)
#             PT_1 PT_3 PT_5 PT_6
#BT_blood      24   26    6  169
#BT_tumor      23   49    6   66
#Not Matched 1816 2529 3243 2701





myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/matched_TCR/seperate_STB_ST_BT/UMAP_matched_clonotypes_STB.pdf"
pdf(myoutf,width=7.5,height=5)
DimPlot(object = PW_Adult_kit.integrated,group.by='STB_match',shape.by="orig.ident",reduction="umap",
        cols = c('#f0f0f0','#2ca25f',"#de2d26","#9ecae1"),pt.size = 0.3)
dev.off()

myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/matched_TCR/seperate_STB_ST_BT/UMAP_matched_clonotypes_ST.pdf"
pdf(myoutf,width=7.5,height=5)
DimPlot(object = PW_Adult_kit.integrated,group.by='ST_match',shape.by="orig.ident",reduction="umap",
        cols = c('#f0f0f0',"#de2d26","#9ecae1"),pt.size = 0.3)
dev.off()


myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/matched_TCR/seperate_STB_ST_BT/UMAP_matched_clonotypes_BT.pdf"
pdf(myoutf,width=7.5,height=5)
DimPlot(object = PW_Adult_kit.integrated,group.by='BT_match',shape.by="orig.ident",reduction="umap",
        cols = c('#2ca25f',"#9ecae1",'#f0f0f0'),pt.size = 0.3)
dev.off()


tar<-unique(PW_Adult_kit.integrated@meta.data$orig.ident)
for(i in 1:length(tar))
{
  cat("\t",i)
  Idents(PW_Adult_kit.integrated)<-"orig.ident"
  PW_Adult_kit.integrated_PT<-subset(PW_Adult_kit.integrated,idents=tar[i])
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/matched_TCR/seperate_STB_ST_BT/UMAP_matched_clonotypes_STB_",tar[i],"_.pdf")
  pdf(myoutf,width=7.5,height=5)
  
  print(DimPlot(object = PW_Adult_kit.integrated_PT,group.by='STB_match',shape.by="label",reduction="umap",
          cols = c('#f0f0f0','#2ca25f',"#de2d26","#9ecae1"),pt.size = 0.3))
  
  dev.off()
  
}



myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/Seruat_PW_Adult_kit.integrated__finished.RDS"
saveRDS(PW_Adult_kit.integrated,myoutf)








####Do chi-square test for the STB matched clonotypes enrichment in different clusters#########
##This part is for the STB_TB_ST seperated plots. Only focusing chi-square on STB match


table(PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.3,PW_Adult_kit.integrated@meta.data$STB_match)


myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/Seruat_PW_Adult_kit.integrated__finished.RDS"
PW_Adult_kit.integrated<-readRDS(myinf)

tar<-unique(PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.3)
tar<-tar[order(tar)]
enrich_all<-matrix(0,length(unique(PW_Adult_kit.integrated@meta.data$orig.ident))*length(tar),3)
colnames(enrich_all)<-c("Cluster","PT","enrich_index")
enrich_all[,"Cluster"]<-as.character(rep(tar,4))
enrich_all[,"PT"]<-as.character(sapply(unique(PW_Adult_kit.integrated@meta.data$orig.ident),function(x){rep(x,10)}))
enrich_all<-as.data.frame(enrich_all)
enrich_all$enrich_index<-as.numeric(levels(enrich_all$enrich_index))[enrich_all$enrich_index]

chi_all<-data.frame()
for(i in 1:length(unique(PW_Adult_kit.integrated@meta.data$orig.ident)))
{
  
  cat("\r",i)
  chi<-matrix(0,2,length(tar))
  PT<-unique(PW_Adult_kit.integrated@meta.data$orig.ident)[i]
  match<-c('matched','not matched')
  row.names(chi)<-paste0(PT,"_",match)
  colnames(chi)<-tar
  
  data<-PW_Adult_kit.integrated@meta.data[PW_Adult_kit.integrated@meta.data$orig.ident==unique(PW_Adult_kit.integrated@meta.data$orig.ident)[i],]
  
  for(j in 1:length(tar))
  {
    data_1<-data[data$integrated_snn_res.0.3==tar[j],]
    chi[1,tar[j]]<-sum(data_1$STB_match=="STB_skin"|data_1$STB_match=="STB_blood"|data_1$STB_match=="STB_tumor")
    #data_1<-data_1[data_1$label!="tumor",]
    chi[2,tar[j]]<-sum(data_1$STB_match=="Not Matched")
    
  }
  #chi_all<-rbind(chi_all,chi)
  chi_square<-chisq.test(chi)
  
  chi_observed<-chi_square$observed
  chi_expected<-chi_square$expected
  ROE<-chi_observed/chi_expected
  
  for(k in 1:length(tar))
  {
    
    enrich_all[enrich_all$PT==PT,][k,"enrich_index"]<-ROE[1,tar[k]]
  }
  
}








tag<-enrich_all$Cluster=="0"|enrich_all$Cluster=="1"
enrich_zero_one<-enrich_all[tag,]


col<-brewer.pal(n = 10, name ="Paired")

myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/matched_TCR/Boxplot_matched_TCR_enrichment_level_cluster0&1.pdf"
pdf(myoutf,width=3,height=5)
print(ggplot(enrich_zero_one, aes(x=Cluster, y=enrich_index,fill=Cluster)) + 
        geom_boxplot()+
        geom_jitter(aes(shape = PT), position=position_jitter(0))+
        theme_classic()+
        scale_fill_manual(values=c("#a1d99b","#9ecae1")))
dev.off()


cluster0<-enrich_zero_one[enrich_zero_one$Cluster=="0","enrich_index"]
cluster1<-enrich_zero_one[enrich_zero_one$Cluster=="1","enrich_index"]


t.test(cluster0,cluster1,alternative = "less",paired = T)

#t.test(cluster0,cluster1,alternative = "less")

#Welch Two Sample t-test

#data:  cluster0 and cluster1
#t = -1.7744, df = 4.0918, p-value = 0.07453
#alternative hypothesis: true difference in means is less than 0
#95 percent confidence interval:
#  -Inf 0.1047329
#sample estimates:
#  mean of x mean of y 
#0.5553939 1.0959650 





####Do percentage of matched cells distribution across all clusters#########
##This part is for the STB_TB_ST seperated plots. Only focusing chi-square on STB match



table(PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.3,PW_Adult_kit.integrated@meta.data$STB_match)


myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/Seruat_PW_Adult_kit.integrated__finished.RDS"
PW_Adult_kit.integrated<-readRDS(myinf)

tar<-unique(PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.3)
tar<-tar[order(tar)]
enrich_all<-matrix(0,length(unique(PW_Adult_kit.integrated@meta.data$orig.ident))*length(tar),3)
colnames(enrich_all)<-c("Cluster","PT","percentage")
enrich_all[,"Cluster"]<-as.character(rep(tar,4))
enrich_all[,"PT"]<-as.character(sapply(unique(PW_Adult_kit.integrated@meta.data$orig.ident),function(x){rep(x,10)}))
enrich_all<-as.data.frame(enrich_all)
enrich_all$percentage<-as.numeric(levels(enrich_all$percentage))[enrich_all$percentage]

chi_all<-data.frame()
for(i in 1:length(unique(PW_Adult_kit.integrated@meta.data$orig.ident)))
{
  
  cat("\r",i)
  chi<-matrix(0,2,length(tar))
  PT<-unique(PW_Adult_kit.integrated@meta.data$orig.ident)[i]
  match<-c('matched','not matched')
  row.names(chi)<-paste0(PT,"_",match)
  colnames(chi)<-tar
  
  data<-PW_Adult_kit.integrated@meta.data[PW_Adult_kit.integrated@meta.data$orig.ident==unique(PW_Adult_kit.integrated@meta.data$orig.ident)[i],]
  
  for(j in 1:length(tar))
  {
    data_1<-data[data$integrated_snn_res.0.3==tar[j],]
    chi[1,tar[j]]<-sum(data_1$STB_match=="STB_skin"|data_1$STB_match=="STB_blood")
    #data_1<-data_1[data_1$label!="tumor",]
    chi[2,tar[j]]<-sum(data_1$STB_match=="Not Matched")
    
  }
  
  
  percent<-chi[1,]/sum(chi[1,])
  chi<-rbind(chi,percent)
  
  
  
  
  
  for(k in 1:length(tar))
  {
    
    enrich_all[enrich_all$PT==PT,][k,"percentage"]<-chi[3,tar[k]]*100
  }
  
}








tag<-enrich_all$Cluster=="0"|enrich_all$Cluster=="1"|enrich_all$Cluster=="3"|enrich_all$Cluster=="6"|enrich_all$Cluster=="7"
enrich_ST<-enrich_all[tag,]
enrich_ST<-enrich_ST[enrich_ST$PT!="PT_3",]


col<-brewer.pal(n = 5, name ="Paired")

myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/matched_TCR/Boxplot_STB_matched_TCR_percentage_level_cluster_in_ST_NoPT3.pdf"
pdf(myoutf,width=3,height=5)
print(ggplot(enrich_ST, aes(x=Cluster, y=percentage,fill=Cluster)) + 
        geom_boxplot()+
        geom_jitter(aes(shape = PT), position=position_jitter(0))+
        theme_classic()+
        scale_fill_manual(values=col))
dev.off()


cluster0<-enrich_ST[enrich_ST$Cluster=="0","percentage"]
cluster1<-enrich_ST[enrich_ST$Cluster=="1","percentage"]


t.test(cluster0,cluster1,alternative = "less",paired = T)























#############Plot single clonotype###########
#Expanded in skin

barcode_all<-character()



tag1<-all_info_ST_skin$frequency.x>10
all_info_ST_skin_alot<-all_info_ST_skin[tag1,]

#clonotype_STB_all<-clonotype_STB_all[order(-as.numeric(clonotype_STB_all$frequency.x)),]
#clonotype_ST_all<-clonotype_ST_all[order(-as.numeric(clonotype_ST_all$frequency.x)),]

#tag2<-clonotype_ST_all$frequency.x>2&clonotype_ST_all$frequency.y>2
#clonotype_ST_all_alot<-clonotype_ST_all[tag2,]


col<-c("#f03b20","#31a354","#756bb1",'#f0f0f0')

for(i in 1:nrow(all_info_ST_skin_alot))
{
  cat("\r",i)
   tag<-barcode_info_all$labeled_clonotype_id%in%all_info_ST_skin_alot[i,c("labeled_clonotype_id.x",
                                                                             "labeled_clonotype_id.y")]
   
   
   barcode<-barcode_info_all[tag,"barcode"]
   
   tag1<-all_info_blood$cdr3s_nt%in%all_info_ST_skin_alot[i,"cdr3s_nt"]
   clonotype_blood<-all_info_blood$labeled_clonotype_id[tag1]
   tag2<-barcode_info_all$labeled_clonotype_id==clonotype_blood
   barcode_blood<-barcode_info_all[tag2,"barcode"]
   barcode<-combine(barcode,barcode_blood)
   barcode_all<-combine(barcode_all,barcode)
   
   
   PW_Adult_kit.integrated@meta.data[,"single_clonotype"]<-rep("Others",nrow(PW_Adult_kit.integrated@meta.data))
   
   tag<-row.names(PW_Adult_kit.integrated@meta.data)%in%barcode
   
   PW_Adult_kit.integrated@meta.data$single_clonotype[tag]<-paste0("Clonotype",i,PW_Adult_kit.integrated@meta.data$label[tag])
   
   
   xx<-length(unique( PW_Adult_kit.integrated@meta.data$single_clonotype))
   
   myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/matched_TCR/single_clonotype_match/clonal_expanded_in_Trm/UMAP_matched_clonotype",i,".pdf")
   pdf(myoutf,width=6.5,height=5)
   print(DimPlot(object = PW_Adult_kit.integrated,group.by='single_clonotype',reduction="umap",
           cols = col[1:xx],pt.size = 0.2))
   dev.off()
   
}


col<-brewer.pal(n = 12, name ="Paired")

for(j in 1:nrow(clonotype_ST_all_alot))
{
   cat("\r",j)
   tag<-barcode_info_all$labeled_clonotype_id%in%clonotype_ST_all_alot[j,c("labeled_clonotype_id.x",
                                                                            "labeled_clonotype_id.y")]
   
   barcode<-barcode_info_all[tag,"barcode"]
   
   PW_Adult_kit.integrated@meta.data[,"single_clonotype"]<-rep("Others",nrow(PW_Adult_kit.integrated@meta.data))
   
   tag<-row.names(PW_Adult_kit.integrated@meta.data)%in%barcode
   
   PW_Adult_kit.integrated@meta.data$single_clonotype[tag]<-paste0("Clonotype",(j+nrow(clonotype_STB_all_alot)))
   
   myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/matched_TCR/single_clonotype_match/UMAP_matched_clonotype",(j+nrow(clonotype_STB_all_alot)),".pdf")
   pdf(myoutf,width=6.5,height=5)
   print(DimPlot(object = PW_Adult_kit.integrated,group.by='single_clonotype',reduction="umap",
                 cols = c("red",'#f0f0f0'),pt.size = 0.2))
   dev.off()
   
}














#Plot single clonotypes that matches between STB and expanded in blood
#Point: IFNG Trm has functional counterpart in blood as well


barcode_STB_blood<-row.names(PW_Adult_kit.integrated@meta.data[PW_Adult_kit.integrated@meta.data$STB_match=="STB_blood",])
tag<-unique(barcode_info_all_brief)$barcode%in%barcode_STB_blood
clonotype<-barcode_info_all_brief[tag,]
clonotype[,"frequency"]<-rep(0,nrow(clonotype))
for(i in 1:length(unique(clonotype$labeled_clonotype_id)))
{
   tag<-clonotype$labeled_clonotype_id==unique(clonotype$labeled_clonotype_id)[i]
   xx<-clonotype[tag,]
   clonotype$frequency[tag]=nrow(xx)
   
}

tag<-clonotype$frequency>2
blood_expanded<-clonotype[tag,]


clonotype_to_plot<-unique(blood_expanded$labeled_clonotype_id)


col<-c('#2ca25f',"#de2d26","#9ecae1",'#f0f0f0')

for(i in 1:length(clonotype_to_plot))
{
   cat("\r",i)
   tag<-all_info_STB$labeled_clonotype_id==clonotype_to_plot[i]
   clonotype_STB<-combine(all_info_STB$labeled_clonotype_id.x[tag],
                          all_info_STB$labeled_clonotype_id.y[tag],
                          clonotype_to_plot[i])
   
   tag<-barcode_info_all_brief$labeled_clonotype_id%in%clonotype_STB
   barcode_STB<-barcode_info_all_brief$barcode[tag]
   barcode_STB_skin<-barcode_STB[grep("skin",barcode_STB)]
   barcode_STB_blood<-barcode_STB[grep("blood",barcode_STB)]
   barcode_STB_tumor<-barcode_STB[grep("tumor",barcode_STB)]
   
   
   PW_Adult_kit.integrated@meta.data[,"single_clonotype"]<-rep("Others",nrow( PW_Adult_kit.integrated@meta.data))
   tag1<-row.names(PW_Adult_kit.integrated@meta.data)%in%barcode_STB_skin
   tag2<-row.names(PW_Adult_kit.integrated@meta.data)%in%barcode_STB_blood
   tag3<-row.names(PW_Adult_kit.integrated@meta.data)%in%barcode_STB_tumor
   
   PW_Adult_kit.integrated@meta.data$single_clonotype[tag1]<-paste0("Clonotype",i,"_Skin")
   PW_Adult_kit.integrated@meta.data$single_clonotype[tag2]<-paste0("Clonotype",i,"_Blood")
   PW_Adult_kit.integrated@meta.data$single_clonotype[tag3]<-paste0("Clonotype",i,"_Tumor")
   
   
   myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/matched_TCR/single_clonotype_match/STB_blood_expanded/UMAP_matched_clonotype",i,".pdf")
   pdf(myoutf,width=6.5,height=5)
   print(DimPlot(object = PW_Adult_kit.integrated,group.by='single_clonotype',shape.by="orig.ident",reduction="umap",
                 cols = col,pt.size = 0.2))
   dev.off()
   
}









#Plot single clonotypes that matches between STB and expanded in both skin and tumor, frequency>2
#Point: IFNG Trm has functional counterpart in blood as well





barcode_STB_blood<-row.names(PW_Adult_kit.integrated@meta.data[PW_Adult_kit.integrated@meta.data$STB_match=="STB_blood",])
tag<-unique(barcode_info_all_brief)$barcode%in%barcode_STB_blood
clonotype<-barcode_info_all_brief[tag,]
clonotype[,"frequency"]<-rep(0,nrow(clonotype))
for(i in 1:length(unique(clonotype$labeled_clonotype_id)))
{
  tag<-clonotype$labeled_clonotype_id==unique(clonotype$labeled_clonotype_id)[i]
  xx<-clonotype[tag,]
  clonotype$frequency[tag]=nrow(xx)
  
}

tag<-clonotype$frequency>0
blood_expanded<-clonotype[tag,]



barcode_STB_tumor<-row.names(PW_Adult_kit.integrated@meta.data[PW_Adult_kit.integrated@meta.data$STB_match=="STB_tumor",])
tag<-unique(barcode_info_all_brief)$barcode%in%barcode_STB_tumor
clonotype<-barcode_info_all_brief[tag,]
clonotype[,"frequency"]<-rep(0,nrow(clonotype))
for(i in 1:length(unique(clonotype$labeled_clonotype_id)))
{
   tag<-clonotype$labeled_clonotype_id==unique(clonotype$labeled_clonotype_id)[i]
   xx<-clonotype[tag,]
   clonotype$frequency[tag]=nrow(xx)
   
}

tag<-clonotype$frequency>2
tumor_expanded<-clonotype[tag,]



barcode_STB_skin<-row.names(PW_Adult_kit.integrated@meta.data[PW_Adult_kit.integrated@meta.data$STB_match=="STB_skin",])
tag<-unique(barcode_info_all_brief)$barcode%in%barcode_STB_skin
clonotype<-barcode_info_all_brief[tag,]
clonotype[,"frequency"]<-rep(0,nrow(clonotype))

for(i in 1:length(unique(clonotype$labeled_clonotype_id)))
{
  tag<-clonotype$labeled_clonotype_id==unique(clonotype$labeled_clonotype_id)[i]
  xx<-clonotype[tag,]
  clonotype$frequency[tag]=nrow(xx)
  
}

tag<-clonotype$frequency>2
skin_expanded<-clonotype[tag,]

tag<-all_info_STB$labeled_clonotype_id%in%blood_expanded$labeled_clonotype_id&all_info_STB$labeled_clonotype_id.x%in%skin_expanded$labeled_clonotype_id&all_info_STB$labeled_clonotype_id.y%in%tumor_expanded$labeled_clonotype_id

all_info_STB_ST_expanded<-all_info_STB[tag,]



enrich_skin_all_clonotypes<-matrix(0,nrow(all_info_STB_ST_expanded),10)
row.names(enrich_skin_all_clonotypes)=paste0("clonotype_",seq(1,nrow(all_info_STB_ST_expanded),1))
colnames(enrich_skin_all_clonotypes)<-seq(0,9,1)

enrich_tumor_all_clonotypes<-matrix(0,nrow(all_info_STB_ST_expanded),10)
row.names(enrich_tumor_all_clonotypes)=paste0("clonotype_",seq(1,nrow(all_info_STB_ST_expanded),1))
colnames(enrich_tumor_all_clonotypes)<-seq(0,9,1)












col<-c('#2ca25f',"#de2d26","#9ecae1",'#f0f0f0')

for(i in 1:nrow(all_info_STB_ST_expanded))
{
  cat("\r",i)
  clonotype_STB<-combine(all_info_STB_ST_expanded$labeled_clonotype_id.x[i],
                         all_info_STB_ST_expanded$labeled_clonotype_id.y[i],
                         all_info_STB_ST_expanded$labeled_clonotype_id[i])
  
  
  tag<-barcode_info_all_brief$labeled_clonotype_id%in%clonotype_STB
  barcode_STB<-barcode_info_all_brief$barcode[tag]
  barcode_STB_skin<-barcode_STB[grep("skin",barcode_STB)]
  barcode_STB_blood<-barcode_STB[grep("blood",barcode_STB)]
  barcode_STB_tumor<-barcode_STB[grep("tumor",barcode_STB)]
  
  
  PW_Adult_kit.integrated@meta.data[,"single_clonotype"]<-rep("Others",nrow( PW_Adult_kit.integrated@meta.data))
  tag1<-row.names(PW_Adult_kit.integrated@meta.data)%in%barcode_STB_skin
  tag2<-row.names(PW_Adult_kit.integrated@meta.data)%in%barcode_STB_blood
  tag3<-row.names(PW_Adult_kit.integrated@meta.data)%in%barcode_STB_tumor
  
  PW_Adult_kit.integrated@meta.data$single_clonotype[tag1]<-paste0("Clonotype",i,"_Skin")
  PW_Adult_kit.integrated@meta.data$single_clonotype[tag2]<-paste0("Clonotype",i,"_Blood")
  PW_Adult_kit.integrated@meta.data$single_clonotype[tag3]<-paste0("Clonotype",i,"_Tumor")
  
  xx<-PW_Adult_kit.integrated@meta.data[tag1,]
  yy<-PW_Adult_kit.integrated@meta.data[tag3,]
  
  chi_STB_skin<-matrix(0,2,10)
  row.names(chi_STB_skin)=c("match","not match")
  colnames(chi_STB_skin)<-seq(0,9,1)
  
  
  chi_STB_tumor<-matrix(0,2,10)
  row.names(chi_STB_tumor)=c("match","not match")
  colnames(chi_STB_tumor)<-seq(0,9,1)
  
  
  
  for (j in 1:length(unique(xx$integrated_snn_res.0.3)))
  {
    chi_STB_skin["match",unique(xx$integrated_snn_res.0.3)[j]]=sum(xx$integrated_snn_res.0.3==unique(xx$integrated_snn_res.0.3)[j])
    chi_STB_skin["not match",unique(xx$integrated_snn_res.0.3)[j]]=sum((PW_Adult_kit.integrated@meta.data[PW_Adult_kit.integrated@meta.data$label=="skin",])$integrated_snn_res.0.3==unique(xx$integrated_snn_res.0.3)[j])-sum(xx$integrated_snn_res.0.3==unique(xx$integrated_snn_res.0.3)[j])
  }
  
  
  
  for (j in 1:length(unique(yy$integrated_snn_res.0.3)))
  {
    chi_STB_tumor["match",unique(yy$integrated_snn_res.0.3)[j]]=sum(yy$integrated_snn_res.0.3==unique(yy$integrated_snn_res.0.3)[j])
    chi_STB_tumor["not match",unique(yy$integrated_snn_res.0.3)[j]]=sum((PW_Adult_kit.integrated@meta.data[PW_Adult_kit.integrated@meta.data$label=="tumor",])$integrated_snn_res.0.3==unique(yy$integrated_snn_res.0.3)[j])-sum(yy$integrated_snn_res.0.3==unique(yy$integrated_snn_res.0.3)[j])
  }
  
  
  chi_skin<-chisq.test(chi_STB_skin)
  chi_observed<-chi_skin$observed
  chi_expected<-chi_skin$expected
  ROE_skin<-chi_observed/chi_expected
  
  
  
  chi_tumor<-chisq.test(chi_STB_tumor)
  chi_observed<-chi_tumor$observed
  chi_expected<-chi_tumor$expected
  ROE_tumor<-chi_observed/chi_expected
  
  enrich_skin_all_clonotypes[i,]<-ROE_skin[1,]
  enrich_tumor_all_clonotypes[i,]<-ROE_tumor[1,]
  
  
  
  #myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/matched_TCR/single_clonotype_match/STB_tumor_skin_expanded/UMAP_matched_clonotype",i,".pdf")
  #pdf(myoutf,width=6.5,height=5)
  #print(DimPlot(object = PW_Adult_kit.integrated,group.by='single_clonotype',shape.by="orig.ident",reduction="umap",
              #  cols = col,pt.size = 0.2))
  #dev.off()
  
}



myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/matched_TCR/single_clonotype_match/STB_tumor_skin_expanded/skin_enrichment_score_each_cluster.xls"
write.table(enrich_skin_all_clonotypes,myoutf,sep='\t')

myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/matched_TCR/single_clonotype_match/STB_tumor_skin_expanded/tumor_enrichment_score_each_cluster.xls"
write.table(enrich_tumor_all_clonotypes,myoutf,sep='\t')






enrich_skin_all_clonotypes[is.na(enrich_skin_all_clonotypes)] <- 0

enrich_tumor_all_clonotypes[is.na( enrich_tumor_all_clonotypes)] <- 0

enrich_all<-matrix(0,300,4)
colnames(enrich_all)<-c("clonotype","tissue","cluster","index")
enrich_all[,"clonotype"]<-rep(paste0("clonotype_",seq(1,nrow(all_info_STB_ST_expanded),1)),20)
enrich_all[,"tissue"]<-combine(rep("skin",150),rep("tumor",150))
enrich_all[,"cluster"]<-rep(seq(0,9,1),each=15)
enrich_all[,"index"]<-combine(combine(enrich_skin_all_clonotypes),combine(enrich_tumor_all_clonotypes))

enrich_all<-as.data.frame(enrich_all)

enrich_all$index=as.numeric(levels(enrich_all$index))[enrich_all$index]

enrich_all_skin<-enrich_all[1:150,]
enrich_all_tumor<-enrich_all[151:300,]





colourCount = 15
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

enrich_skin_plot<-enrich_all_skin[as.character(enrich_all_skin$cluster)%in%c("0","1","3","6","7"),]
enrich_tumor_plot<-enrich_all_tumor[as.character(enrich_all_tumor$cluster)%in%as.character(c(0,1,3,6,7)),]



myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/matched_TCR/single_clonotype_match/STB_tumor_skin_expanded/line_plot_enrichment_cluster_in_skin_tumor_skin.pdf"
pdf(myoutf,width=30,height=20)
print(ggplot(data=enrich_skin_plot, aes(x=cluster, y=index, group=clonotype)) + 
         scale_x_discrete(limits=unique(enrich_skin_plot$cluster))+
        scale_y_continuous(expand = c(0,0.025))+
         geom_line(aes(color=clonotype),size=1.5)+
         geom_hline(yintercept=1, linetype="dashed", 
                    color = "red", size=2)+
         scale_color_manual(values=getPalette(colourCount))+
         geom_point(aes(color=clonotype,shape=tissue),size=8)
      +theme(axis.text.x = element_text(angle = 60, hjust = 1,size=6),axis.title.x = element_blank(),
             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black")),
      legend.text=element_text(size=3),
      legend.key.size = unit(0.1, "cm"),
      legend.key.width = unit(0.1,"cm"))

dev.off()








myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/matched_TCR/single_clonotype_match/STB_tumor_skin_expanded/line_plot_enrichment_cluster_in_skin_tumor_tumor.pdf"
pdf(myoutf,width=30,height=20)
print(ggplot(data=enrich_tumor_plot, aes(x=cluster, y=index, group=clonotype)) + 
         scale_x_discrete(limits=unique(enrich_tumor_plot$cluster))+
        scale_y_continuous(expand = c(0,0.025))+
         geom_line(aes(color=clonotype),size=1.5)+
         geom_hline(yintercept=1, linetype="dashed", 
                    color = "red", size=2)+
         scale_color_manual(values=getPalette(colourCount))+
         geom_point(aes(color=clonotype,shape=tissue),size=8)
      +theme(axis.text.x = element_text(angle = 60, hjust = 1,size=6),axis.title.x = element_blank(),
             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black")),
      legend.text=element_text(size=3),
      legend.key.size = unit(0.1, "cm"),
      legend.key.width = unit(0.1,"cm"))

dev.off()










myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/matched_TCR/single_clonotype_match/STB_tumor_skin_expanded/line_plot_enrichment_each_cluster.pdf"
pdf(myoutf,width=15,height=30)
print(ggplot(data=enrich_all, aes(x=cluster, y=index, group=clonotype)) + 
        scale_x_discrete(limits=unique(enrich_all$cluster))+
        geom_line(aes(color=clonotype),size=1.5)+
         geom_hline(yintercept=1, linetype="dashed", 
                    color = "red", size=2)+
        scale_color_manual(values=getPalette(colourCount))+
        geom_point(aes(color=clonotype,shape=tissue),size=8)
      +theme(axis.text.x = element_text(angle = 60, hjust = 1,size=6),axis.title.x = element_blank(),
             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black")),
      legend.text=element_text(size=3),
      legend.key.size = unit(0.1, "cm"),
      legend.key.width = unit(0.1,"cm"))

dev.off()






myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/matched_TCR/single_clonotype_match/STB_tumor_skin_expanded/line_plot_enrichment_each_cluster_skin.pdf"
pdf(myoutf,width=15,height=20)
print(ggplot(data=enrich_all_skin, aes(x=cluster, y=index, group=clonotype)) + 
        scale_x_discrete(limits=unique(enrich_all$cluster))+
        geom_line(aes(color=clonotype),size=1.5)+
        scale_color_manual(values=getPalette(colourCount))+
        geom_hline(yintercept=1, linetype="dashed", 
                    color = "red", size=2)+
        geom_point(aes(color=clonotype,shape=tissue),size=8)
      +theme(axis.text.x = element_text(angle = 60, hjust = 1,size=6),axis.title.x = element_blank(),
             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black")),
      legend.text=element_text(size=3),
      legend.key.size = unit(0.1, "cm"),
      legend.key.width = unit(0.1,"cm"))

dev.off()


t.test(enrich_all_skin$index[1:15],enrich_all_skin$index[16:30],alternative = "less")


myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/matched_TCR/single_clonotype_match/STB_tumor_skin_expanded/line_plot_enrichment_each_cluster_tumor.pdf"
pdf(myoutf,width=15,height=20)
print(ggplot(data=enrich_all_tumor, aes(x=cluster, y=index, group=clonotype)) + 
        scale_x_discrete(limits=unique(enrich_all$cluster))+
        geom_line(aes(color=clonotype),size=1.5)+
        scale_color_manual(values=getPalette(colourCount))+
        geom_point(aes(color=clonotype,shape=tissue),size=8)
        +geom_hline(yintercept=1, linetype="dashed", 
                    color = "red", size=2)
      +theme(axis.text.x = element_text(angle = 60, hjust = 1,size=6),axis.title.x = element_blank(),
             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black")),
      legend.text=element_text(size=3),
      legend.key.size = unit(0.1, "cm"),
      legend.key.width = unit(0.1,"cm"))
dev.off()


t.test(enrich_all_tumor$index[1:15],enrich_all_tumor$index[16:30],alternative = "less")
















enrich_skin_all_clonotypes[is.na(enrich_skin_all_clonotypes)] <- 0

enrich_tumor_all_clonotypes[is.na( enrich_tumor_all_clonotypes)] <- 0


enrich_skin<-enrich_skin_all_clonotypes[,c(1,2,4,7,8)]
enrich_tumor<-enrich_tumor_all_clonotypes[,c(1,2,4,7,8)]

xx<-max.col(enrich_skin)
yy<-max.col(enrich_tumor)

same<-sum(xx==yy)/length(xx)
not_same<-1-same

library(plotrix)
slices <- c(same,not_same)
lbls <- c("Most enriched in the same cluster", "Most enriched in different clusters")

myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/matched_TCR/single_clonotype_match/STB_tumor_skin_expanded/pie_chart_if_enriched_same_cluster.pdf"
pdf(myoutf,width=10,height=10)
print(pie3D(slices,labels=lbls,explode=0.1,
            main=("Pie Chart of Enrichments"),
            col=c("#de2d26","#31a354"),
      labelpos=c(1.5,4.8)))
dev.off()



















#Plot single clonotypes that matches between ST and expanded in both skin and tumor
#Point: IFNG Trm has functional counterpart in blood as well


barcode_ST_skin<-row.names(PW_Adult_kit.integrated@meta.data[PW_Adult_kit.integrated@meta.data$ST_match=="ST_skin",])
tag<-unique(barcode_info_all_brief)$barcode%in%barcode_ST_skin
clonotype<-barcode_info_all_brief[tag,]
clonotype[,"frequency"]<-rep(0,nrow(clonotype))
for(i in 1:length(unique(clonotype$labeled_clonotype_id)))
{
  tag<-clonotype$labeled_clonotype_id==unique(clonotype$labeled_clonotype_id)[i]
  xx<-clonotype[tag,]
  clonotype$frequency[tag]=nrow(xx)
  
}

tag<-clonotype$frequency>2
skin_expanded<-clonotype[tag,]


barcode_ST_tumor<-row.names(PW_Adult_kit.integrated@meta.data[PW_Adult_kit.integrated@meta.data$ST_match=="ST_tumor",])
tag<-unique(barcode_info_all_brief)$barcode%in%barcode_ST_tumor
clonotype<-barcode_info_all_brief[tag,]
clonotype[,"frequency"]<-rep(0,nrow(clonotype))

for(i in 1:length(unique(clonotype$labeled_clonotype_id)))
{
  tag<-clonotype$labeled_clonotype_id==unique(clonotype$labeled_clonotype_id)[i]
  xx<-clonotype[tag,]
  clonotype$frequency[tag]=nrow(xx)
  
}

tag<-clonotype$frequency>2
tumor_expanded<-clonotype[tag,]

tag<-all_info_ST$labeled_clonotype_id.y%in%tumor_expanded$labeled_clonotype_id & all_info_ST$labeled_clonotype_id.x%in%skin_expanded$labeled_clonotype_id

all_info_ST_ST_expanded<- all_info_ST[tag,]



enrich_skin_all_clonotypes<-matrix(0,nrow(all_info_ST_ST_expanded),10)
row.names(enrich_skin_all_clonotypes)=paste0("clonotype_",seq(1,nrow(all_info_ST_ST_expanded),1))
colnames(enrich_skin_all_clonotypes)<-seq(0,9,1)

enrich_tumor_all_clonotypes<-matrix(0,nrow(all_info_ST_ST_expanded),10)
row.names(enrich_tumor_all_clonotypes)=paste0("clonotype_",seq(1,nrow(all_info_ST_ST_expanded),1))
colnames(enrich_tumor_all_clonotypes)<-seq(0,9,1)







col<-c("#de2d26","#9ecae1",'#f0f0f0')

for(i in 1:nrow(all_info_ST_ST_expanded))
{
  cat("\r",i)
  clonotype_ST<-combine(all_info_ST_ST_expanded$labeled_clonotype_id.x[i],
                         all_info_ST_ST_expanded$labeled_clonotype_id.y[i])
  
  
  tag<-barcode_info_all_brief$labeled_clonotype_id%in%clonotype_ST
  barcode_ST<-barcode_info_all_brief$barcode[tag]
  barcode_ST_skin<-barcode_ST[grep("skin",barcode_ST)]
 
  barcode_ST_tumor<-barcode_ST[grep("tumor",barcode_ST)]
  
  
  PW_Adult_kit.integrated@meta.data[,"single_clonotype"]<-rep("Others",nrow( PW_Adult_kit.integrated@meta.data))
  tag1<-row.names(PW_Adult_kit.integrated@meta.data)%in%barcode_ST_skin
  tag3<-row.names(PW_Adult_kit.integrated@meta.data)%in%barcode_ST_tumor
  
  PW_Adult_kit.integrated@meta.data$single_clonotype[tag1]<-paste0("Clonotype",i,"_Skin")
  PW_Adult_kit.integrated@meta.data$single_clonotype[tag3]<-paste0("Clonotype",i,"_Tumor")
  
  xx<-PW_Adult_kit.integrated@meta.data[tag1,]
  yy<-PW_Adult_kit.integrated@meta.data[tag3,]
  
  chi_ST_skin<-matrix(0,2,10)
  row.names(chi_ST_skin)=c("match","not match")
  colnames(chi_ST_skin)<-seq(0,9,1)
  
  
  chi_ST_tumor<-matrix(0,2,10)
  row.names(chi_ST_tumor)=c("match","not match")
  colnames(chi_ST_tumor)<-seq(0,9,1)
  
  for (j in 1:length(unique(xx$integrated_snn_res.0.3)))
  {
    chi_ST_skin["match",unique(xx$integrated_snn_res.0.3)[j]]=sum(xx$integrated_snn_res.0.3==unique(xx$integrated_snn_res.0.3)[j])
    chi_ST_skin["not match",unique(xx$integrated_snn_res.0.3)[j]]=sum((PW_Adult_kit.integrated@meta.data[PW_Adult_kit.integrated@meta.data$label=="skin",])$integrated_snn_res.0.3==unique(xx$integrated_snn_res.0.3)[j])-sum(xx$integrated_snn_res.0.3==unique(xx$integrated_snn_res.0.3)[j])
  }
  
  
  
  for (j in 1:length(unique(yy$integrated_snn_res.0.3)))
  {
    chi_ST_tumor["match",unique(yy$integrated_snn_res.0.3)[j]]=sum(yy$integrated_snn_res.0.3==unique(yy$integrated_snn_res.0.3)[j])
    chi_ST_tumor["not match",unique(yy$integrated_snn_res.0.3)[j]]=sum((PW_Adult_kit.integrated@meta.data[PW_Adult_kit.integrated@meta.data$label=="tumor",])$integrated_snn_res.0.3==unique(yy$integrated_snn_res.0.3)[j])-sum(yy$integrated_snn_res.0.3==unique(yy$integrated_snn_res.0.3)[j])
  }
  
  
  chi_skin<-chisq.test(chi_ST_skin)
  chi_observed<-chi_skin$observed
  chi_expected<-chi_skin$expected
  ROE_skin<-chi_observed/chi_expected
  
  
  
  chi_tumor<-chisq.test(chi_ST_tumor)
  chi_observed<-chi_tumor$observed
  chi_expected<-chi_tumor$expected
  ROE_tumor<-chi_observed/chi_expected
  
  enrich_skin_all_clonotypes[i,]<-ROE_skin[1,]
  enrich_tumor_all_clonotypes[i,]<-ROE_tumor[1,]
  
  
  
  
  
  
  #myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/matched_TCR/single_clonotype_match/ST_skin_tumor_expanded/UMAP_matched_clonotype",i,".pdf")
  #pdf(myoutf,width=6.5,height=5)
  #print(DimPlot(object = PW_Adult_kit.integrated,group.by='single_clonotype',shape.by="orig.ident",reduction="umap",
  #              cols = col,pt.size = 0.2))
  #dev.off()
  
}



myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/matched_TCR/single_clonotype_match/ST_skin_tumor_expanded/skin_enrichment_score_each_cluster.xls"
write.table(enrich_skin_all_clonotypes,myoutf,sep='\t')

myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/matched_TCR/single_clonotype_match/ST_skin_tumor_expanded/tumor_enrichment_score_each_cluster.xls"
write.table(enrich_tumor_all_clonotypes,myoutf,sep='\t')





enrich_skin_all_clonotypes[is.na(enrich_skin_all_clonotypes)] <- 0

enrich_tumor_all_clonotypes[is.na( enrich_tumor_all_clonotypes)] <- 0

enrich_all<-matrix(0,360,4)
colnames(enrich_all)<-c("clonotype","tissue","cluster","index")
enrich_all[,"clonotype"]<-rep(paste0("clonotype_",seq(1,nrow(all_info_ST_ST_expanded),1)),20)
enrich_all[,"tissue"]<-combine(rep("skin",180),rep("tumor",180))
enrich_all[,"cluster"]<-rep(seq(0,9,1),each=18)
enrich_all[,"index"]<-combine(combine(enrich_skin_all_clonotypes),combine(enrich_tumor_all_clonotypes))

enrich_all<-as.data.frame(enrich_all)

enrich_all$index=as.numeric(levels(enrich_all$index))[enrich_all$index]

enrich_all_skin<-enrich_all[1:180,]
enrich_all_tumor<-enrich_all[181:360,]





colourCount = 18
getPalette = colorRampPalette(brewer.pal(9, "Set1"))


myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/matched_TCR/single_clonotype_match/ST_skin_tumor_expanded/line_plot_enrichment_each_cluster.pdf"
pdf(myoutf,width=15,height=30)
print(ggplot(data=enrich_all, aes(x=cluster, y=index, group=clonotype)) + 
        scale_x_discrete(limits=unique(enrich_all$cluster))+
        geom_line(aes(color=clonotype),size=1.5)+
        scale_color_manual(values=getPalette(colourCount))+
        geom_point(aes(color=clonotype,shape=tissue),size=8)
      +theme(axis.text.x = element_text(angle = 60, hjust = 1,size=6),axis.title.x = element_blank(),
             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black")),
      legend.text=element_text(size=3),
      legend.key.size = unit(0.1, "cm"),
      legend.key.width = unit(0.1,"cm"))

dev.off()






myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/matched_TCR/single_clonotype_match/ST_skin_tumor_expanded/line_plot_enrichment_each_cluster_skin.pdf"
pdf(myoutf,width=15,height=20)
print(ggplot(data=enrich_all_skin, aes(x=cluster, y=index, group=clonotype)) + 
        scale_x_discrete(limits=unique(enrich_all$cluster))+
        geom_line(aes(color=clonotype),size=1.5)+
        scale_color_manual(values=getPalette(colourCount))+
        geom_hline(yintercept=1, linetype="dashed", 
                   color = "red", size=2)+
        geom_point(aes(color=clonotype,shape=tissue),size=8)
      +theme(axis.text.x = element_text(angle = 60, hjust = 1,size=6),axis.title.x = element_blank(),
             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black")),
      legend.text=element_text(size=3),
      legend.key.size = unit(0.1, "cm"),
      legend.key.width = unit(0.1,"cm"))

dev.off()


t.test(enrich_all_skin$index[1:15],enrich_all_skin$index[16:30],alternative = "less")


myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/matched_TCR/single_clonotype_match/ST_skin_tumor_expanded/line_plot_enrichment_each_cluster_tumor.pdf"
pdf(myoutf,width=15,height=20)
print(ggplot(data=enrich_all_tumor, aes(x=cluster, y=index, group=clonotype)) + 
        scale_x_discrete(limits=unique(enrich_all$cluster))+
        geom_line(aes(color=clonotype),size=1.5)+
        scale_color_manual(values=getPalette(colourCount))+
        geom_point(aes(color=clonotype,shape=tissue),size=8)
      +geom_hline(yintercept=1, linetype="dashed", 
                  color = "red", size=2)
      +theme(axis.text.x = element_text(angle = 60, hjust = 1,size=6),axis.title.x = element_blank(),
             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black")),
      legend.text=element_text(size=3),
      legend.key.size = unit(0.1, "cm"),
      legend.key.width = unit(0.1,"cm"))
dev.off()


t.test(enrich_all_tumor$index[1:15],enrich_all_tumor$index[16:30],alternative = "less")













enrich_skin_plot<-enrich_all_skin[as.character(enrich_all_skin$cluster)%in%c("0","1","3","6","7"),]
enrich_tumor_plot<-enrich_all_tumor[as.character(enrich_all_tumor$cluster)%in%as.character(c(0,1,3,6,7)),]



myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/matched_TCR/single_clonotype_match/ST_skin_tumor_expanded/line_plot_enrichment_cluster_in_skin_tumor_skin.pdf"
pdf(myoutf,width=30,height=20)
print(ggplot(data=enrich_skin_plot, aes(x=cluster, y=index, group=clonotype)) + 
         scale_x_discrete(limits=unique(enrich_skin_plot$cluster))+
        scale_y_continuous(expand = c(0,0.025))+
         geom_line(aes(color=clonotype),size=1.5)+
         geom_hline(yintercept=1, linetype="dashed", 
                    color = "red", size=2)+
         scale_color_manual(values=getPalette(colourCount))+
         geom_point(aes(color=clonotype,shape=tissue),size=8)
      +theme(axis.text.x = element_text(angle = 60, hjust = 1,size=6),axis.title.x = element_blank(),
             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black")),
      legend.text=element_text(size=3),
      legend.key.size = unit(0.1, "cm"),
      legend.key.width = unit(0.1,"cm"))

dev.off()








myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/matched_TCR/single_clonotype_match/ST_skin_tumor_expanded/line_plot_enrichment_cluster_in_skin_tumor_tumor.pdf"
pdf(myoutf,width=30,height=20)
print(ggplot(data=enrich_tumor_plot, aes(x=cluster, y=index, group=clonotype)) + 
         scale_x_discrete(limits=unique(enrich_tumor_plot$cluster))+
        scale_y_continuous(expand = c(0,0.025))+
         geom_line(aes(color=clonotype),size=1.5)+
         geom_hline(yintercept=1, linetype="dashed", 
                    color = "red", size=2)+
         scale_color_manual(values=getPalette(colourCount))+
         geom_point(aes(color=clonotype,shape=tissue),size=8)
      +theme(axis.text.x = element_text(angle = 60, hjust = 1,size=6),axis.title.x = element_blank(),
             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black")),
      legend.text=element_text(size=3),
      legend.key.size = unit(0.1, "cm"),
      legend.key.width = unit(0.1,"cm"))

dev.off()













enrich_skin_all_clonotypes[is.na(enrich_skin_all_clonotypes)] <- 0

enrich_tumor_all_clonotypes[is.na( enrich_tumor_all_clonotypes)] <- 0


enrich_skin<-enrich_skin_all_clonotypes[,c(1,2,4,7,8)]
enrich_tumor<-enrich_tumor_all_clonotypes[,c(1,2,4,7,8)]

xx<-max.col(enrich_skin)
yy<-max.col(enrich_tumor)

same<-sum(xx==yy)/length(xx)
not_same<-1-same

library(plotrix)
slices <- c(same,not_same)
lbls <- c("Most enriched in the same cluster", "Most enriched in different clusters")

myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/matched_TCR/single_clonotype_match/ST_skin_tumor_expanded/pie_chart_if_enriched_same_cluster.pdf"
pdf(myoutf,width=10,height=10)
print(pie3D(slices,labels=lbls,explode=0.1,
            main=("Pie Chart of Enrichments"),
            col=c("#de2d26","#31a354"),
            labelpos=c(1.5,4.8)))
dev.off()




























#Plot single clonotypes that matches between BT and expanded in both blood and tumor
#Point: IFNG Trm has functional counterpart in blood as well


barcode_BT_blood<-row.names(PW_Adult_kit.integrated@meta.data[PW_Adult_kit.integrated@meta.data$BT_match=="BT_blood",])
tag<-unique(barcode_info_all_brief)$barcode%in%barcode_BT_blood
clonotype<-barcode_info_all_brief[tag,]
clonotype[,"frequency"]<-rep(0,nrow(clonotype))
for(i in 1:length(unique(clonotype$labeled_clonotype_id)))
{
  tag<-clonotype$labeled_clonotype_id==unique(clonotype$labeled_clonotype_id)[i]
  xx<-clonotype[tag,]
  clonotype$frequency[tag]=nrow(xx)
  
}

tag<-clonotype$frequency>0
blood_expanded<-clonotype[tag,]


barcode_BT_tumor<-row.names(PW_Adult_kit.integrated@meta.data[PW_Adult_kit.integrated@meta.data$BT_match=="BT_tumor",])
tag<-unique(barcode_info_all_brief)$barcode%in%barcode_BT_tumor
clonotype<-barcode_info_all_brief[tag,]
clonotype[,"frequency"]<-rep(0,nrow(clonotype))

for(i in 1:length(unique(clonotype$labeled_clonotype_id)))
{
  tag<-clonotype$labeled_clonotype_id==unique(clonotype$labeled_clonotype_id)[i]
  xx<-clonotype[tag,]
  clonotype$frequency[tag]=nrow(xx)
  
}

tag<-clonotype$frequency>2
tumor_expanded<-clonotype[tag,]

tag<-all_info_BT$labeled_clonotype_id.y%in%tumor_expanded$labeled_clonotype_id & all_info_BT$labeled_clonotype_id.x%in%blood_expanded$labeled_clonotype_id

all_info_BT_BT_expanded<- all_info_BT[tag,]



enrich_blood_all_clonotypes<-matrix(0,nrow(all_info_BT_BT_expanded),10)
row.names(enrich_blood_all_clonotypes)=paste0("clonotype_",seq(1,nrow(all_info_BT_BT_expanded),1))
colnames(enrich_blood_all_clonotypes)<-seq(0,9,1)

enrich_tumor_all_clonotypes<-matrix(0,nrow(all_info_BT_BT_expanded),10)
row.names(enrich_tumor_all_clonotypes)=paste0("clonotype_",seq(1,nrow(all_info_BT_BT_expanded),1))
colnames(enrich_tumor_all_clonotypes)<-seq(0,9,1)







col<-c('#2ca25f',"#3333cc",'#f0f0f0')

for(i in 1:nrow(all_info_BT_BT_expanded))
{
  cat("\r",i)
  clonotype_BT<-combine(all_info_BT_BT_expanded$labeled_clonotype_id.x[i],
                        all_info_BT_BT_expanded$labeled_clonotype_id.y[i])
  
  
  tag<-barcode_info_all_brief$labeled_clonotype_id%in%clonotype_BT
  barcode_BT<-barcode_info_all_brief$barcode[tag]
  barcode_BT_blood<-barcode_BT[grep("blood",barcode_BT)]
  
  barcode_BT_tumor<-barcode_BT[grep("tumor",barcode_BT)]
  
  
  PW_Adult_kit.integrated@meta.data[,"single_clonotype"]<-rep("Others",nrow( PW_Adult_kit.integrated@meta.data))
  tag1<-row.names(PW_Adult_kit.integrated@meta.data)%in%barcode_BT_blood
  tag3<-row.names(PW_Adult_kit.integrated@meta.data)%in%barcode_BT_tumor
  
  PW_Adult_kit.integrated@meta.data$single_clonotype[tag1]<-paste0("Clonotype",i,"_blood")
  PW_Adult_kit.integrated@meta.data$single_clonotype[tag3]<-paste0("Clonotype",i,"_Tumor")
  
  xx<-PW_Adult_kit.integrated@meta.data[tag1,]
  yy<-PW_Adult_kit.integrated@meta.data[tag3,]
  
  chi_BT_blood<-matrix(0,2,10)
  row.names(chi_BT_blood)=c("match","not match")
  colnames(chi_BT_blood)<-seq(0,9,1)
  
  
  chi_BT_tumor<-matrix(0,2,10)
  row.names(chi_BT_tumor)=c("match","not match")
  colnames(chi_BT_tumor)<-seq(0,9,1)
  
  for (j in 1:length(unique(xx$integrated_snn_res.0.3)))
  {
    chi_BT_blood["match",unique(xx$integrated_snn_res.0.3)[j]]=sum(xx$integrated_snn_res.0.3==unique(xx$integrated_snn_res.0.3)[j])
    chi_BT_blood["not match",unique(xx$integrated_snn_res.0.3)[j]]=sum((PW_Adult_kit.integrated@meta.data[PW_Adult_kit.integrated@meta.data$label=="blood",])$integrated_snn_res.0.3==unique(xx$integrated_snn_res.0.3)[j])-sum(xx$integrated_snn_res.0.3==unique(xx$integrated_snn_res.0.3)[j])
  }
  
  
  
  for (j in 1:length(unique(yy$integrated_snn_res.0.3)))
  {
    chi_BT_tumor["match",unique(yy$integrated_snn_res.0.3)[j]]=sum(yy$integrated_snn_res.0.3==unique(yy$integrated_snn_res.0.3)[j])
    chi_BT_tumor["not match",unique(yy$integrated_snn_res.0.3)[j]]=sum((PW_Adult_kit.integrated@meta.data[PW_Adult_kit.integrated@meta.data$label=="tumor",])$integrated_snn_res.0.3==unique(yy$integrated_snn_res.0.3)[j])-sum(yy$integrated_snn_res.0.3==unique(yy$integrated_snn_res.0.3)[j])
  }
  
  
  chi_blood<-chisq.test(chi_BT_blood)
  chi_observed<-chi_blood$observed
  chi_expected<-chi_blood$expected
  ROE_blood<-chi_observed/chi_expected
  
  
  
  chi_tumor<-chisq.test(chi_BT_tumor)
  chi_observed<-chi_tumor$observed
  chi_expected<-chi_tumor$expected
  ROE_tumor<-chi_observed/chi_expected
  
  enrich_blood_all_clonotypes[i,]<-ROE_blood[1,]
  enrich_tumor_all_clonotypes[i,]<-ROE_tumor[1,]
  
  
  
  
  
  
  #myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/matched_TCR/single_clonotype_match/BT_blood_tumor_expanded/UMAP_matched_clonotype",i,".pdf")
  #pdf(myoutf,width=6.5,height=5)
  #print(DimPlot(object = PW_Adult_kit.integrated,group.by='single_clonotype',shape.by="orig.ident",reduction="umap",
  #              cols = col,pt.size = 0.2))
  #dev.off()
  
}




myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/matched_TCR/single_clonotype_match/BT_blood_tumor_expanded/blood_enrichment_score_each_cluster.xls"
write.table( enrich_blood_all_clonotypes,myoutf,sep='\t')

myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/matched_TCR/single_clonotype_match/BT_blood_tumor_expanded/tumor_enrichment_score_each_cluster.xls"
write.table(enrich_tumor_all_clonotypes,myoutf,sep='\t')








enrich_blood_all_clonotypes[is.na(enrich_blood_all_clonotypes)] <- 0

enrich_tumor_all_clonotypes[is.na( enrich_tumor_all_clonotypes)] <- 0

enrich_all<-matrix(0,220,4)
colnames(enrich_all)<-c("clonotype","tissue","cluster","index")
enrich_all[,"clonotype"]<-rep(paste0("clonotype_",seq(1,nrow(all_info_BT_BT_expanded),1)),20)
enrich_all[,"tissue"]<-combine(rep("blood",110),rep("tumor",110))
enrich_all[,"cluster"]<-rep(seq(0,9,1),each=11)
enrich_all[,"index"]<-combine(combine(enrich_blood_all_clonotypes),combine(enrich_tumor_all_clonotypes))

enrich_all<-as.data.frame(enrich_all)

enrich_all$index=as.numeric(levels(enrich_all$index))[enrich_all$index]

enrich_all_blood<-enrich_all[1:110,]
enrich_all_tumor<-enrich_all[111:220,]





colourCount = 11
getPalette = colorRampPalette(brewer.pal(9, "Set1"))


myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/matched_TCR/single_clonotype_match/BT_blood_tumor_expanded/line_plot_enrichment_each_cluster.pdf"
pdf(myoutf,width=15,height=30)
print(ggplot(data=enrich_all, aes(x=cluster, y=index, group=clonotype)) + 
        scale_x_discrete(limits=unique(enrich_all$cluster))+
        geom_line(aes(color=clonotype),size=1.5)+
        scale_color_manual(values=getPalette(colourCount))+
        geom_point(aes(color=clonotype,shape=tissue),size=8)
      +theme(axis.text.x = element_text(angle = 60, hjust = 1,size=6),axis.title.x = element_blank(),
             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black")),
      legend.text=element_text(size=3),
      legend.key.size = unit(0.1, "cm"),
      legend.key.width = unit(0.1,"cm"))

dev.off()






myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/matched_TCR/single_clonotype_match/BT_blood_tumor_expanded/line_plot_enrichment_each_cluster_blood.pdf"
pdf(myoutf,width=15,height=20)
print(ggplot(data=enrich_all_blood, aes(x=cluster, y=index, group=clonotype)) + 
        scale_x_discrete(limits=unique(enrich_all$cluster))+
        geom_line(aes(color=clonotype),size=1.5)+
        scale_color_manual(values=getPalette(colourCount))+
        geom_hline(yintercept=1, linetype="dashed", 
                   color = "red", size=2)+
        geom_point(aes(color=clonotype,shape=tissue),size=8)
      +theme(axis.text.x = element_text(angle = 60, hjust = 1,size=6),axis.title.x = element_blank(),
             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black")),
      legend.text=element_text(size=3),
      legend.key.size = unit(0.1, "cm"),
      legend.key.width = unit(0.1,"cm"))

dev.off()


t.test(enrich_all_blood$index[1:10],enrich_all_blood$index[11:20],alternative = "less")


myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/matched_TCR/single_clonotype_match/BT_blood_tumor_expanded/line_plot_enrichment_each_cluster_tumor.pdf"
pdf(myoutf,width=15,height=20)
print(ggplot(data=enrich_all_tumor, aes(x=cluster, y=index, group=clonotype)) + 
        scale_x_discrete(limits=unique(enrich_all$cluster))+
        geom_line(aes(color=clonotype),size=1.5)+
        scale_color_manual(values=getPalette(colourCount))+
        geom_point(aes(color=clonotype,shape=tissue),size=8)
      +geom_hline(yintercept=1, linetype="dashed", 
                  color = "red", size=2)
      +theme(axis.text.x = element_text(angle = 60, hjust = 1,size=6),axis.title.x = element_blank(),
             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black")),
      legend.text=element_text(size=3),
      legend.key.size = unit(0.1, "cm"),
      legend.key.width = unit(0.1,"cm"))
dev.off()


t.test(enrich_all_tumor$index[1:15],enrich_all_tumor$index[16:30],alternative = "less")













enrich_blood_plot<-enrich_all_blood[as.character(enrich_all_blood$cluster)%in%c("0","1","3","6","7"),]
enrich_tumor_plot<-enrich_all_tumor[as.character(enrich_all_tumor$cluster)%in%as.character(c(0,1,3,6,7)),]



myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/matched_TCR/single_clonotype_match/BT_blood_tumor_expanded/line_plot_enrichment_cluster_in_blood_tumor_blood.pdf"
pdf(myoutf,width=30,height=20)
print(ggplot(data=enrich_blood_plot, aes(x=cluster, y=index, group=clonotype)) + 
        scale_x_discrete(limits=unique(enrich_blood_plot$cluster))+
        scale_y_continuous(expand = c(0,0.025))+
        geom_line(aes(color=clonotype),size=1.5)+
        geom_hline(yintercept=1, linetype="dashed", 
                   color = "red", size=2)+
        scale_color_manual(values=getPalette(colourCount))+
        geom_point(aes(color=clonotype,shape=tissue),size=8)
      +theme(axis.text.x = element_text(angle = 60, hjust = 1,size=6),axis.title.x = element_blank(),
             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black")),
      legend.text=element_text(size=3),
      legend.key.size = unit(0.1, "cm"),
      legend.key.width = unit(0.1,"cm"))

dev.off()








myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/matched_TCR/single_clonotype_match/BT_blood_tumor_expanded/line_plot_enrichment_cluster_in_blood_tumor_tumor.pdf"
pdf(myoutf,width=30,height=20)
print(ggplot(data=enrich_tumor_plot, aes(x=cluster, y=index, group=clonotype)) + 
        scale_x_discrete(limits=unique(enrich_tumor_plot$cluster))+
        scale_y_continuous(expand = c(0,0.025))+
        geom_line(aes(color=clonotype),size=1.5)+
        geom_hline(yintercept=1, linetype="dashed", 
                   color = "red", size=2)+
        scale_color_manual(values=getPalette(colourCount))+
        geom_point(aes(color=clonotype,shape=tissue),size=8)
      +theme(axis.text.x = element_text(angle = 60, hjust = 1,size=6),axis.title.x = element_blank(),
             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black")),
      legend.text=element_text(size=3),
      legend.key.size = unit(0.1, "cm"),
      legend.key.width = unit(0.1,"cm"))

dev.off()













enrich_blood_all_clonotypes[is.na(enrich_blood_all_clonotypes)] <- 0

enrich_tumor_all_clonotypes[is.na( enrich_tumor_all_clonotypes)] <- 0


enrich_blood<-enrich_blood_all_clonotypes[,c(1,2,4,7,8)]
enrich_tumor<-enrich_tumor_all_clonotypes[,c(1,2,4,7,8)]

xx<-max.col(enrich_blood)
yy<-max.col(enrich_tumor)

same<-sum(xx==yy)/length(xx)
not_same<-1-same

library(plotrix)
slices <- c(same,not_same)
lbls <- c("Most enriched in the same cluster", "Most enriched in different clusters")

myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/matched_TCR/single_clonotype_match/BT_blood_tumor_expanded/pie_chart_if_enriched_same_cluster.pdf"
pdf(myoutf,width=10,height=10)
print(pie3D(slices,labels=lbls,explode=0.1,
            main=("Pie Chart of Enrichments"),
            col=c("#de2d26","#31a354"),
            labelpos=c(1.5,4.8)))
dev.off()


























































































###########Find clonal matches level between different clusters##################

PW_Adult_kit.integrated@meta.data[,'clonotype']<-rep("NA",nrow(PW_Adult_kit.integrated@meta.data))

clonotype_info<-barcode_info_all[,c("barcode","labeled_clonotype_id")]

for(i in 1:nrow(PW_Adult_kit.integrated@meta.data))
{
  cat("\r",i)
  tag<-grep(row.names(PW_Adult_kit.integrated@meta.data)[i],clonotype_info$barcode)
  xx<-unique(clonotype_info[tag,])
  PW_Adult_kit.integrated@meta.data[i,"clonotype"]=xx$labeled_clonotype_id
}




tar<-c(1,3,5,6)

data_all<-data.frame()

for(i in 1:length(tar))
{
  cat('\r',i)
  myinf1<-paste0("/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scTCRseq/Run_",tar[i],"/TCR/data/skin/skin_clonotypes.csv")
  myinf2<-paste0("/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scTCRseq/Run_",tar[i],"/TCR/data/blood/blood_clonotypes.csv")
  myinf3<-paste0("/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scTCRseq/Run_",tar[i],"/TCR/data/tumor/tumor_clonotypes.csv")
  data_skin<-read.csv(myinf1)
  data_blood<-read.csv(myinf2)
  data_tumor<-read.csv(myinf3)
  data_skin[,"label"]<-rep("skin",nrow(data_skin))
  data_blood[,"label"]<-rep("blood",nrow(data_blood))
  data_tumor[,"label"]<-rep("tumor",nrow(data_tumor))
  
  data_skin[,"labeled_clonotype_id"]<-paste0(data_skin$clonotype_id,"_",data_skin$label,tar[i])
  data_blood[,"labeled_clonotype_id"]<-paste0(data_blood$clonotype_id,"_",data_blood$label,tar[i])
  data_tumor[,"labeled_clonotype_id"]<-paste0(data_tumor$clonotype_id,"_",data_tumor$label,tar[i])
  
  data<-rbind(data_skin,data_blood,data_tumor)
  
  data_all<-rbind(data_all,data)

}

PW_Adult_kit.integrated@meta.data[,'cdr3']<-rep("NA",nrow(PW_Adult_kit.integrated@meta.data))

for(i in 1:nrow(PW_Adult_kit.integrated@meta.data))
{
  cat('\r',i)
  tag<-PW_Adult_kit.integrated@meta.data$clonotype[i]==data_all$labeled_clonotype_id
  
  if(sum(tag)!=0){PW_Adult_kit.integrated@meta.data$cdr3[i]<-as.character(data_all$cdr3s_nt[tag])}
  else{PW_Adult_kit.integrated@meta.data$cdr3[i]="NA"}
  
  
}


library(gtools)

tar<-unique(PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.3)


all_comb<-permutations(length(as.character(tar)), 2, v=as.character(tar),
             set=T, repeats.allowed=F)


match_different_cluster<-matrix(0,10,10)
row.names(match_different_cluster)<-as.character(seq(0,9,1))
colnames(match_different_cluster)<-as.character(seq(0,9,1))


  for (j in 1:nrow(all_comb))
  {
    
    cat("\r",j)
    pair<-all_comb[j,]
    tag1<-pair[1]
    tag2<-pair[2]
    
    xx1<-PW_Adult_kit.integrated@meta.data[PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.3==tag1,]
    xx2<-PW_Adult_kit.integrated@meta.data[PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.3==tag2,]
    
    total_xx1<-length(unique(xx1$cdr3))
    total_xx2<-length(unique(xx2$cdr3))
    match_xx1_xx2<-sum(unique(xx1$cdr3)%in%unique(xx2$cdr3))
    
    match_xx1<-match_xx1_xx2/total_xx1
    match_xx2<-match_xx1_xx2/total_xx2
    match_different_cluster[tag1,tag2]<-match_xx1
    match_different_cluster[tag2,tag1]<-match_xx2
    
    
  }
 

myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/matched_TCR/clonal_match_between_clusters/percentage_clonal_match_between_clusters.xls"
write.table(match_different_cluster,myoutf,quote=F,sep='\t')

library(ggcorrplot)
library(reshape2)
library(ggpubr)



reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}


myinf2<-"/dartfs-hpc/rc/lab/T/TurkLab/Jichang/bulkTCR/match/different_tissue_each_PT_match_number/"
files1<-list.files(myinf2)
files1<-files1[grep(".xls",files)]



breaklist<-seq(0,3,0.1)


  myoutf<-paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/matched_TCR/clonal_match_between_clusters/cor_plot_clonal_match_between_clusters.pdf")
  pdf(myoutf,width=6,height=5)
  print(corrplot(match_different_cluster, method="color",is.corr=F))
  
  dev.off()
 



















###########Do chi-square test for the matched clonotypes enrichment in different clusters###########



table(PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.3,PW_Adult_kit.integrated@meta.data$match)

xx<-matrix(0,10,2)
row.names(xx)<-seq(0,9,1)
colnames(xx)<-c("Not matched","Matched")
tar<-unique(PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.3)

for (i in 1:length(tar))
{
  xx[tar[i],"Not matched"]<-sum(PW_Adult_kit.integrated@meta.data[PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.3==tar[i],]$match=="Not Matched")
  xx[tar[i],"Matched"]<-sum(PW_Adult_kit.integrated@meta.data[PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.3==tar[i],]$match=="Skin Tumor Specifically Matched")+
    sum(PW_Adult_kit.integrated@meta.data[PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.3==tar[i],]$match=="Skin Tumor Blood Matched")
}


chi_square<-chisq.test(xx)

chi_observed<-chi_square$observed
chi_expected<-chi_square$expected
ROE<-chi_observed/chi_expected

myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/Seruat_PW_Adult_kit.integrated__finished.RDS"
saveRDS(PW_Adult_kit.integrated, file = myoutf)











myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/Seruat_PW_Adult_kit.integrated__finished.RDS"
PW_Adult_kit.integrated<-readRDS(myinf)

tar<-unique(PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.3)
tar<-tar[order(tar)]
enrich_all<-matrix(0,length(unique(PW_Adult_kit.integrated@meta.data$orig.ident))*length(tar),3)
colnames(enrich_all)<-c("Cluster","PT","enrich_index")
enrich_all[,"Cluster"]<-as.character(rep(tar,4))
enrich_all[,"PT"]<-as.character(sapply(unique(PW_Adult_kit.integrated@meta.data$orig.ident),function(x){rep(x,10)}))
enrich_all<-as.data.frame(enrich_all)
enrich_all$enrich_index<-as.numeric(levels(enrich_all$enrich_index))[enrich_all$enrich_index]

chi_all<-data.frame()
for(i in 1:length(unique(PW_Adult_kit.integrated@meta.data$orig.ident)))
{
   
   cat("\r",i)
   chi<-matrix(0,2,length(tar))
   PT<-unique(PW_Adult_kit.integrated@meta.data$orig.ident)[i]
   match<-c('matched','not matched')
   row.names(chi)<-paste0(PT,"_",match)
   colnames(chi)<-tar
   
   data<-PW_Adult_kit.integrated@meta.data[PW_Adult_kit.integrated@meta.data$orig.ident==unique(PW_Adult_kit.integrated@meta.data$orig.ident)[i],]
   
   for(j in 1:length(tar))
   {
      data_1<-data[data$integrated_snn_res.0.3==tar[j],]
      chi[1,tar[j]]<-sum(data_1$match=="Matched with TILs")
      chi[2,tar[j]]<-sum(data_1$match=="Not Matched")
      
   }
   #chi_all<-rbind(chi_all,chi)
   chi_square<-chisq.test(chi)
   
   chi_observed<-chi_square$observed
   chi_expected<-chi_square$expected
   ROE<-chi_observed/chi_expected
   
   for(k in 1:length(tar))
   {
      
      enrich_all[enrich_all$PT==PT,][k,"enrich_index"]<-ROE[1,tar[k]]
   }
   
}
   


tag<-enrich_all$Cluster=="0"|enrich_all$Cluster=="1"
enrich_zero_one<-enrich_all[tag,]


col<-brewer.pal(n = 10, name ="Paired")

myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/matched_TCR/Boxplot_matched_TCR_enrichment_level_cluster0&1.pdf"
pdf(myoutf,width=3,height=5)
print(ggplot(enrich_zero_one, aes(x=Cluster, y=enrich_index,fill=Cluster)) + 
         geom_boxplot()+
         geom_jitter(aes(shape = PT), position=position_jitter(0))+
         theme_classic()+
        scale_fill_manual(values=c("#a1d99b","#9ecae1")))
dev.off()


cluster0<-enrich_zero_one[enrich_zero_one$Cluster=="0","enrich_index"]
cluster1<-enrich_zero_one[enrich_zero_one$Cluster=="1","enrich_index"]


wilcox.test(cluster0,cluster1,alternative = "less")

#t.test(cluster0,cluster1,alternative = "less")

#Welch Two Sample t-test

#data:  cluster0 and cluster1
#t = -1.7744, df = 4.0918, p-value = 0.07453
#alternative hypothesis: true difference in means is less than 0
#95 percent confidence interval:
#  -Inf 0.1047329
#sample estimates:
#  mean of x mean of y 
#0.5553939 1.0959650 







































#############Clonality each cluster##################


rm(list=ls())


tar<-c(1,3,5,6)
barcode_info_all<-data.frame()

for(i in 1:length(tar))
{
  myinf1<-paste0("/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scTCRseq/Run_",tar[i],"/TCR/data/skin/skin_filtered_contig_annotations.csv")
  myinf2<-paste0("/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scTCRseq/Run_",tar[i],"/TCR/data/blood/blood_filtered_contig_annotations.csv")
  myinf3<-paste0("/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scTCRseq/Run_",tar[i],"/TCR/data/tumor/tumor_filtered_contig_annotations.csv")
  data_skin<-read.csv(myinf1)
  data_blood<-read.csv(myinf2)
  data_tumor<-read.csv(myinf3)
  data_skin[,"label"]<-rep("skin",nrow(data_skin))
  data_blood[,"label"]<-rep("blood",nrow(data_blood))
  data_tumor[,"label"]<-rep("tumor",nrow(data_tumor))
  
  data_skin[,"labeled_clonotype_id"]<-paste0(data_skin$raw_clonotype_id,"_",data_skin$label,tar[i])
  data_blood[,"labeled_clonotype_id"]<-paste0(data_blood$raw_clonotype_id,"_",data_blood$label,tar[i])
  data_tumor[,"labeled_clonotype_id"]<-paste0(data_tumor$raw_clonotype_id,"_",data_tumor$label,tar[i])
  
  barcode_info<-rbind(data_skin,data_blood,data_tumor)
  
  
  tag<-barcode_info$is_cell=="True"
  barcode_info<-barcode_info[tag,]
  tag<-barcode_info$high_confidence=="True"
  barcode_info<-barcode_info[tag,]
  tag<-barcode_info$productive=="True"
  barcode_info<-barcode_info[tag,]
  
  barcode_info$barcode<-gsub("-1","",barcode_info$barcode)
  barcode_info$barcode<-paste0(barcode_info$label,"_",barcode_info$barcode)
  barcode_info$barcode<-gsub("_",paste0("_",tar[i],"_"),barcode_info$barcode)
  barcode_info[,"PT"]<-rep(paste0("PT_",tar[i]),nrow(barcode_info))
  barcode_info_all<-rbind(barcode_info,barcode_info_all)
}

barcode_info_all_brief<-unique(barcode_info_all[,c("barcode","label","labeled_clonotype_id","PT")])



myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/Seruat_PW_Adult_kit.integrated__finished.RDS"
PW_Adult_kit.integrated<-readRDS(myinf)

RNAinfo<-PW_Adult_kit.integrated@meta.data


tag<-barcode_info_all_brief$barcode%in%row.names(RNAinfo)

barcode_info<-barcode_info_all_brief[tag,]
barcode_info[,"Cluster"]<-rep("0",nrow(barcode_info))


tar<-unique(PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.3)
tar<-tar[order(tar)]

for(i in 1:length(tar))
{
  barcode<-row.names(RNAinfo[RNAinfo$integrated_snn_res.0.3==tar[i],])
  tag<-barcode_info$barcode%in%barcode
  barcode_info$Cluster[tag]<-as.character(tar[i])
}

barcode_info[,"frequency"]<-rep(0,nrow(barcode_info))




tar<-unique(PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.3)
tar<-tar[order(tar)]


for(i in 1:length(tar))

{
  
  tar1<-unique(barcode_info[barcode_info$Cluster==tar[i],]$labeled_clonotype_id)
  for(j in 1:length(tar1))
  {
    barcode_info[barcode_info$Cluster==tar[i],][barcode_info[barcode_info$Cluster==tar[i],]$labeled_clonotype_id==tar1[j],"frequency"]<-sum(barcode_info[barcode_info$Cluster==tar[i],]$labeled_clonotype_id==tar1[j])
  }
  
}


ts<-unique(barcode_info[,c("label","labeled_clonotype_id","PT","Cluster","frequency")])


sc_clonality<-matrix(0,40,3)
colnames(sc_clonality)<-c("Cluster","PT","clonality_index")
sc_clonality<-as.data.frame(sc_clonality)

#sc_clonality$clonality_index<-as.numeric(levels(sc_clonality$clonality_index))[sc_clonality$clonality_index]
tar<-unique(PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.3)
tar<-tar[order(tar)]
sc_clonality$Cluster<-rep(tar,4)
sc_clonality$PT<-as.character(sapply(unique(PW_Adult_kit.integrated@meta.data$orig.ident),function(x){rep(x,10)}))


tar<-unique(PW_Adult_kit.integrated@meta.data$orig.ident)
tar1<-unique(ts$Cluster)

for(i in 1:length(tar))
{
  xx<-ts[ts$PT==tar[i],]
  
  for(j in 1:length(tar1))
  {
    yy<-xx[xx$Cluster==tar1[j],]
    sum_all<-sum(yy$frequency)
    
    index<-1-((sum(-(yy$frequency/sum_all)*log2(yy$frequency/sum_all)))/log2(nrow(yy)))
    
    sc_clonality[sc_clonality$PT==tar[i],][sc_clonality[sc_clonality$PT==tar[i],]$Cluster==tar1[j],"clonality_index"]<-index
  }
  
}


tag<-sc_clonality$Cluster=="0"|sc_clonality$Cluster=="1"
sc_clonality_zero_one<-sc_clonality[tag,]




col<-brewer.pal(n = 10, name ="Paired")

myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/matched_TCR/Boxplot_clonality_cluster0&1.pdf"
pdf(myoutf,width=3,height=5)
print(ggplot(sc_clonality_zero_one, aes(x=Cluster, y=clonality_index,fill=Cluster)) + 
        geom_boxplot()+
        geom_jitter(aes(shape = PT), position=position_jitter(0))+
        theme_classic()+
        scale_fill_manual(values=c("#a1d99b","#9ecae1")))
dev.off()


cluster0<-sc_clonality_zero_one[sc_clonality_zero_one$Cluster=="0","clonality_index"]
cluster1<-sc_clonality_zero_one[sc_clonality_zero_one$Cluster=="1","clonality_index"]

t.test(cluster0,cluster1,alternative = "greater")

#t.test(cluster0,cluster1,alternative = "greater")

#Welch Two Sample t-test

#data:  cluster0 and cluster1
#t = 0.64076, df = 5.9365, p-value = 0.2728
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  -0.04080973         Inf
#sample estimates:
#  mean of x  mean of y 
#0.09712415 0.07710393 






























































#########Gene list tumor-targeting Trm VS not tumor-targeting Trm#########
#   Not matched    Matched
#0  1.1158121 0.42392592
#1   0.9422768 1.28712768
#2   0.7585106 2.20122041
#3   1.0227658 0.88675799
#4   1.1648609 0.17994624
#5   0.6598938 2.69176156
#6   1.1949194 0.03042894
#7   1.1320951 0.34293078
#8   1.1970199 0.01998065
#9   1.1216293 0.39498944

myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/Seruat_PW_Adult_kit.integrated__finished.RDS"
PW_Adult_kit.integrated<-readRDS(myinf)



Idents(PW_Adult_kit.integrated)<-'integrated_snn_res.0.3'

tmp.markers<-FindMarkers(PW_Adult_kit.integrated,
                         ident.1=c('0'),logfc.threshold = 0.01)
myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/matched_TCR/non_tumor_targeting_Trm_list.xls"
write.table(tmp.markers,myoutf,sep="\t",quote=F)

tmp.markers<-FindMarkers(PW_Adult_kit.integrated,
                         ident.1=c('7'),logfc.threshold = 0.01)
myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/matched_TCR/Tex_list.xls"
write.table(tmp.markers,myoutf,sep="\t",quote=F)

tmp.markers<-FindMarkers(PW_Adult_kit.integrated,
                         ident.1=c('6'),logfc.threshold = 0.01)
myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/matched_TCR/Tgzmk_list.xls"
write.table(tmp.markers,myoutf,sep="\t",quote=F)



tmp.markers<-FindMarkers(PW_Adult_kit.integrated,
                         ident.1=c('2'),logfc.threshold = 0.01)
myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/matched_TCR/tumor_targeting_Tcir_list.xls"
write.table(tmp.markers,myoutf,sep="\t",quote=F)


tmp.markers<-FindMarkers(PW_Adult_kit.integrated,
                         ident.1=c('1'),logfc.threshold = 0.01)
myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/matched_TCR/tumor_targeting_Trm_list.xls"
write.table(tmp.markers,myoutf,sep="\t",quote=F)


tmp.markers<-FindMarkers(PW_Adult_kit.integrated,
                         ident.1=c('5'),logfc.threshold = 0.01)
myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/matched_TCR/tumor_targeting_Tcyto_list.xls"
write.table(tmp.markers,myoutf,sep="\t",quote=F)



tmp.markers<-FindMarkers(PW_Adult_kit.integrated,
                         ident.1=c('4',"8"),logfc.threshold = 0.01)
myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/matched_TCR/non_tumor_targeting_Tcirc_list.xls"
write.table(tmp.markers,myoutf,sep="\t",quote=F)

tmp.markers<-FindMarkers(PW_Adult_kit.integrated,
                         ident.1=c('3'),logfc.threshold = 0.01)
myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/matched_TCR/tumor_targeting_skin_non_trm_list.xls"
write.table(tmp.markers,myoutf,sep="\t",quote=F)

tmp.markers<-FindMarkers(PW_Adult_kit.integrated,
                         ident.1=c('0',"1","7"),logfc.threshold = 0.01)
myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/survival/signatures/all_trm_list.xls"
write.table(tmp.markers,myoutf,sep="\t",quote=F)

tmp.markers<-FindMarkers(PW_Adult_kit.integrated,
                         ident.1=c('0',"1"),logfc.threshold = 0.01)
myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/survival/signatures/all_trm_list_no_Tex.xls"
write.table(tmp.markers,myoutf,sep="\t",quote=F)


#########Gene list tumor-targeting Trm VS tumor targeting circulating#########

Idents(PW_Adult_kit.integrated)<-"match"
tumor_targeting<-subset(PW_Adult_kit.integrated,idents=c("Skin Tumor specifically matched clonotypes","Skin Tumor Blood matched clonotypes"))
table(tumor_targeting@meta.data$integrated_snn_res.0.3,tumor_targeting@meta.data$match)

Idents(PW_Adult_kit.integrated)<-'integrated_snn_res.0.3'

tmp.markers<-FindMarkers(PW_Adult_kit.integrated, ident.2 = c('1','5'),
                         ident.1='0',logfc.threshold = 0.01)
myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/matched_TCR/differentially_genes_tumor_targetingTrm_vs_tumortargeting_Tem.xls"
write.table(tmp.markers,myoutf,sep="\t",quote=F)

########Gene list ST specific VS STB matched tumor targeting in skin#####
Idents(tumor_targeting)<-"label"
skin_match<-subset(tumor_targeting,idents=c("skin"))

table(skin_match@meta.data$orig.ident,skin_match@meta.data$match)

#                Skin Tumor Blood matched clonotypes
#PT_1                                   5
#PT_3                                   3
#PT_5                                 122
#PT_6                                  48

#                        Skin Tumor specifically matched clonotypes
#PT_1                                         10
#PT_3                                          8
#PT_5                                        130
#PT_6                                         63





























###get z transformed data
myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/expression_data/CD8_normalized_data.xls"
data_ztransformed<-read.table(myinf,sep="\t",stringsAsFactors = FALSE,header=TRUE)
avg <- apply(data_ztransformed,1, mean)
tag <- avg > 0
data_ztransformed <- data_ztransformed[tag,]
res <- apply(data_ztransformed,1, function(arg) (arg-mean(arg))/sd(arg))
res <- t(res)
myoutf2 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/expression_data/PW_Adult_kit.integrated_ztransformed_normalized_data.xls"
write.table(res,myoutf2,sep="\t",quote=F)

















# Average and median expression of each cluster based on z transformed data ----------

myinf1<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/expression_data/PW_Adult_kit.integrated_ztransformed_normalized_data.xls"
data <- read.table(myinf1,sep="\t",quote=NULL)
myinf2<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/expression_data/CD8_normalized_cell_info.xls"
info<-read.table(myinf2,sep="\t",quote=NULL)

#Just get avg expression on selected genes

gene<-c("CD69","NR4A1","IFNG","TNF","GZMB","RUNX3","RGS1","JUN","JUNB","FOS","FOSB","JUND",
        "SELL","CCR7","TCF7","LEF1","KLF2","S1PR1","CX3CR1","KLRG1","NKG7","FGFBP2","PRF1","GNLY",
        "TOX","CXCL13","PDCD1","CTLA4","EOMES","LAG3","TIGIT","HAVCR2")


data<-data[row.names(data)%in%gene,]

tar <- unique(info$integrated_snn_res.0.3)
#tar <- as.numeric(tar)
tar <- tar[order(tar)]

for(i in 1 : length(tar))
{
  res <- matrix(0, nrow(data),4)
  row.names(res) <- row.names(data)
  colnames(res) <- c("C1_avg_expression","C1_median_expression","C2_avg_expression","C2_median_expression")
  
  tag <- info$integrated_snn_res.0.3 == tar[i]
  
  sam1 <- row.names(info)[tag==1]
  sam2 <- row.names(info)[tag==0]
  
  for(k in 1 : nrow(data))
  {
    cat("\r",i,"----->",k)
    xx1 <- as.numeric(data[k, sam1])
    xx2 <- as.numeric(data[k, sam2])
    
    res[k,"C1_avg_expression"] <- mean(xx1)
    res[k,"C2_avg_expression"] <- mean(xx2)
    res[k,"C1_median_expression"]<-median(xx1)
    res[k,"C2_median_expression"] <- median(xx2)
  }
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/Trm_signature/For_selected_genes_cluster_",tar[i],"_vs_others_ztransformed_expression.xls")
  write.table(res,myoutf,sep="\t",quote=F)
}


#Heatmap

res<-matrix(0,nrow(data),length(tar))
row.names(res)<-row.names(data)
colnames(res)<-tar[order(tar)]

for (i in 1:length(tar))
{
  myinf<-paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/Trm_signature/For_selected_genes_cluster_",tar[i],"_vs_others_ztransformed_expression.xls")
  xx<-read.table(myinf,sep='\t')
  res[,i]<-xx$C1_avg_expression
}

col.pal <- rev(RColorBrewer::brewer.pal(11, "Spectral"))
col.pal


myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/cluster_identity/heatmap_identify_everycluster.pdf")
pdf(myoutf,width=20,height=10)
pheatmap(mat=res,
         color             = col.pal,
         border_color      = NA,
         show_colnames     = T,
         show_rownames     = F,
         drop_levels       = TRUE,
         scale='none',
         main=paste0("Expression of Trm signature defined by Mackey Trm self modified signature"),
         fontsize=20)
dev.off()







#Get avg expression on all genes
tar <- unique(info$integrated_snn_res.0.3)
#tar <- as.numeric(tar)
tar <- tar[order(tar)]

for(i in 1 : length(tar))
{
  res <- matrix(0, nrow(data),4)
  row.names(res) <- row.names(data)
  colnames(res) <- c("C1_avg_expression","C1_median_expression","C2_avg_expression","C2_median_expression")
  
  tag <- info$integrated_snn_res.0.3 == tar[i]
  
  sam1 <- row.names(info)[tag==1]
  sam2 <- row.names(info)[tag==0]
  
  for(k in 1 : nrow(data))
  {
    cat("\r",i,"----->",k)
    xx1 <- as.numeric(data[k, sam1])
    xx2 <- as.numeric(data[k, sam2])
    
    res[k,"C1_avg_expression"] <- mean(xx1)
    res[k,"C2_avg_expression"] <- mean(xx2)
    res[k,"C1_median_expression"]<-median(xx1)
    res[k,"C2_median_expression"] <- median(xx2)
  }
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/Trm_signature/cluster_",tar[i],"_vs_others_ztransformed_expression.xls")
  write.table(res,myoutf,sep="\t",quote=F)
}






































tar <- unique(info$integrated_snn_res.0.3)
#tar <- as.numeric(tar)
tar <- tar[order(tar)]

for(i in 1 : length(tar))
{
  data <- matrix(0, nrow(res),4)
  row.names(data) <- row.names(res)
  colnames(data) <- c("C1_avg_expression","C1_median_expression","C2_avg_expression","C2_median_expression")
  
  tag <- info$integrated_snn_res.0.3 == tar[i]
  
  sam1 <- row.names(info)[tag==1]
  sam2 <- row.names(info)[tag==0]
  
  for(k in 1 : nrow(res))
  {
    cat("\r",i,"----->",k)
    xx1 <- as.numeric(res[k, sam1])
    xx2 <- as.numeric(res[k, sam2])
    
    data[k,"C1_avg_expression"] <- mean(xx1)
    data[k,"C2_avg_expression"] <- mean(xx2)
    data[k,"C1_median_expression"]<-median(xx1)
    data[k,"C2_median_expression"] <- median(xx2)
  }
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/GSEA/cluster_",tar[i],"_vs_others_normalized_expression.xls")
  write.table(data,myoutf,sep="\t",quote=F)
}


tar<-0:10
tmpinf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/GSEA/cluster_1_vs_others_normalized_expression.xls"
res <- read.table(tmpinf,sep="\t",quote=NULL)
GSEA_file<-matrix(0,nrow(res),11)
row.names(GSEA_file) <- row.names(res)
colnames(GSEA_file) <- c("cluster_0","cluster_1","cluster_2","cluster_3","cluster_4","cluster_5","cluster_6","cluster_7","cluster_8",
                         "cluster_9","cluster_10")

for(i in 1 : length(tar))
{
  
  tmpinf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/GSEA/cluster_",tar[i],"_vs_others_normalized_expression.xls")
  res <- read.table(tmpinf,sep="\t",quote=NULL)
  exp <- res$C1_avg_expression
  GSEA_file[,i]<-exp
}

myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/GSEA/PW_Adult_kit.integrated_noregress_res.0.3_GSEAinput.gct")
write.table(GSEA_file,file=myoutf,quote=F,sep="\t",row.names=T)





#info<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/GSEA/cluster_0_vs_others_normalized_expression.xls"
#ts<-read.table(info,sep="\t",quote=NULL)















































############Make TCR matched clone UMAP plot##########
myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/Seruat_PW_Adult_kit.integrated__finished.RDS"
PW_Adult_kit.integrated<-readRDS(myinf)


myinf1<-"/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scTCRseq/Run_6/TCR/data/skin/skin_clonotypes.csv"
myinf2<-"/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scTCRseq/Run_6/TCR/data/tumor/tumor_clonotypes.csv"
myinf3<-"/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scTCRseq/Run_6/TCR/data/tumor/tumor_clonotypes.csv"
data_skin<-read.csv(myinf1)
data_tumor<-read.csv(myinf2)
data_tumor<-read.csv(myinf3)
data_skin[,"label"]<-rep("skin",nrow(data_skin))
data_tumor[,"label"]<-rep("tumor",nrow(data_tumor))
data_tumor[,"label"]<-rep("tumor",nrow(data_tumor))
data_skin[,"clonotype_id"]<-paste0(data_skin$clonotype_id,"_",data_skin$label)
data_tumor[,"clonotype_id"]<-paste0(data_tumor$clonotype_id,"_",data_tumor$label)
data_tumor[,"clonotype_id"]<-paste0(data_tumor$clonotype_id,"_",data_tumor$label)

data_tumor<-data_tumor[data_tumor$frequency>2,]

data_skin_matched_expanded_tumor<-data_skin[data_skin$cdr3s_nt%in%data_tumor$cdr3s_nt,]
data_tumor_matched_expanded_tumor<-data_tumor[data_tumor$cdr3s_nt%in%data_tumor$cdr3s_nt,]
data_tumor_targeting<-rbind(data_skin_matched_expanded_tumor,
                            data_tumor_matched_expanded_tumor)
#data_match<-rbind(data_skin,data_tumor,data_tumor)




myinf4<-"/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scTCRseq/Run_6/TCR/data/skin/skin_filtered_contig_annotations.csv"
myinf5<-"/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scTCRseq/Run_6/TCR/data/tumor/tumor_filtered_contig_annotations.csv"
myinf6<-"/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scTCRseq/Run_6/TCR/data/tumor/tumor_filtered_contig_annotations.csv"

data1<-read.csv(myinf4)
data2<-read.csv(myinf5)
data3<-read.csv(myinf6)

data1[,"label"]<-rep("skin",nrow(data1))
data2[,"label"]<-rep("tumor",nrow(data2))
data3[,"label"]<-rep("tumor",nrow(data3))

data1[,"labeled_clonotype_id"]<-paste0(data1$raw_clonotype_id,"_",data1$label)
data2[,"labeled_clonotype_id"]<-paste0(data2$raw_clonotype_id,"_",data2$label)
data3[,"labeled_clonotype_id"]<-paste0(data3$raw_clonotype_id,"_",data3$label)

data<-rbind(data1,data2,data3)
tag<-data$is_cell=="True"
data<-data[tag,]
tag<-data$high_confidence=="True"
data<-data[tag,]
tag<-data$productive=="True"
data<-data[tag,]


table(data$label)

data$barcode<-gsub("-1","",data$barcode)
data$barcode<-paste0(data$label,"_",data$barcode)
data$barcode<-gsub("_","_6_",data$barcode)

barcode_matched<-data[data$labeled_clonotype_id%in%data_tumor_targeting$clonotype_id,"barcode"]
barcode_RNAseq<-row.names(PW_Adult_kit.integrated@meta.data)

tag<-(barcode_RNAseq%in%barcode_matched)

PW_Adult_kit.integrated@meta.data[,'TCR_matched_with_expanded_tumor']<-rep("Not matched",nrow(PW_Adult_kit.integrated@meta.data))
PW_Adult_kit.integrated@meta.data$TCR_matched_with_expanded_tumor[tag]="Cells with clonal expanded TIL clonotypes"



myoutf <- "/dartfs-hpc/rc/lab/T/TurkLab/Jichang/qual/UMAP_tumor_targeting_clonotypes.pdf"
pdf(myoutf,width=10,height=5)
DimPlot(object = PW_Adult_kit.integrated,group.by='TCR_matched_with_expanded_tumor',reduction="umap",cols=c("#de2d26",'#efeded'),
        cells=PW_Adult_kit.integrated@meta.data$orig.ident==c("skin_6","tumor_6"))
dev.off()
































































###############All clusters' marker genes' feature plot##########


myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/Seruat_PW_Adult_kit.integrated__finished.RDS"
PW_Adult_kit.integrated<-readRDS(myinf)

cluster_0_feature<-c("GZMK",'GZMA','NKG7','EOMES',
                     'LAG3','KLRB1','KLF2','S1PR1',
                     'SPON2','SELL','KLRG1')
cluster_1_feature<-c('IGKV3-20','IGKC','IGHG1','CXCR4',
                     'CCL4','DUSP2','BTG2','NR4A2',
                     'JUN','COTL1','CD69')
cluster_2_feature<-c('GZMK','CCL4','CCL4L2','IFNG',
                     'CD74','TNF','FOS','CD69','NR4A2',
                     'XCL2','CCL5','FOSB','CCL3','NR4A1',
                     'NR4A3','CD160','JUN')
#cluster_3_feature<-c('IL7R','TIMP1','MAL','FOS')

cluster_4_feature<-c('CXCL13','GEM','GK','CD4',
                     'CTLA4','HAVCR2','TIGIT','PDCD1')

#cluster_5_feature<-c('GZMK','RGS2','NR4A2','JUN',
'DNAJA1','CXCR4','RGS1','CCL5',
'SERTAD1','XCL1','FOSB','XCL2',
'IFNG')

#cluster6 same with cluster0

cluster_7_feature<-c('FOS','FXYD7','CD40LG','FOSB',
                     'NR4A1','IL7R','JUNB','CD40LG','C1orf162',
                     'ZNF683','CD28')

cluster_8_feature<-c('HAVCR2','CXCL13','CTLA4','LAYN',
                     'GEM','RGS1','DUSP4','TNFRSF9','LAG3',
                     'GZMB','CXCR6','TOX','PDCD1',
                     'ZNF331','ITGAE','CD69','IL2RB',
                     'CCL5','C10orf54')


#cluster_9_feature<-c('LTB','EEF1G','TRADD','NOSIP',
'NME2','PTPRCAP','KLF2','S1PR1','AQP3',
'TCF7','LEF1','SELL','IL7R')

#cluster_10_feature<-c('FGFBP2','CX3CR1','FGR','GZMH',
'S1PR5','GNLY','KLRD1','NKG7','HOPX',
'KLRG1','ABI3','CCL5')






cluster_Tex_feature<-c('SLAMF6','TCF7','CXCR3','ENTPD1','HAVCR2','MKI67','PCBD2','TBX21','CXCR5','TOX')

myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/feature_plot_res0.7/cluster0.pdf"
pdf(myoutf,width=22,height=16)
FeaturePlot(object = PW_Adult_kit.integrated, features= cluster_0_feature,cols = c("grey", "red"),
            reduction = "umap",ncol = 4,pt.size=1)
dev.off()

myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/feature_plot_res0.7/cluster1.pdf"
pdf(myoutf,width=22,height=12)
FeaturePlot(object = PW_Adult_kit.integrated, features= cluster_1_feature,cols = c("grey", "red"),
            reduction = "umap",ncol = 4,pt.size=1)
dev.off()

myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/feature_plot_res0.7/cluster2.pdf"
pdf(myoutf,width=22,height=16)
FeaturePlot(object = PW_Adult_kit.integrated, features= cluster_2_feature,cols = c("grey", "red"),
            reduction = "umap",ncol = 4,pt.size=1)
dev.off()

myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/feature_plot_res0.7/cluster3.pdf"
pdf(myoutf,width=22,height=12)
FeaturePlot(object = PW_Adult_kit.integrated, features= cluster_3_feature,cols = c("grey", "red"),
            reduction = "umap",ncol = 4,pt.size=1)
dev.off()

myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/feature_plot_res0.7/cluster4.pdf"
pdf(myoutf,width=22,height=12)
FeaturePlot(object = PW_Adult_kit.integrated, features= cluster_4_feature,cols = c("grey", "red"),
            reduction = "umap",ncol = 4,pt.size=1)
dev.off()

myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/feature_plot_res0.7/cluster5.pdf"
pdf(myoutf,width=22,height=16)
FeaturePlot(object = PW_Adult_kit.integrated, features= cluster_5_feature,cols = c("grey", "red"),
            reduction = "umap",ncol = 4,pt.size=1)
dev.off()

myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/feature_plot_res0.7/cluster6.pdf"
pdf(myoutf,width=22,height=16)
FeaturePlot(object = PW_Adult_kit.integrated, features= cluster_6_feature,cols = c("grey", "red"),
            reduction = "umap",ncol = 4,pt.size=1)
dev.off()

myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/feature_plot_res0.7/cluster7.pdf"
pdf(myoutf,width=22,height=12)
FeaturePlot(object = PW_Adult_kit.integrated, features= cluster_7_feature,cols = c("grey", "red"),
            reduction = "umap",ncol = 4,pt.size=1)
dev.off()

myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/feature_plot_res0.7/cluster8.pdf"
pdf(myoutf,width=22,height=20)
FeaturePlot(object = PW_Adult_kit.integrated, features= cluster_8_feature,cols = c("grey", "red"),
            reduction = "umap",ncol = 4,pt.size=1)
dev.off()

myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/feature_plot_res0.7/cluster9.pdf"
pdf(myoutf,width=22,height=12)
FeaturePlot(object = PW_Adult_kit.integrated, features= cluster_9_feature,cols = c("grey", "red"),
            reduction = "umap",ncol = 4,pt.size=1)
dev.off()

myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/feature_plot_res0.7/cluster10.pdf"
pdf(myoutf,width=22,height=12)
FeaturePlot(object = PW_Adult_kit.integrated, features= cluster_10_feature,cols = c("grey", "red"),
            reduction = "umap",ncol = 4,pt.size=1)
dev.off()

myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/feature_plot_res0.7/cluster11.pdf"
pdf(myoutf,width=22,height=4)
FeaturePlot(object = PW_Adult_kit.integrated, features= cluster_11_feature,cols = c("grey", "red"),
            reduction = "umap",ncol = 4,pt.size=1)
dev.off()

myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/feature_plot_res0.7/cluster12.pdf"
pdf(myoutf,width=22,height=16)
FeaturePlot(object = PW_Adult_kit.integrated, features= cluster_12_feature,cols = c("grey", "red"),
            reduction = "umap",ncol = 4,pt.size=1)
dev.off()

myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/feature_plot_res0.7/cluster13.pdf"
pdf(myoutf,width=22,height=8)
FeaturePlot(object = PW_Adult_kit.integrated, features= cluster_13_feature,cols = c("grey", "red"),
            reduction = "umap",ncol = 4,pt.size=1)
dev.off()







































































#[2]Get the moSTB representative samples
rm(list=ls())
myinf1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/expression_data/CD8_normalized_data.xls"
myinf2 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/expression_data/CD8_normalized_cell_info.xls"

data <- read.table(myinf1,sep="\t",quote=NULL)
info <- read.table(myinf2,sep="\t",quote=NULL)

#Judge by fractons tumor_5,tumor_5,skin_5


tag1 <- grep("tumor_5", info$orig.ident)
tag2 <- grep("tumor_5",info$orig.ident)
tag3 <- grep("skin_5",info$orig.ident)

info1 <- info[tag1,]
info2 <- info[tag2,]
info3 <- info[tag3,]
info <- rbind(info1,info2,info3)

com <- intersect(row.names(info),colnames(data))

data <- data[,com]

myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/expression_data/rep/CD8_normalized_data.xls"
write.table(data,myoutf,sep="\t",quote=F)
myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/expression_data/rep/CD8_normalized_data.xls"
teSTB<-read.table(myinf,sep="\t",STBringsAsFactors = FALSE,header=TRUE)






# FGSEA noregress ---------------------------------------------------------


myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/Seruat_PW_Adult_kit.integrated__finished.RDS"
PW_Adult_kit.integrated<-readRDS(myinf)
Idents(object = PW_Adult_kit.integrated) <- "integrated_snn_res.0.3"
tar<-unique(PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.3)
tar<-tar[order(tar)]
for(i in 1:length(tar))
{
  cluster.markers <- FindMarkers(PW_Adult_kit.integrated, ident.1 =tar[i],logfc.threshold = 0.01,min.pct = 0.01)
  myoutf<-paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/GSEA/marker/cluster_",tar[i],"differentially_expressed_markers.xls")
  write.table(cluster.markers,myoutf,quote=F,sep="\t")
}


filename<-"/dartfs-hpc/rc/lab/T/TurkLab/Jichang/pub_data/c2.all.v6.0.symbols.gmt.txt"

data = file(filename,open="r")
n = 1
pathway = list()

while(TRUE){
  line = readLines(data,n=1)
  if(length(line) == 0){
    break
  }
  pathway[n] = line
  n = n+1
}
close(data)


all_pathway = list()
names = c()
for(i in 1:length(pathway)){
  one_pathway = STBrsplit(as.character(pathway[i]),"\t")[[1]]
  names[i] = one_pathway[1]
  one_pathway = one_pathway[-1]
  one_pathway = one_pathway[-1]
  all_pathway[[i]] = one_pathway
}
names(all_pathway) = names

dir1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/GSEA/"
files1 <- list.files(dir1)
dir<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/GSEA/marker/"
files <- list.files(dir)
file_nam <- gsub(".xls","",files)
file_nam1<-gsub(".xls","",files1)
files1<-files1[order(files1)]
files<-files[order(files)]

for(i in 1 : length(files))
{
  cat("\r",i)
  tmpinf1 <- paste0(dir1,files1[i])
  tmpinf<-paste0(dir,files[i])
  res <- read.table(tmpinf,sep="\t",quote=NULL)
  res1<- read.table(tmpinf1,sep="\t",quote=NULL)
  
  tag<-res1$C1_avg_expression>0.02
  nm<-row.names(res1)[tag]
  
  res<-res[nm,]
  
  
  res <- res[order(res$avg_logFC,decreasing=T),]
  
  xx <- res$avg_logFC
  names(xx) <- row.names(res)
  names(xx) <- toupper(names(xx))
  tag = duplicated(names(xx))
  xx = xx[tag==0]
  data <- fgsea(pathways = all_pathway, 
                STBats = xx,
                minSize=15,
                maxSize=500,
                nperm=10000)
  
  
  tag <- data$pval < 0.05
  data <- data[tag,]
  data <- data[order(data$NES,decreasing=T),]
  data <- as.data.frame(data)
  data <- data[,1:7]
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/GSEA/fgsea/fgsea_cluster_compare_",file_nam[i],"_GSEA.xls")
  write.table(data,myoutf,sep="\t",quote=F)
}





































# Name each cluster -------------------------------------------------------
current.cluster.ids <- c("ExhauSTBed Trm-like CD8 T cells", "Naive-like CD8 T cells","Central memory CD8 T cells","Effector memory CD8 T cells", "Effector CD8 T cells",
                         "mucosal-associated invariant T cells", "Tissue resident memory CD8")
new.cluster.ids<- c("ExhauSTBed Trm-like CD8 T cells", "Naive-like CD8 T cells","Central memory CD8 T cells","Effector memory CD8 T cells", "Effector CD8 T cells",
                    "Mucosal-associated invariant T cells", "Tissue resident memory CD8")
PW_Adult_kit.integrated@ident <- plyr::mapvalues(x = PW_Adult_kit.integrated@ident, from = current.cluster.ids, to = new.cluster.ids)
myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human//PCx/cluster_defined/TSNE_defined.pdf"
pdf(myoutf,width=8,height=5)
TSNEPlot(object = PW_Adult_kit.integrated, pt.size = 0.7,label.size = 2)
dev.off()

PW_Adult_kit.integrated <- stashIdent(object = PW_Adult_kit.integrated, save.name = "JichangDefined")
myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/cluster_defined/Seruat_PW_Adult_kit.integrated_finished_cluster_defined.RDS"
saveRDS(PW_Adult_kit.integrated, file = myoutf)




# stacked plot show tissue contribution to each cluster -------------------


###STBacked plot show tissue contribution to each cluster
myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/Seruat_PW_Adult_kit.integrated__finished.RDS"
PW_Adult_kit.integrated<-readRDS(myinf)

Idents(PW_Adult_kit.integrated)<-"integrated_snn_res.0.3"

data<-table(PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.3,PW_Adult_kit.integrated@meta.data$label)
#     tumor_5 skin_5 tumor_5
#0       0     19     786
#1     573      0     215
#2       0      6     439
#3     362      5       3
#4      34      3     328
#5     117    145      86
#6     233      3       4
#7      53      1       7
#8       0     21       0
data<-as.data.frame(data)
colnames(data)<-c('cluster','tissue','cell_number')
data[,'fraction']<-0
contribution<-matrix(0,9,3)
row.names(contribution)<-c('0','1','2','3','4','5','6','7','8')
colnames(contribution)<-c('tumor','tumor','skin')
tar1<-c('tumor_5','tumor_5','skin_5')
tar<-c('0','1','2','3','4','5','6','7','8')
for (i in 1:length(tar))
{
  tag<-grep(tar[i],data$cluster)
  xx<-data[tag,]
  totalcell<-as.numeric(sum(xx$cell_number))
  for (j in 1:length(tar1))
  {
    data$fraction[tag[j]]<-xx$cell_number[j]/totalcell
    
  }
}

myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/tissue_contribution/tissue_contribution_each_cluster.xls"
write.table(data,myoutf,sep="\t",quote=F)

myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/tissue_contribution/STBackplot_tissue_contribution_each_cluster.pdf"
pdf(myoutf,width=8,height=8)
ggplot() + geom_bar(aes(y = fraction, x = cluster, fill = tissue), data = data,
                    STBat="identity")+theme(axis.text.x = element_text(angle = 60, hjuSTB = 1,size=10),
                                            axis.title.x = element_blank())

dev.off()

# Go enrichment analysis --------------------------------------------------


# GO enrichment analysis --------------------------------------------------


###Go enrichment analysis


tar1<-seq(0,9,1)
Idents(PW_Adult_kit.integrated)<-"integrated_snn_res.0.3"


for(i in 1 : length(tar1))
{
  cat("\r",i)
  tmp <- PW_Adult_kit.integrated
  tmp.markers <- FindMarkers(object = tmp, ident.1=tar1[i],min.pct = 0.2, logfc.threshold = 0.2)
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/GO/cluster_markers_",tar1[i],"vs_others.xls")
  write.table(tmp.markers,myoutf,sep="\t",quote=F)
}

myinf1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/GO/cluster_markers_1vs_others.xls"
myinf2 <- "/dartfs-hpc/rc/lab/T/TurkLab/Jichang/pub_data/homo_Go.gaf"
myinf3 <- "/dartfs-hpc/rc/lab/T/TurkLab/Jichang/pub_data/Goterms.txt"
myinf4 <- "/dartfs-hpc/rc/lab/T/TurkLab/Jichang/pub_data/Go_domain.txt"

myoutf1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/GO/cluster0_up_gene_enrichment.xls"
myoutf2 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/GO/cluster0_dn_gene_enrichment.xls"

res <- read.table(myinf1,sep="\t",quote=NULL,header=TRUE)
res[,"HGNC.symbol"] <- toupper(row.names(res))


info <- read.delim(myinf2,sep="\t",skip=34,header= F,stringsAsFactors = FALSE)
terms <- read.table(myinf3,sep="\t",stringsAsFactors = FALSE)
go <- read.table(myinf4,sep="\t",quote=NULL,stringsAsFactors=F,header=T)


colnames(info)[3] <- "HGNC.symbol"
colnames(info)[5] <- "GO.Term.Accession"

tag_up <- res$avg_logFC > 0
tag_dn <- res$avg_logFC < 0

dat1 <-res[tag_up==1,]#114
dat1[,"gene"] <- dat1[,"HGNC.symbol"]
row.names(dat1) <- seq(1,nrow(dat1),1)

dat2 <- res[tag_dn==1,]#79
dat2[,"gene"] <- dat2[,"HGNC.symbol"]
row.names(dat2) <- seq(1,nrow(dat2),1)

useinfo <-c("GO.Term.Accession","HGNC.symbol")
info <- info[,useinfo]
terms[,"GO.Term.Accession"]<- row.names(terms)
row.names(terms) <- seq(1,nrow(terms),1)

info <- merge(info,terms,by.x="GO.Term.Accession",by.y="GO.Term.Accession")
colnames(info)[3] <- "GO.Term.Name"
allgene <- unique(info[,"HGNC.symbol"]) #19050
up_info <- merge(dat1,info,by.x="gene",by.y="HGNC.symbol")

up_domain <- unique(up_info[,"GO.Term.Name"]) #1457
up_domain <- up_domain[up_domain!=""]
res <- matrix(nrow=length(up_domain),ncol=5)
colnames(res) <- c("enrichment","p.value","odds.ratio","shared_gene","all_gene")

for (i in 1 : length(up_domain))
{
  cat("\r",i)
  tag <- up_info[,"GO.Term.Name"] == up_domain[i]
  tmp <- up_info
  useinfo<- c("gene","GO.Term.Name")
  tmp <- tmp[,useinfo]
  tmp[tag==0,"GO.Term.Name"] <- "Others"
  tmp[,"GO.Term.Name"] <- factor(tmp[,"GO.Term.Name"])
  tmp_1<- tmp[tag==1,]
  tmp_1 <- unique(tmp_1)
  i_gene <- tmp_1[,"gene"]
  xx1 <-length(i_gene)
  tmp_2 <- tmp[tag==0,]
  tmp_2 <- unique(tmp_2)
  tag_r <- is.na(pmatch(tmp_2[,"gene"],i_gene))
  tmp_2 <- tmp_2[tag_r==1,]
  tmp <- rbind(tmp_1,tmp_2)
  xx2 <- nrow(tmp)
  white <- info[which(info[,"GO.Term.Name"] %in% up_domain[i]),"HGNC.symbol"]
  white <- unique(white)
  black <- is.na(pmatch(allgene,white))
  black <- unique(allgene[black])
  white <- length(white)
  black <- length(black)
  ER = (xx1/xx2)/(white/(white+black))
  myp = sum(dhyper(xx1:xx2, white, black, xx2))
  res[i,"p.value"] <- myp
  res[i,"odds.ratio"] <- ER
  res[i,"shared_gene"] <- xx1
  res[i,"all_gene"] <- white
  res[i,"enrichment"] <-up_domain[i]
  
  
}
cutoff <- 10
tag_cut <- as.numeric(res[,"shared_gene"]) > 5
res <- res[tag_cut==1,]
p.value <- res[,"p.value"]
p.value <- as.numeric(p.value)
p.adjust <- p.adjust(p.value,method = "BH")
res <- as.data.frame(res)
res[,"p.adjust"] <- p.adjust
res <- as.matrix(res)
rank <- order(as.numeric(res[,"p.value"]),decreasing=FALSE)
res <- res[rank,]
write.table(res,myoutf1,sep="\t", quote=F)


down_info <- merge(dat2,info,by.x="gene",by.y="HGNC.symbol")
down_domain <- unique(down_info[,"GO.Term.Name"])
down_domain <- down_domain[down_domain!=""] #9328
res <- matrix(nrow=length(down_domain),ncol=5)
colnames(res) <- c("enrichment","p.value","odds.ratio","shared_gene","all_gene")
value1 <- rep(length(down_domain),0)
value2 <- rep(length(down_domain),0)
for (i in 1 : length(down_domain))
{
  cat("\r",i)
  tag <- down_info[,"GO.Term.Name"] == down_domain[i]
  tmp <- down_info
  useinfo<- c("gene","GO.Term.Name")
  tmp <- tmp[,useinfo]
  tmp[tag==0,"GO.Term.Name"] <- "Others"
  tmp[,"GO.Term.Name"] <- factor(tmp[,"GO.Term.Name"])
  tmp_1<- tmp[tag==1,]
  tmp_1 <- unique(tmp_1)
  i_gene <- tmp_1[,"gene"]
  xx1 <-length(i_gene)
  tmp_2 <- tmp[tag==0,]
  tmp_2 <- unique(tmp_2)
  tag_r <- is.na(pmatch(tmp_2[,"gene"],i_gene))
  tmp_2 <- tmp_2[tag_r==1,]
  tmp <- rbind(tmp_1,tmp_2)
  xx2 <- nrow(tmp)
  white <- info[which(info[,"GO.Term.Name"] %in% down_domain[i]),"HGNC.symbol"]
  white <- unique(white)
  black <- is.na(pmatch(allgene,white))
  black <- unique(allgene[black])
  white <- length(white)
  black <- length(black)
  ER = (xx1/xx2)/(white/(white+black))
  myp = sum(dhyper(xx1:xx2, white, black, xx2))
  res[i,"p.value"] <- myp
  res[i,"odds.ratio"] <- ER
  res[i,"shared_gene"] <- xx1
  res[i,"all_gene"] <- white
  res[i,"enrichment"] <-down_domain[i]
  value1[i] <- xx1
  value2[i] <- white
}

cutoff <- 10
tag_cut <- as.numeric(res[,"shared_gene"]) > 5
res <- res[tag_cut==1,]
p.value <- res[,"p.value"]
p.value <- as.numeric(p.value)
p.adjust <- p.adjust(p.value,method = "BH")
res <- as.data.frame(res)
res[,"p.adjust"] <- p.adjust
res <- as.matrix(res)
rank <- order(as.numeric(res[,"p.value"]),decreasing=FALSE)
res <- res[rank,]

write.table(res,myoutf2,sep="\t", quote=F)

###Heatmap of GO terms

myinf1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/GO/cluster6_up_gene_enrichment.xls"
myinf2 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/GO/cluster6_dn_gene_enrichment.xls"
up_pathway<-read.table(myinf1,sep="\t",quote=NULL)
dn_pathway<-read.table(myinf2,sep="\t",quote=NULL)
tag1<-up_pathway$p.adjuSTB<0.05
up_pathway<-up_pathway[tag1,]
tag2<-dn_pathway$p.adjuSTB<0.05
dn_pathway<-dn_pathway[tag2,]
up_pathway<-up_pathway %>% select("enrichment","odds.ratio")
dn_pathway<-dn_pathway %>% select("enrichment","odds.ratio")

xx1 <- up_pathway[,-1]
xx1<-as.matrix(xx1)
row.names(xx1)<-up_pathway$enrichment
colnames(xx1)<-'Enrichment_score'
xx1<-cbind(xx1,xx1)



col.pal <- rev(RColorBrewer::brewer.pal(11, "Spectral"))
col.pal

myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/GO/heatmap_up_pathways.pdf"
pdf(myoutf,width=10,height=50)
pheatmap(mat=xx1,
         color             = col.pal,
         border_color      = NA,
         scale='column',
         main="up regulated genes enriched GO terms",
         show_colnames = T,
         show_rownames = T,
         cluster_row=T,
         clutser_col=FALSE,
         fontsize=12)
dev.off()


# Trm gene signature and GSVA ---------------------------------------------


###Trm_gene_signature and GSVA###


# Do Trm gene signature heatmap using Trm signatures -----------------



#AKM D45 skin signature prepare
myinf<-"/dartfs-hpc/rc/lab/T/TurkLab/Jichang/pub_data/Aleksey_D45_Skin_Trm.txt"
AKM_D45<-read.table(myinf,sep="\t",STBringsAsFactors = FALSE,header=TRUE)
tag<-AKM_D45$avg_logFC>0
AKM_UP<-row.names(AKM_D45[tag,])
AKM_DN<-row.names(AKM_D45[!tag,])
AKM_UP<-toupper(AKM_UP)
AKM_DN<-toupper(AKM_DN)
AKM_UP<-append(AKM_UP,rep("",103))








#Other Trm signature prepare
myinf<-"/dartfs-hpc/rc/lab/T/TurkLab/Jichang/pub_data/Trm_signature.txt"
Trm<-read.table(myinf,sep='\t',STBringsAsFactors = FALSE,header=TRUE)
#Trm[,"AKM_D45_SKIN_CD8_TRM_UP"]<-AKM_UP

#MACKAY_UP<-Trm[,"MACKAY_SKIN_TRM_UP"]
#MACKAY_DN<-Trm[,"MACKAY_SKIN_TRM_DN"]
#ZEMIN_ZNF683<-Trm[,"ZEMIN_CD8_C5_ZNF683"]
#ZEMIN_LAYN<-Trm[,"ZEMIN_CD8_C6_LAYN"]
#BC_UP<-Trm[,"SINGLECELL_BC_TRM"]

tar<-0:10

for (j in 1:ncol(Trm))
{
  gene<-Trm[,j]
  gene<-gene[gene != ""]
  res <- matrix(0, length(gene),11)
  row.names(res) <- gene
  colnames(res) <- c("cluster_0","cluster_1","cluster_2","cluster_3","cluster_4","cluster_5","cluster_6","cluster_7","cluster_8",
                     "cluster_9","cluster_10")
  for (i in 1:length(tar))
  {
    myinf1<-paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/Trm_signature/cluster_",tar[i],"_vs_others_ztransformed_expression.xls")
    data<-read.table(myinf1,sep="\t",header=TRUE,quote = NULL)
    data_use<-data[gene,]$C1_avg_expression
    res[,i]<-data_use
    myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/Trm_signature/Trm_signature",colnames(Trm)[j],"expression_everycluster.xls")
    write.table(res,myoutf,sep="\t",quote=F)
  }
}


#Define Trm score of each cluster using all Trm signatures
score<-matrix(0,ncol(Trm),11)
row.names(score) <- colnames(Trm)
colnames(score) <- c("cluster_0","cluster_1","cluster_2","cluster_3","cluster_4","cluster_5","cluster_6","cluster_7","cluster_8",
                     "cluster_9","cluster_10")
for (i in 1:ncol(Trm))
{
  myinf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/Trm_signature/Trm_signature",colnames(Trm)[i],"expression_everycluster.xls")
  res<-read.table(myinf,sep='\t',quote=NULL)
  tag<-!is.na(res[,1])
  res<-res[tag,]
  avg<-apply(res,2,sum)
  score[colnames(Trm)[i],]<-avg
}

xx1<-matrix(0,2,11)
row.names(xx1) <- c("AKM_score","Mackay_score")
colnames(xx1) <- c("cluster_0","cluster_1","cluster_2","cluster_3","cluster_4","cluster_5","cluster_6","cluster_7","cluster_8",
                   "cluster_9","cluster_10")
xx1["AKM_score",]<-as.numeric(score["AKM_D45_SKIN_CD8_TRM_UP",]-score["AKM_D45_SKIN_CD8_TRM_DN",])
xx1["Mackay_score",]<-as.numeric(score["MACKAY_SKIN_TRM_UP",]-score["MACKAY_SKIN_TRM_DN",])
score<-rbind(score,xx1)

myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/Trm_signature/Trm_score_all_signatures.xls")
write.table(score,myoutf,sep="\t",quote=F)






myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/expression_data/PW_Adult_kit.integrated_ztransformed_normalized_data.xls"
data<-read.table(myinf,sep='\t',STBringsAsFactors = FALSE,header=TRUE)
myinf2 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/expression_data/CD8_normalized_cell_info.xls"
info<-read.table(myinf2,sep="\t",quote=NULL)
myinf3<-"/dartfs-hpc/rc/lab/T/TurkLab/Jichang/pub_data/Trm_signature.txt"
Trm<-read.table(myinf3,sep='\t',STBringsAsFactors = FALSE,header=TRUE)
tag1<-row.names(data) %in% Trm$AKM_D45_SKIN_CD8_TRM_UP
tag2<-row.names(data) %in% Trm$AKM_D45_SKIN_CD8_TRM_DN
tag3<-row.names(data) %in% Trm$MACKAY_SKIN_TRM_UP
tag4<-row.names(data) %in% Trm$MACKAY_SKIN_TRM_DN

data_score<-matrix(0,7,ncol(data))
colnames(data_score)<-colnames(data)
row.names(data_score)<-c("AKM_UP","AKM_DN","MAC_UP","MAC_DN","AKM_score","MACKAY_score","cluster")
data_score["AKM_UP",]<-apply(data[tag1,],2,sum)
data_score["AKM_DN",]<-apply(data[tag2,],2,sum)
data_score["MAC_UP",]<-apply(data[tag3,],2,sum)
data_score["MAC_DN",]<-apply(data[tag4,],2,sum)
data_score["AKM_score",]<-data_score["AKM_UP",]-data_score["AKM_DN",]
data_score["MACKAY_score",]<-data_score["MAC_UP",]-data_score["MAC_DN",]
tar<-unique(info$integrated_snn_res.0.3)
tar<-tar[order(tar)]
for(i in 1:length(tar))
{
  tag<-info$integrated_snn_res.0.3==tar[i]
  nm<-row.names(info[tag,])
  tag1<-colnames(data_score) %in% nm
  data_score["cluster",tag1]<-tar[i]
}

myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/Trm_signature/SCORE_AKM_Trm_signature_boxplot.pdf"
pdf(myoutf,15,5)
boxplot(data_score["AKM_score",]~data_score["cluster",],data_score)
dev.off()

myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/Trm_signature/SCORE_AKM_UP_Trm_signature_boxplot.pdf"
pdf(myoutf,15,5)
boxplot(data_score["AKM_UP",]~data_score["cluster",],data_score)
dev.off()

myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/Trm_signature/SCORE_MACKAY_Trm_signature_boxplot.pdf"
pdf(myoutf,15,5)
boxplot(data_score["MACKAY_score",]~data_score["cluster",],data_score)
dev.off()

myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/Trm_signature/SCORE_MACKAY_UP_Trm_signature_boxplot.pdf"
pdf(myoutf,15,5)
boxplot(data_score["MAC_UP",]~data_score["cluster",],data_score)
dev.off()











for (j in 1:ncol(Trm))
{
  myinf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/Trm_signature/Trm_signature",colnames(Trm)[j],"expression_everycluster.xls")
  res<-read.table(myinf,sep='\t',quote=NULL)
  tag<-!is.na(res[,1])
  res<-res[tag,]
  
  
  col.pal <- rev(RColorBrewer::brewer.pal(11, "Spectral"))
  col.pal
  
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/Trm_signature/heatmap_signature_",colnames(Trm)[j],"_everycluster.pdf")
  pdf(myoutf,width=10,height=50)
  pheatmap(mat=res,
           color             = col.pal,
           border_color      = NA,
           show_colnames     = T,
           show_rownames     = F,
           drop_levels       = TRUE,
           scale='none',
           main=paste0("Expression of Trm signature defined by", colnames(Trm)[j],"signature"),
           fontsize=20)
  dev.off()
  
}
















# GSVA score ------------------------------------------------------------
rm(list=ls())

#[3]Discovery7 (pathway activity calculation)
command:
  #qsub -I -l nodes=1:ppn=20 -l walltime=200:00:00 -l feature=n05  
  #qsub -I -l nodes=1:ppn=40 -l walltime=12:00:00 -l feature=celln
  #/usr/bin/R --vanilla <<EOF
  library(GSVA)
filename = "/dartfs-hpc/rc/lab/T/TurkLab/Jichang/pub_data/c2.all.v6.0.symbols.gmt.txt"
data = file(filename,open="r")
n = 1
pathway = list()

while(TRUE){
  line = readLines(data,n=1)
  if(length(line) == 0){
    break
  }
  pathway[n] = line
  n = n+1
}
close(data)


all_pathway = list()
names = c()
for(i in 1:length(pathway)){
  one_pathway = STBrsplit(as.character(pathway[i]),"\t")[[1]]
  names[i] = one_pathway[1]
  one_pathway = one_pathway[-1]
  one_pathway = one_pathway[-1]
  all_pathway[[i]] = one_pathway
}
names(all_pathway) = names

#get the given amount of pathways

#Set up the data set for running
myinf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/expression_data/rep/CD8_normalized_data.xls"
res <- read.table(myinf,sep="\t",quote=NULL)

#screening out genes
var <- apply(res, 1, function(arg) var(arg))
tag <- var >0
res <- res[tag,]
tag <- row.names(res) == ""
res <- res[tag==0,]

#Make the gene names identical
gene <- row.names(res)
gene_num <- rep(0,length(all_pathway))

for(i in 1 : length(all_pathway))
{
  cat("\r",i)
  
  com <- intersect(as.vector(unlist(all_pathway[[i]])), gene)
  all_pathway[[i]] <- com
  gene_num[i] <- length(com)
}

tag <- gene_num > 20
all_pathway <- all_pathway[tag]
gene_num <- gene_num[tag]

#teSTB <- all_pathway[3]
#gsva_es <- gsva(as.matrix(res), gset.idx.list=teSTB, method="gsva")

#tmp.res<-res[,1:5]
#gsva_es <- gsva(as.matrix(tmp.res), gset.idx.list=all_pathway, method="gsva")




myoutf <-  "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/GSVA/SB3_noregress_pathway_activity.xls"
gsva_es <- gsva(as.matrix(res), gset.idx.list=all_pathway, method="gsva")
write.table(gsva_es,file = myoutf,row.names = TRUE,col.names = TRUE,sep = "\t",quote=F)


#[4]Calculate the p value for each cluster
myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/GSVA/SB3_noregress_pathway_activity.xls"
teSTB<-read.table(myinf,sep='\t',header=TRUE)

# calculate P value for each cluster --------------------------------------
rm(list=ls())
myinf1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/GSVA/SB3_noregress_pathway_activity.xls"
myinf2 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/Seruat_PW_Adult_kit.integrated_noregress_finished.RDS"

data <- read.table(myinf1,sep="\t",quote=NULL)
info <- readRDS(myinf2)

info <- info@meta.data

tar <- unique(info$res.0.9)
tar <- as.numeric(tar)
tar <- tar[order(tar)]

for(i in 1 : length(tar))
{
  res <- matrix(0, nrow(data),5)
  row.names(res) <- row.names(data)
  colnames(res) <- c("C1_path_score","C2_path_score","T.score","T.pval","W.pval")
  
  tag <- info$res.0.9 == tar[i]
  
  sam1 <- row.names(info)[tag==1]
  sam2 <- row.names(info)[tag==0]
  
  for(k in 1 : nrow(data))
  {
    cat("\r",i,"----->",k)
    xx1 <- as.numeric(data[k, sam1])
    xx2 <- as.numeric(data[k, sam2])
    
    fit1 <- t.teSTB(xx1,xx2)
    fit2 <- wilcox.teSTB(xx1,xx2)
    
    res[k,"C1_path_score"] <- mean(xx1)
    res[k,"C2_path_score"] <- mean(xx2)
    
    res[k,"W.pval"] <- wilcox.teSTB(xx1, xx2)$p.value
    res[k,"T.score"] <- t.teSTB(xx1, xx2)$STBatiSTBic
    res[k,"T.pval"] <- t.teSTB(xx1, xx2)$p.value
  }
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.3/GSVA/cluster_",tar[i],"_vs_others.xls")
  write.table(res,myoutf,sep="\t",quote=F)
}





# do GSVA heatmap --------------------------------------------------------------
myinf<-"/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/combine/RNA/vitiligo/STB1356/tumor_523/CD8/clean/PCx/GSVA/GSVA_comparison/cluster_2_vs_others.xls"
data<-read.table(myinf,sep="\t",header=TRUE,quote = NULL)
data<-data[order(data$T.pval),]
data_moSTB_sig<-data[1:30,]
path<-row.names(data_moSTB_sig)
tar<-0:5
res <- matrix(0, length(path),6)
row.names(res) <- path
colnames(res) <- c("tumor_cluster_0","tumor_cluster_1","tumor_cluster_2","tumor_cluster_3","tumor_cluster_4","tumor_cluster_5")
for (i in 1:length(tar))
{
  myinf3<-paste0("/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/combine/RNA/vitiligo/STB1356/tumor_523/CD8/clean/PCx/GSVA/GSVA_comparison/cluster_",tar[i],"_vs_others.xls")
  data<-read.table(myinf3,sep="\t",header=TRUE,quote = NULL)
  data_use<-data[path,]$C1_path_score
  res[,i]<-data_use
}
myoutf <- paste0("/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/combine/RNA/vitiligo/STB1356/tumor_523/CD8/clean/PCx/GSVA/GSVA_comparison/heatmap_c1path_score.xls")
write.table(res,myoutf,sep="\t",quote=F)
res<-res[c(1:3,5:30),]
res_all<-cbind(res,res_tumor)
myoutf <- "/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/combine/RNA/vitiligo/STB1356/tumor_523/CD8/clean/PCx/GSVA/GSVA_comparison/heatmap_compare_tumor_tumor_523_moSTBsig_30path.pdf"
pdf(myoutf,width=30,height=20)
pheatmap(mat=res_all,
         color             = inferno(10),
         border_color      = NA,
         show_colnames     = T,
         show_rownames     = T,
         drop_levels       = TRUE,
         scale="row",
         main="moSTB significant pathways differentiated in Trm-like cluster",
         fontsize=20)
dev.off()





















############################################################################################################			
rm(list=ls())
#The up regulated
myinf1 <- "/lorax/chenglab/yanding/CoLab/Randy_VISTBA/Data/CD4_Polyclonal/Seruat_data_object_finished_CD4.RDS"
myinf3 <- "/lorax/chenglab/cc59/PubDat/Dataset/GSEA/GSEA_C2_all_geneset.txt"



data <- readRDS(myinf1)
raw.data <- data

data <- data@data
data <- as.matrix(data)
info <- raw.data@meta.data
row.names(data) <- toupper(row.names(data))

pathway <- read.table(myinf3,sep="\t",quote=NULL)
com <- intersect(row.names(data), row.names(pathway))
data <- data[com,]
pathway <- pathway[com,]
raw.pathway <- pathway

#Up genes
myinf <- "/lorax/chenglab/yanding/CoLab/Randy_VISTBA/Data/CD4_Polyclonal/Internal_cluster_C2_enrichment/"
files <- list.files(myinf)
tag <- grep("up",files)
files <- files[tag]
for(k in 1 : length(files))
{
  tmpinf <- paste0(myinf,files[k])
  
  res <- read.table(tmpinf, sep="\t",quote=NULL)
  tag <- res$p.adjuSTB < 0.06
  res <- res[tag,]
  xx1 <- row.names(res)
  
  mydir <- "/lorax/chenglab/yanding/CoLab/Randy_VISTBA/Data/CD4_Polyclonal/Internal_cluster_comparison/"
  myfiles <- list.files(mydir)
  pathway <- raw.pathway[,xx1]
  
  
  tmpinf <- paste0(mydir, myfiles[k])
  
  gene <- read.table(tmpinf,sep="\t",quote=NULL)
  tag <- is.na(gene$Pvalue)
  gene <- gene[tag==0,]
  tag <- gene$Pvalue < 0.06
  gene <- gene[tag==1,]
  tag <- gene$Avg_logFC > 0
  gene <- gene[tag,]
  row.names(gene) <- toupper(row.names(gene))
  
  tmp_pathway <- pathway[intersect(row.names(gene),row.names(pathway)),]
  
  label <- STBrsplit(myfiles[k],"_")
  cluster <- label[[1]][2]
  label <- paste0(label[[1]][1],"_",label[[1]][2])
  
  
  tag <- info$res.0.9 == cluster
  tmp_info <- info[tag,]
  tmp_data <- data[intersect(row.names(gene),row.names(data)),row.names(tmp_info)]
  tmp_data <- apply(tmp_data, 1, function(arg) (arg-mean(arg))/sd(arg))
  tmp_data <- t(tmp_data)
  
  for(j in 1 : ncol(tmp_pathway))
  {
    cat("\r",k,"--->--->",j)
    target <- tmp_pathway[,j]
    names(target) <- row.names(tmp_pathway)
    target <- target[target==1]
    
    tmp_res <- tmp_data[intersect(row.names(tmp_data),names(target)),]
    
    tag1 <- grep("KO",colnames(tmp_res))
    tag2 <- grep("WT",colnames(tmp_res))
    
    tmp_res1 <- tmp_res[,tag1]
    tmp_res2 <- tmp_res[,tag2]
    KO <- apply(tmp_res1,1, mean)
    WT <- apply(tmp_res2,1,mean)
    
    tmp_res <- cbind(KO,WT)
    #tmp_res <- apply(tmp_res, 1, function(arg) scale(arg))
    #tmp_res <- t(tmp_res)
    
    tag <- tmp_res > 0.1
    tmp_res[tag] <- 0.1
    
    tag <- tmp_res < -0.1
    tmp_res[tag] <- (-1)*0.1
    
    myoutf <- paste0("/lorax/chenglab/yanding/CoLab/Randy_VISTBA/Figures/CD4_Polyclonal/Sp_internal_cluster_Pathway/",label,"_Up_",colnames(tmp_pathway)[j],".pdf")		
    pdf(myoutf)
    pheatmap(tmp_res)
    dev.off()
  }		
}		























############################################################################
############################################################################

current.cluster.ids <- c("0","1","2","3","4","5","6")
new.cluster.ids<- c("effector memory like PW_Adult_kit.integrated T cells", "activated PW_Adult_kit.integrated T cells","exhauSTBed PW_Adult_kit.integrated T cells","memory like PW_Adult_kit.integrated T cells to be defined","memory like PW_Adult_kit.integrated T cells", 
                    "resident memory PW_Adult_kit.integrated T cells", "CD4 monocytes/dendritic cells")
PW_Adult_kit.integrated@ident <- plyr::mapvalues(x = PW_Adult_kit.integrated@ident, from = current.cluster.ids, to = new.cluster.ids)
myoutf <- "/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/combine/RNA/vitiligo/STB1356/skin_tumor/PW_Adult_kit.integrated/PCx/res.0.3/output/res.1/TSNE1_defined.pdf"
pdf(myoutf,width=8,height=5)
TSNEPlot(object = PW_Adult_kit.integrated, pt.size = 0.7,label.size = 2)
dev.off()

PW_Adult_kit.integrated <- STBashIdent(object = PW_Adult_kit.integrated, save.name = "CellType")
for(i in 1:length(nrow(PW_Adult_kit.integrated@meta.data)))
{
  PW_Adult_kit.integrated@meta.data$CellType<-paSTBe(PW_Adult_kit.integrated@meta.data$orig.ident,PW_Adult_kit.integrated@meta.data$CellType,sep="_")
}

myoutf <- "/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/combine/RNA/vitiligo/STB1356/skin_tumor/PW_Adult_kit.integrated/PCx/res.0.3/output/res.1/TSNE2_defined.pdf"
pdf(myoutf,width=8,height=5)
TSNEPlot(object = PW_Adult_kit.integrated, pt.size = 0.7,label.size = 2,group.by="CellType")
dev.off()

PW_Adult_kit.integrated <- SetAllIdent(PW_Adult_kit.integrated, id = "CellType")
plot<-SubsetData(object=PW_Adult_kit.integrated,ident.use=c("tumor_effector memory like PW_Adult_kit.integrated T cells",
                                                "tumor_activated PW_Adult_kit.integrated T cells",
                                                "tumor_exhauSTBed PW_Adult_kit.integrated T cells",
                                                "skin_exhauSTBed PW_Adult_kit.integrated T cells",
                                                "tumor_memory like PW_Adult_kit.integrated T cells to be defined",
                                                "skin_memory like PW_Adult_kit.integrated T cells",
                                                "tumor_memory like PW_Adult_kit.integrated T cells",
                                                "skin_resident memory PW_Adult_kit.integrated T cells",
                                                "tumor_resident memory PW_Adult_kit.integrated T cells",
                                                "skin_CD4 monocytes/dendritic cells"
))
myoutf <- "/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/combine/RNA/vitiligo/STB1356/skin_tumor/PW_Adult_kit.integrated/PCx/res.0.3/output/res.1/TSNE_defined_by_tissue.pdf"
pdf(myoutf,width=8.3,height=5)
TSNEPlot(object = plot, pt.size = 0.7,label.size = 2.5)
dev.off()

myoutf <- "/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/combine/RNA/vitiligo/STB1356/skin_tumor/PW_Adult_kit.integrated/PCx/res.0.3/output/TSNE_defined_unlabel.pdf"
pdf(myoutf,width=10,height=6.5)
TSNEPlot(object = plot, pt.size = 0.5,label.size = 2.5)
dev.off()

myoutf <- "/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/combine/RNA/vitiligo/STB1356/skin_tumor/PW_Adult_kit.integrated/Seruat_PW_Adult_kit.integrated_undefined_finished.RDS"
saveRDS(PW_Adult_kit.integrated, file = myoutf)
myinf<-
  
  
  
  
  [6]Validation bias
myoutf <- "/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/combine/RNA/vitiligo/STB1356/skin_tumor/PW_Adult_kit.integrated/TSNE_plot_UMI.pdf"
pdf(myoutf,width=5,height=5)
FeaturePlot(skin_tumor, features.plot=c('nUMI'), pt.size=0.5)
dev.off()

myoutf <- "/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/combine/RNA/vitiligo/STB1356/skin_tumor/PW_Adult_kit.integrated/TSNE_plot_nGene.pdf"
pdf(myoutf,width=5,height=5)
FeaturePlot(skin_tumor, features.plot=c('nGene'), pt.size=0.5)
dev.off()

myoutf <- "/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/combine/RNA/vitiligo/STB1356/skin_tumor/PW_Adult_kit.integrated/TSNE_plot_Mitochondira.pdf"
pdf(myoutf,width=5,height=5)
FeaturePlot(skin_tumor, features.plot=c('percent.mito'), pt.size=0.5)
dev.off()

myoutf <- "/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/combine/RNA/vitiligo/STB1356/skin_tumor/PW_Adult_kit.integrated/TSNE_plot_undefined.pdf"
pdf(myoutf,width=5,height=5)
TSNEPlot(object = skin_tumor)
dev.off()

#Do unsupervised clustering
skin_tumor <- BuildclusterTree(
  skin_tumor,
  pcs.use = 1:12,
  do.reorder = F,
  reorder.numeric = F,
  do.plot=F)

myoutf <- "/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/combine/RNA/vitiligo/STB1356/skin_tumor/PW_Adult_kit.integrated/cluster_tree_validation.pdf"
pdf(myoutf,width=5,height=5)
PlotclusterTree(skin_tumor)
dev.off()

myoutf <- "/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/combine/RNA/vitiligo/STB1356/skin_tumor/PW_Adult_kit.integrated/TSNE_plot_merge.pdf"
pdf(myoutf,width=5,height=5)
TSNEPlot(object = skin_tumor,do.label=T,group.by="res.0.2")
dev.off()


myoutf <- "/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/combine/RNA/vitiligo/STB1356/skin_tumor/PW_Adult_kit.integrated/Seruat_data_object_finished.RDS"
saveRDS(skin_tumor, file = myoutf)

[8]Plot the diSTBribution of marker genes
PW_Adult_kit.integratedNaiveMarkers <- c("Ccr7", "Sell", "Tcf7", "Lef1", "Ccl3", "Ccl5", "Erdr1")
T_CD4_TRegMarkers <- c("Il2ra", "Entpd1", "Cd40", "Nt5e", "Itgae", "Tnfrsf18", "Foxp3")
T_CD4_memory <- c("Il2ra","Cd44", "Ptprc", "Sell", "Il7r", "Tcra", "Tcrb","Ccr9","Ccr7")
T_CD4_Th9 <- c("Il9")
T_CD4_Tfh <- c("TBKBP1", "HOPX", "TMEM106A", "MORN4", "CHSTB11", "TAP1", "IFITM3", "IFITM2", "GBP2", "DNAJC15", "GIMAP4", "SLAMF7", "IL12RB2")
T_CD4_Th2 <- c("Gata3", "Rap1b", "Trat1", "Wdfy2")
T_CD4_Th1 <- c("Ifngr1","Fasl", "Cxcr3", "Ccr5", "Il12rb1", "Il18r1")
T_CD4_Th17 <- c("Ptprc","Il6ra", "Klrb1c", "Ccr4", "Ccr6", "Il21r", "Il17a", "Il17f")
T_CD4_Th22 <- c("Ccr10", "Pdgfra","Ccr4", "Ccr6")

## PW_Adult_kit.integrated T cells
T_PW_Adult_kit.integrated_cytMarkers <- c("PW_Adult_kit.integrateda", "PW_Adult_kit.integratedb1", "Prf1", "Nkg7", "Gzma", "Gzmb", "Gzmh", "Gzmk", "Cxcr5", "Klrk1", "Tcra", "Tcrb")
T_PW_Adult_kit.integrated_exhMarkers <- c("Lag3", "Pdcd1", "Cblb", "Cxcl13", "Tim-3", "Ctla4", "Tigit", "Prdm1", "Icos", "Sell")
T_PW_Adult_kit.integrated_Mem <- c("Cd3e", "PW_Adult_kit.integrateda", "PW_Adult_kit.integratedb1", "Cd44", "Ptprc", "Tcra", "Tcrb")

# Other lymphoid cells
T_GammaDelta <- c("Cd2", "Cd3e", "Cd4", "PW_Adult_kit.integrateda", "PW_Adult_kit.integratedb1", "Tcrg-C1", "Tcrg-C2", "Tcrg-C3", "Tcrg-C4", "Tcrg-V3", "Tcrg-V4", "Tcrg-V5", "Tcrg-V6", "Tcrg-V7")
NKT_cell <- c("Cd1d1", "Cd3e", "Cd160", "Zbtb16")
ILC_cells <- c("Cd4", "Il7r", "Ccr6", "Rorc", "Ptprc", "Klrg1", "Ly6a", "Il1rl1", "Klrd1", "Itgae", "Klrb1c")
NK_cells <- c("Itga2", "Klrd1", "Cd96", "Cd161", "Ncr1")

#NKT_cell
myoutf <-"/lorax/chenglab/yanding/CoLab/Randy_VISTBA/Figures/2w1s_WT_vs_KO_scRNA/Cell_identify/NKT_cell_Naive.pdf"
pdf(myoutf,width=10,height=10)
FeaturePlot(object = data, features.plot = NKT_cell,cols.use = c("grey", "red"),
            reduction.use = "tsne")
dev.off()

#NK_cells
myoutf <-"/lorax/chenglab/yanding/CoLab/Randy_VISTBA/Figures/2w1s_WT_vs_KO_scRNA/Cell_identify/NK_cell_Naive.pdf"
pdf(myoutf,width=10,height=10)
FeaturePlot(object = data, features.plot = intersect(row.names(data@data),NK_cells),cols.use = c("grey", "red"),
            reduction.use = "tsne")
dev.off()

#PW_Adult_kit.integrated T cells
myoutf <-"/lorax/chenglab/yanding/CoLab/Randy_VISTBA/Figures/2w1s_WT_vs_KO_scRNA/Cell_identify/T_PW_Adult_kit.integrated_Mem_cell.pdf"
pdf(myoutf,width=10,height=10)
FeaturePlot(object = data, features.plot = intersect(row.names(data@data),T_PW_Adult_kit.integrated_Mem),cols.use = c("grey", "red"),
            reduction.use = "tsne")
dev.off()


#Begin to annotate T cell subype

#Naive T
myoutf <-"/lorax/chenglab/yanding/CoLab/Randy_VISTBA/Figures/2w1s_WT_vs_KO_scRNA/Cell_identify/T_cell_Naive.pdf"
pdf(myoutf,width=10,height=10)
FeaturePlot(object = data, features.plot = PW_Adult_kit.integratedNaiveMarkers,cols.use = c("grey", "red"),
            reduction.use = "tsne")
dev.off()

#CD4Treg
myoutf <-"/lorax/chenglab/yanding/CoLab/Randy_VISTBA/Figures/2w1s_WT_vs_KO_scRNA/Cell_identify/CD4_Treg_cell.pdf"
pdf(myoutf,width=10,height=10)
FeaturePlot(object = data, features.plot = intersect(row.names(data@data),T_CD4_TRegMarkers),cols.use = c("grey", "red"),
            reduction.use = "tsne")
dev.off()

#CD4 memory
myoutf <-"/lorax/chenglab/yanding/CoLab/Randy_VISTBA/Figures/2w1s_WT_vs_KO_scRNA/Cell_identify/CD4_Mem_cell.pdf"
pdf(myoutf,width=10,height=10)
FeaturePlot(object = data, features.plot = intersect(row.names(data@data),T_CD4_memory),cols.use = c("grey", "red"),
            reduction.use = "tsne")
dev.off()

#T_CD4_Th2
myoutf <-"/lorax/chenglab/yanding/CoLab/Randy_VISTBA/Figures/2w1s_WT_vs_KO_scRNA/Cell_identify/T_CD4_Th2_cell.pdf"
pdf(myoutf,width=10,height=10)
FeaturePlot(object = data, features.plot = intersect(row.names(data@data),T_CD4_Th2),cols.use = c("grey", "red"),
            reduction.use = "tsne")
dev.off()

#T_CD4_Th1

myoutf <-"/lorax/chenglab/yanding/CoLab/Randy_VISTBA/Figures/2w1s_WT_vs_KO_scRNA/Cell_identify/T_CD4_Th1_cell.pdf"
pdf(myoutf,width=10,height=10)
FeaturePlot(object = data, features.plot = intersect(row.names(data@data),T_CD4_Th1),cols.use = c("grey", "red"),
            reduction.use = "tsne")
dev.off()

#T_CD4_Th17
myoutf <-"/lorax/chenglab/yanding/CoLab/Randy_VISTBA/Figures/2w1s_WT_vs_KO_scRNA/Cell_identify/T_CD4_Th17_cell.pdf"
pdf(myoutf,width=10,height=10)
FeaturePlot(object = data, features.plot = intersect(row.names(data@data),T_CD4_Th17),cols.use = c("grey", "red"),
            reduction.use = "tsne")
dev.off()

[9]T cell STBate difference

#T_Cell_Anergy
T_Cell_Anergy <- c("Cblb","Cbl","NFAT","Lag3","Ctla4","NT5e","Izumo1r","Dgka","Rap1","Itch","Rnf128","Dtx1")
myoutf <-"/lorax/chenglab/yanding/CoLab/Randy_VISTBA/Figures/2w1s_WT_vs_KO_scRNA/Cell_identify/T_Cell_Anergy_cell.pdf"
pdf(myoutf,width=10,height=10)
FeaturePlot(object = data, features.plot = intersect(row.names(data@data),T_Cell_Anergy),cols.use = c("grey", "red"),
            reduction.use = "tsne")
dev.off()

#Type_I_IFN
#change the data gene name

tmp_data <- data
row.names(tmp_data@data) <- toupper(tmp_data@data)

Type_I_IFN_1 <- c("IRF1","IFIH1","IFITM3","DDX58","IFI44L","IFI6","IFITM2","NAMPT","OASL","RTP4")
Type_I_IFN_2 <- c("TREX1","ADAR","FAM46C","LY6E","MCOLN2","APOBEC3G","IL15",
                  "ISG15","MX1","TLR3")
myoutf <-"/lorax/chenglab/yanding/CoLab/Randy_VISTBA/Figures/2w1s_WT_vs_KO_scRNA/Cell_identify/Type_I_IFN_1_cell.pdf"
pdf(myoutf,width=10,height=10)

FeaturePlot(object = tmp_data, features.plot = intersect(row.names(tmp_data@data),Type_I_IFN_1),cols.use = c("grey", "red"),
            reduction.use = "tsne")
dev.off()

[10]
#find markers
data.markers <- FindAllMarkers(object = data, only.pos = TRUE, min.pct = 0.25,
                               thresh.use = 0.25)

myoutf <- "/lorax/chenglab/yanding/CoLab/Randy_VISTBA/Data/Single_Cell_RNA_CD4_CD44_Lo_Poly/Marker_for_each_cluster.xls"
write.table(data.markers,myoutf,sep="\t",quote=F)


[11]Difference expression between WT and KO
myinf <- "/lorax/chenglab/yanding/CoLab/Randy_VISTBA/Data/Single_Cell_RNA_CD4_CD44_Lo_Poly/Seruat_data_object_finished.RDS"
data <- readRDS(myinf)
data@meta.data[,"label"] <- rep("WT",nrow(data@meta.data))
tag <- grep("KO",data@meta.data[,"orig.ident"])
data@meta.data[tag,"label"] <- "KO"

data <- SetAllIdent(data, id = "label")

res <- FindMarkers(data, ident.1 = "KO", ident.2 = "WT",logfc.threshold = 0.01,min.pct = 0.01)
tag <- res$p_val_adj <0.05
res <- res[tag,]

myoutf <- "/lorax/chenglab/yanding/CoLab/Randy_VISTBA/Data/Single_Cell_RNA_CD4_CD44_Lo_Poly/Wilcox_KO_vs_WT.xls"
write.table(res,myoutf,sep="\t",quote=F)



[12]Go enrichment analysis
rm(list=ls())
myinf1 <- "/lorax/chenglab/yanding/CoLab/Randy_VISTBA/Data/Single_Cell_RNA_CD4_CD44_Lo_Poly/Wilcox_KO_vs_WT.xls"
myinf2 <- "/lorax/chenglab/yanding/Pub_Dat/homo_Go.gaf"
myinf3 <- "/lorax/chenglab/yanding/Pub_Dat/Goterms.txt"
myinf4 <- "/lorax/chenglab/yanding/Pub_Dat/Go_domain.txt"

myoutf1 <- "/lorax/chenglab/yanding/CoLab/Randy_VISTBA/Data/Single_Cell_RNA_CD4_CD44_Lo_Poly/CD4_CD44_Lo_Poly_K__vs__W__up_gene_enrichment.xls"
myoutf2 <- "/lorax/chenglab/yanding/CoLab/Randy_VISTBA/Data/Single_Cell_RNA_CD4_CD44_Lo_Poly/CD4_CD44_Lo_Poly_K__vs__W__dn_gene_enrichment.xls"

res <- read.table(myinf1,sep="\t",quote=NULL,header=TRUE)
res[,"HGNC.symbol"] <- toupper(row.names(res))


info <- read.delim(myinf2,sep="\t",skip=34,header= FALSE,STBringsAsFactors = FALSE)
terms <- read.table(myinf3,sep="\t",STBringsAsFactors = FALSE)
go <- read.table(myinf4,sep="\t",quote=NULL,STBringsAsFactors=F,header=T)


colnames(info)[3] <- "HGNC.symbol"
colnames(info)[5] <- "GO.Term.Accession"

tag_up <- res$avg_logFC > 0
tag_dn <- res$avg_logFC < 0

dat1 <-res[tag_up==1,]#142
dat1[,"gene"] <- dat1[,"HGNC.symbol"]
row.names(dat1) <- seq(1,nrow(dat1),1)

dat2 <- res[tag_dn==1,]#127
dat2[,"gene"] <- dat2[,"HGNC.symbol"]
row.names(dat2) <- seq(1,nrow(dat2),1)

useinfo <-c("GO.Term.Accession","HGNC.symbol")
info <- info[,useinfo]
terms[,"GO.Term.Accession"]<- row.names(terms)
row.names(terms) <- seq(1,nrow(terms),1)

info <- merge(info,terms,by.x="GO.Term.Accession",by.y="GO.Term.Accession")
colnames(info)[3] <- "GO.Term.Name"
allgene <- unique(info[,"HGNC.symbol"]) #19050
up_info <- merge(dat1,info,by.x="gene",by.y="HGNC.symbol")

up_domain <- unique(up_info[,"GO.Term.Name"]) #9670
up_domain <- up_domain[up_domain!=""]
res <- matrix(nrow=length(up_domain),ncol=5)
colnames(res) <- c("enrichment","p.value","odds.ratio","shared_gene","all_gene")

for (i in 1 : length(up_domain))
{
  cat("\r",i)
  tag <- up_info[,"GO.Term.Name"] == up_domain[i]
  tmp <- up_info
  useinfo<- c("gene","GO.Term.Name")
  tmp <- tmp[,useinfo]
  tmp[tag==0,"GO.Term.Name"] <- "Others"
  tmp[,"GO.Term.Name"] <- factor(tmp[,"GO.Term.Name"])
  tmp_1<- tmp[tag==1,]
  tmp_1 <- unique(tmp_1)
  i_gene <- tmp_1[,"gene"]
  xx1 <-length(i_gene)
  tmp_2 <- tmp[tag==0,]
  tmp_2 <- unique(tmp_2)
  tag_r <- is.na(pmatch(tmp_2[,"gene"],i_gene))
  tmp_2 <- tmp_2[tag_r==1,]
  tmp <- rbind(tmp_1,tmp_2)
  xx2 <- nrow(tmp)
  white <- info[which(info[,"GO.Term.Name"] %in% up_domain[i]),"HGNC.symbol"]
  white <- unique(white)
  black <- is.na(pmatch(allgene,white))
  black <- unique(allgene[black])
  white <- length(white)
  black <- length(black)
  ER = (xx1/xx2)/(white/(white+black))
  myp = sum(dhyper(xx1:xx2, white, black, xx2))
  res[i,"p.value"] <- myp
  res[i,"odds.ratio"] <- ER
  res[i,"shared_gene"] <- xx1
  res[i,"all_gene"] <- white
  res[i,"enrichment"] <-up_domain[i]
  
  
}
cutoff <- 10
tag_cut <- as.numeric(res[,"shared_gene"]) > 5
res <- res[tag_cut==1,]
p.value <- res[,"p.value"]
p.value <- as.numeric(p.value)
p.adjuSTB <- p.adjuSTB(p.value,method = "BH")
res <- as.data.frame(res)
res[,"p.adjuSTB"] <- p.adjuSTB
res <- as.matrix(res)
rank <- order(as.numeric(res[,"p.value"]),decreasing=FALSE)
res <- res[rank,]
write.table(res,myoutf1,sep="\t", quote=F)


[]Deletion
down_info <- merge(dat2,info,by.x="gene",by.y="HGNC.symbol")
down_domain <- unique(down_info[,"GO.Term.Name"])
down_domain <- down_domain[down_domain!=""] #9328
res <- matrix(nrow=length(down_domain),ncol=5)
colnames(res) <- c("enrichment","p.value","odds.ratio","shared_gene","all_gene")
value1 <- rep(length(down_domain),0)
value2 <- rep(length(down_domain),0)
for (i in 1 : length(down_domain))
{
  cat("\r",i)
  tag <- down_info[,"GO.Term.Name"] == down_domain[i]
  tmp <- down_info
  useinfo<- c("gene","GO.Term.Name")
  tmp <- tmp[,useinfo]
  tmp[tag==0,"GO.Term.Name"] <- "Others"
  tmp[,"GO.Term.Name"] <- factor(tmp[,"GO.Term.Name"])
  tmp_1<- tmp[tag==1,]
  tmp_1 <- unique(tmp_1)
  i_gene <- tmp_1[,"gene"]
  xx1 <-length(i_gene)
  tmp_2 <- tmp[tag==0,]
  tmp_2 <- unique(tmp_2)
  tag_r <- is.na(pmatch(tmp_2[,"gene"],i_gene))
  tmp_2 <- tmp_2[tag_r==1,]
  tmp <- rbind(tmp_1,tmp_2)
  xx2 <- nrow(tmp)
  white <- info[which(info[,"GO.Term.Name"] %in% down_domain[i]),"HGNC.symbol"]
  white <- unique(white)
  black <- is.na(pmatch(allgene,white))
  black <- unique(allgene[black])
  white <- length(white)
  black <- length(black)
  ER = (xx1/xx2)/(white/(white+black))
  myp = sum(dhyper(xx1:xx2, white, black, xx2))
  res[i,"p.value"] <- myp
  res[i,"odds.ratio"] <- ER
  res[i,"shared_gene"] <- xx1
  res[i,"all_gene"] <- white
  res[i,"enrichment"] <-down_domain[i]
  value1[i] <- xx1
  value2[i] <- white
}

cutoff <- 10
tag_cut <- as.numeric(res[,"shared_gene"]) > 5
res <- res[tag_cut==1,]
p.value <- res[,"p.value"]
p.value <- as.numeric(p.value)
p.adjuSTB <- p.adjuSTB(p.value,method = "BH")
res <- as.data.frame(res)
res[,"p.adjuSTB"] <- p.adjuSTB
res <- as.matrix(res)
rank <- order(as.numeric(res[,"p.value"]),decreasing=FALSE)
res <- res[rank,]

write.table(res,myoutf2,sep="\t", quote=F)

[13]Pathway analysis
rm(list=ls())
myinf1 <-"/lorax/chenglab/yanding/CoLab/Randy_VISTBA/Data/Single_Cell_RNA_CD4_CD44_Lo_Poly/Wilcox_KO_vs_WT.xls"
myinf2 <- "/lorax/chenglab/cc59/PubDat/Dataset/GSEA/GSEA_C2_all_geneset.txt"

data <- read.table(myinf1,sep="\t",quote=NULL)
pathway.data <- read.delim(myinf2)
dn_enrichment_results <- matrix(0, ncol(pathway.data), 2)
row.names(dn_enrichment_results) <- colnames(pathway.data)
colnames(dn_enrichment_results) <- c("Pvalue","Odd.ratio")
up_enrichment_results <- dn_enrichment_results

tar <- toupper(row.names(data))	
tag <- duplicated(tar)
data <- data[tag==0,]
row.names(data) <- toupper(row.names(data))

int.gene <- intersect(row.names(data), row.names(pathway.data))
myres <- data[int.gene, ]
pathway.data.match <- pathway.data[int.gene, ]

################################################################################################
################################################################################################

tag1 <- data$avg_logFC > 0
tag2 <- data$avg_logFC < 0


up.gene <- row.names(data)[tag1]
dn.gene <- row.names(data)[tag2]

de.gene <- c(up.gene, dn.gene)
tag <- which(row.names(pathway.data) %in% de.gene)
non.de.gene <- row.names(pathway.data)[-tag]

for(j in 1:ncol(pathway.data.match)) {
  
  
  target.index <- pathway.data.match[, j] == 1
  target.gene <- row.names(pathway.data.match)[target.index]
  non.target.gene <- row.names(pathway.data.match)[!target.index]
  
  ## Enrichment of up-regulated genes
  
  n1 <- length(intersect(up.gene, target.gene))
  n2 <- length(intersect(up.gene, non.target.gene))
  n3 <- length(intersect(non.de.gene, target.gene))
  n4 <- length(intersect(non.de.gene, non.target.gene))
  
  confusion.matrix <- matrix(c(n1, n2, n3, n4), 2, 2, byrow = TRUE)
  
  fisher.obj <- fisher.teSTB(confusion.matrix, alternative = "greater")
  
  up_enrichment_results[j, "Pvalue"] <- fisher.obj[["p.value"]]
  up_enrichment_results[j, "Odd.ratio"]	<- as.numeric(fisher.obj[3])
  
  ## Enrichment of dn-regulated genes
  
  n1 <- length(intersect(dn.gene, target.gene))
  n2 <- length(intersect(dn.gene, non.target.gene))
  n3 <- length(intersect(non.de.gene, target.gene))
  n4 <- length(intersect(non.de.gene, non.target.gene))
  
  confusion.matrix <- matrix(c(n1, n2, n3, n4), 2, 2, byrow = TRUE)
  
  fisher.obj <- fisher.teSTB(confusion.matrix, alternative = "greater")
  
  dn_enrichment_results[j, "Pvalue"] <- fisher.obj[["p.value"]]
  dn_enrichment_results[j, "Odd.ratio"]	<- as.numeric(fisher.obj[3])
}



myoutf1 = "/lorax/chenglab/yanding/CoLab/Randy_VISTBA/Results/up_enrichments_results_C2.txt"
myoutf2 = "/lorax/chenglab/yanding/CoLab/Randy_VISTBA/Results/dn_enrichments_results_C2.txt"

up_enrichment_results <- up_enrichment_results[order(up_enrichment_results[,"Act_vs_WT"],decreasing=F),]
up_enrichment_results <- as.data.frame(up_enrichment_results)
up_enrichment_results_FDR <- apply(up_enrichment_results, 2, function(x) p.adjuSTB(x, method="BH"))

colnames(up_enrichment_results) <- paste0(colnames(up_enrichment_results),"_Pval")
colnames(up_enrichment_results_FDR) <- paste0(colnames(up_enrichment_results),"_FDR")

dn_enrichment_results <- dn_enrichment_results[order(dn_enrichment_results[,"Act_vs_WT"],decreasing=F),]
dn_enrichment_results <- as.data.frame(dn_enrichment_results)
dn_enrichment_results_FDR <- apply(dn_enrichment_results, 2, function(x) p.adjuSTB(x, method="BH"))

colnames(dn_enrichment_results) <- paste0(colnames(dn_enrichment_results),"_Pval")
colnames(dn_enrichment_results_FDR) <- paste0(colnames(dn_enrichment_results),"_FDR")

up_enrichment_results <- cbind(up_enrichment_results, up_enrichment_results_FDR)
dn_enrichment_results <- cbind(dn_enrichment_results, dn_enrichment_results_FDR)

write.table(up_enrichment_results, myoutf1,sep="\t",quote=F)
write.table(dn_enrichment_results, myoutf2,sep="\t",quote=F)





#compare with WT and KO
tag <- grep("KO",data@meta.data[,"orig.ident"])
data@meta.data[tag,"label"] <- "KO"

TSNEPlot(object = data,group.by="label")


tag <- grep("KO",data@meta.data[,"orig.ident"])
sam1 <- row.names(data@meta.data)[tag]

tag <- grep("WT",data@meta.data[,"orig.ident"])
sam2 <- row.names(data@meta.data)[tag]

TSNEPlot(object = data,cells.use=sam1)


info <- data@meta.data

tag <- grep("KO",info$label)
info1 <- info[tag,]

tag <- grep("WT",info$label)
info2 <- info[tag,]

summary(as.factor(info1$res.0.3))
0    1    2    3    4    5
1627 7076 6316 1709  728  425

summary(as.factor(info2$res.0.3))
0    1    2    3    4    5
9218 2081 1545  444  551  462

###############################################
#T cell STBate annotation
###############################################
[1]Validation of measure of CD4T cells

[2]Annotation different T cell STBatus
myinf1 <- "/lorax/chenglab/yanding/CoLab/Randy_VISTBA/Data/T_cell_STBate_marker/


###################trash teSTB#######try to validate there is no matched skin resident population in GE and TCR data######
a=0
z=0
for(i in 1:length(resident_skin))
{
   z<-grep(resident_skin[i],barcode_skin)==0
   a<-a+z
}














































