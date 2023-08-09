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
library(GSVA)
library(clustree)
library(rlang)
rm(list=ls()) 

myinf1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/raw/Nan_PW_Normal"
PW_C.data <- Read10X(data.dir=myinf1)
PW_C <- CreateSeuratObject(counts = PW_C.data, project = "PW_C",min.cells=3,min.features=200)
PW_C<-RenameCells(object=PW_C,add.cell.id=unique(PW_C@meta.data$orig.ident))
PW_C<-RenameCells(object=PW_C,new.names=gsub('-1','',Cells(PW_C)))

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

xx<- merge(PW_C, y = c(PW_Q,PW_S,PW_U,PW_W,PW_Adu1,PW_Adu2), project = "original")
PW_Adult_only.list<-list(PW_C,PW_Q,PW_S,PW_U,PW_W,PW_Adu1,PW_Adu2)
PW_Adult_only.list <- lapply(X = PW_Adult_only.list, FUN = function(x) {
    x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")
    x <- subset(x, subset = nCount_RNA<40000 & nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})

 features <- SelectIntegrationFeatures(object.list = PW_Adult_only.list,nfeatures = 3000)
 PW_Adult_only.list <- lapply(X = PW_Adult_only.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
 })

 immune.anchors <- FindIntegrationAnchors(object.list = PW_Adult_only.list, anchor.features = features, reduction = "rpca")
 PW_Adult_only.integrated <- IntegrateData(anchorset = immune.anchors)
 DefaultAssay(PW_Adult_only.integrated) <- "integrated"

 PW_Adult_only.integrated <- ScaleData(PW_Adult_only.integrated, verbose = FALSE)
  PW_Adult_only.integrated <- RunPCA(PW_Adult_only.integrated, npcs = 30, verbose = FALSE)


 x<-c(1:13)

 PW_Adult_only.integrated<-FindNeighbors(PW_Adult_only.integrated, dims = x)
 PW_Adult_only.integrated <- FindClusters(PW_Adult_only.integrated,resolution = seq(0.1,2.0,0.1))
 PW_Adult_only.integrated <- RunUMAP(PW_Adult_only.integrated, dims= x)



#######Subset only MNPs for downstream analysis###########################

 tar<-"integrated_snn_res.0.7"
 Idents(object = PW_Adult_only.integrated) <- tar
 subcluster<-c('0','2','3','5','8','9','11','12','13','14','17','20')
 tag<-PW_Adult_only.integrated@meta.data[,tar]%in%subcluster
 tag1<-GetAssayData(PW_Adult_only.integrated,assay='RNA', slot='data')['CD3E',]>0.05
 cell_ID<-row.names(PW_Adult_only.integrated@meta.data)[tag&!tag1]
 
 Myeloid_PW_Adult_only.integrated<-subset(x=xx,cells=cell_ID,nCount_RNA>3000)

#PW_Adu1 PW_Adu2    PW_C    PW_Q    PW_S    PW_U    PW_W 
#   2858    2438    5196    2566    3575    1176    4083 

 Myeloid_PW_Adult_only.integrated.list <- SplitObject(Myeloid_PW_Adult_only.integrated, split.by = "orig.ident")

 Myeloid_PW_Adult_only.integrated.list <- lapply(X = Myeloid_PW_Adult_only.integrated.list, FUN = function(x) {
    x <- NormalizeData(x,assay='RNA')
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000,assay='RNA')
 })

 features <- SelectIntegrationFeatures(object.list = Myeloid_PW_Adult_only.integrated.list,nfeatures = 3000)
 Myeloid_PW_Adult_only.integrated.list <- lapply(X = Myeloid_PW_Adult_only.integrated.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = T)
    x <- RunPCA(x, features = features, verbose = T)
 })

 immune.anchors <- FindIntegrationAnchors(object.list = Myeloid_PW_Adult_only.integrated.list, anchor.features = features, reduction = "rpca")
 Myeloid_PW_Adult_only.integrated <- IntegrateData(anchorset = immune.anchors)
 
 DefaultAssay(Myeloid_PW_Adult_only.integrated) <- "integrated"

 Myeloid_PW_Adult_only.integrated <- ScaleData(Myeloid_PW_Adult_only.integrated, verbose = T)
  Myeloid_PW_Adult_only.integrated <- RunPCA(Myeloid_PW_Adult_only.integrated, npcs = 30, verbose = T)
 x<-c(1:15)

 Myeloid_PW_Adult_only.integrated<-FindNeighbors(Myeloid_PW_Adult_only.integrated, dims = x)
 Myeloid_PW_Adult_only.integrated <- FindClusters(Myeloid_PW_Adult_only.integrated,resolution = seq(0.1,1.0,0.1))
 Myeloid_PW_Adult_only.integrated <- RunUMAP(Myeloid_PW_Adult_only.integrated, dims= x)


 myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_only_Human_nCountRNA_modified/Seruat_Myeloid_PW_Adult_only.integrated__finished.RDS"
 saveRDS(Myeloid_PW_Adult_only.integrated, file = myoutf)

#############ClusterTree###################################

myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_only_Human_nCountRNA_modified/Seruat_Myeloid_PW_Adult_only.integrated__finished.RDS"
Myeloid_PW_Adult_only.integrated<-readRDS(myinf)
 myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_only_Human_nCountRNA_modified/PCx/clustertree/cluster_tree_res_0.1-1.pdf"
  pdf(myoutf,width=8,height=12)
  clustree(Myeloid_PW_Adult_only.integrated, prefix = "integrated_snn_res.")
 dev.off()

##############The best resolution is 0.6 and recolor the clusters####################

Idents(object = Myeloid_PW_Adult_only.integrated) <- 'integrated_snn_res.0.6'
levels(Myeloid_PW_Adult_only.integrated@meta.data$integrated_snn_res.0.6)<-seq(0,15,1)
  myoutf1 <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_only_Human_nCountRNA_modified/PCx/res.0.6/UMAP_recolored.pdf"
  pdf(myoutf1,7,5)
  print(DimPlot(Myeloid_PW_Adult_only.integrated, reduction = "umap",cols=c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
            "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
            "#8dd3c7", "#bebada", "#fb8072", "#80b1d3", "#fdb462",
            "#b3de69")))
  dev.off()

############GSVA score for the converting macrophage gene signature################

load("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/PW_Mu_Allen/GSE189031_seurat.combined.Rdata")
table(sce.combined$label)
Idents(sce.combined)='label'
tmp.markers <- FindAllMarkers(object = sce.combined, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.2,assay='integrated')

B6_LCM_markers<-FindMarkers(object = sce.combined, only.pos = TRUE, ident.1='C57BL/6-like LCM',ident.2=c('B Cells','Converting SCM','DC-like SCM','Monocyte-like SCM','Naive Intermediate','Proliferating'),min.pct = 0.2, logfc.threshold = 0.2,assay='integrated')
BalbC_LCM_markers<-FindMarkers(object = sce.combined, only.pos = TRUE, ident.1='BALB/c-like LCM',ident.2=c('B Cells','Converting SCM','DC-like SCM','Monocyte-like SCM','Naive Intermediate','Proliferating'),min.pct = 0.2, logfc.threshold = 0.2,assay='integrated')
B6_LCM_sig<-row.names(subset(B6_LCM_markers,p_val_adj<0.05&avg_log2FC>0.5))
BalbC_LCM_sig<-row.names(subset(BalbC_LCM_markers,p_val_adj<0.05&avg_log2FC>0.5))

sig<-matrix(NA,300,length(unique(tmp.markers$cluster)))
colnames(sig)<-unique(tmp.markers$cluster)
for (i in 1:ncol(sig))
{
    tmp<-tmp.markers[tmp.markers$cluster==colnames(sig)[i],]
    tmp<-subset(tmp,p_val_adj<0.05&avg_log2FC>0.5)
    sig[1:nrow(tmp),i]<-tmp$gene
}

sig[1:length(B6_LCM_sig),'C57BL/6-like LCM']<-B6_LCM_sig
sig[1:length(BalbC_LCM_sig),'BALB/c-like LCM']<-BalbC_LCM_sig


Mouse_to_human <- function(x){

 human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
 mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

   genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T,verbose=F)

   humanx <- unique(genesV2[, 2])

   return(genesV2)
}

CM_sig<-apply(sig,2,Mouse_to_human)
CM_sig_GSVA<-lapply(CM_sig,function(x) x[,'HGNC.symbol'])

exp<-as.matrix(GetAssayData(Myeloid_PW_Adult_only.integrated,assay='integrated',slot='data'))
gsva.es_integrated <- gsva(exp, CM_sig_GSVA, verbose=FALSE)
Myeloid_PW_Adult_only.integrated[['CM_GSVA']]<-CreateAssayObject(gsva.es)
DefaultAssay(Myeloid_PW_Adult_only.integrated)='CM_GSVA'
Idents(Myeloid_PW_Adult_only.integrated)<-"integrated_snn_res.0.6"
Myeloid_PW_Adult_only.integrated<-ScaleData(Myeloid_PW_Adult_only.integrated)

myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_only_Human_nCountRNA_modified/PCx/res.0.6/GSVA/CM/CM_score_violin.pdf"
  pdf(myoutf,width=12,height=5)
 VlnPlot(object = Myeloid_PW_Adult_only.integrated,assay='CM_GSVA',slot='scale.data',features=row.names(Myeloid_PW_Adult_only.integrated[['CM_GSVA']]),ncol=2,pt.size=0.3)
 dev.off()


p1 <- (FeaturePlot(Myeloid_PW_Adult_only.integrated, features = c("GATA6-dependency"),slot='scale.data') & 
  scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red')) +
  ggtitle('GATA6-dependency activity')
p2 <- (FeaturePlot(Myeloid_PW_Adult_only.integrated, features = c("CM"),slot='scale.data',pt.size=0.5,min.cutoff=2, max.cutoff=2) & 
  scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red')) +
  ggtitle('Converting Mac Score')

myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_only_Human_nCountRNA_modified/PCx/res.0.6/GSVA/CM/converting_mac_score.pdf"
pdf(myoutf,width=20,height=10)
print(cowplot::plot_grid(p1,p2,ncol=2))
dev.off()


table(Myeloid_PW_Adult_only.integrated$integrated_snn_res.0.6)
CM_score<-matrix(NA,as.numeric(table(Myeloid_PW_Adult_only.integrated$integrated_snn_res.0.6)[1]),length(unique(Myeloid_PW_Adult_only.integrated$integrated_snn_res.0.6)))
colnames(CM_score)<-sort(as.character(unique(Myeloid_PW_Adult_only.integrated$integrated_snn_res.0.6)))

for (i in 1:ncol(CM_score))
{
 cell_ID<-row.names(Myeloid_PW_Adult_only.integrated@meta.data[Myeloid_PW_Adult_only.integrated$integrated_snn_res.0.6==colnames(CM_score)[i],])
 CM_score[1:length(cell_ID),colnames(CM_score)[i]]<-GetAssayData(Myeloid_PW_Adult_only.integrated,slot='scale.data')['CM',cell_ID]
}

myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_only_Human_nCountRNA_modified/PCx/res.0.6/GSVA/CM/CM_score_by_cluster.xls"
write.table(CM_score,myoutf,sep="\t",quote=F)




#############Feature plots######################
DefaultAssay(Myeloid_PW_Adult_only.integrated) <- "RNA"
paper_1<-c('GATA6','SELP','ITGAM','CD14')
paper_2<-c('CCR2','MRC1','LYVE1','TIMD4')
p<-FeaturePlot(object = Myeloid_PW_Adult_only.integrated, features = paper_1,cols= c('#deebf7', "#bd0026"),
               reduction = "umap",slot = "data",max.cutoff=4,pt.size=0.5,combine=F,order = TRUE)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + theme(plot.title = element_text(size=50),
                          axis.title = element_text(family="Helvetica",face = "bold", color = "black",size=40), 
                           axis.text=element_text(family="Helvetica",face = "bold", color = "black",size=30),
                           axis.line = element_line(colour = 'black', size = 3),
                           axis.ticks = element_line(colour = "black", size = 3),
                           axis.ticks.length=unit(1, "cm"))
}
 
myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_only_Human_nCountRNA_modified/PCx/paper_marker_Fig1B.pdf"
pdf(myoutf,width=22,height=20)
print(cowplot::plot_grid(plotlist = p,ncol=2))
dev.off()

p<-FeaturePlot(object = Myeloid_PW_Adult_only.integrated, features = paper_2,cols= c('#deebf7', "#bd0026"),
               reduction = "umap",slot = "data",max.cutoff=4,pt.size=0.5,combine=F,order = TRUE)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + theme(plot.title = element_text(size=50),
                          axis.title = element_text(family="Helvetica",face = "bold", color = "black",size=40), 
                           axis.text=element_text(family="Helvetica",face = "bold", color = "black",size=30),
                           axis.line = element_line(colour = 'black', size = 3),
                           axis.ticks = element_line(colour = "black", size = 3),
                           axis.ticks.length=unit(1, "cm"))
}
 
myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_only_Human_nCountRNA_modified/PCx/paper_marker_Fig1F.pdf"
pdf(myoutf,width=22,height=20)
print(cowplot::plot_grid(plotlist = p,ncol=2))
dev.off()


#############Cluster contribution by each different patient########## 
myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_only_Human_nCountRNA_modified/Seruat_Myeloid_PW_Adult_only.integrated__finished.RDS"
Myeloid_PW_Adult_only.integrated<-readRDS(myinf)

   Idents(Myeloid_PW_Adult_only.integrated)<-"integrated_snn_res.0.6"
 
 xx<-table(Myeloid_PW_Adult_only.integrated@meta.data$integrated_snn_res.0.6,Myeloid_PW_Adult_only.integrated@meta.data$orig.ident)
 yy<-t(apply(xx,1,function(x){x/sum(x)})   )
 myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_only_Human_nCountRNA_modified/PCx/res.0.6/cluster_patient_distribution/patient_distribution_to_each_cluster.xls"
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
 data$patient<-factor(data$patient,levels=unique(data$patient))

  col<-c('#ffffcc','#d9f0a3','#addd8e','#78c679','#41ab5d','#238443','#005a32')

 myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_only_Human_nCountRNA_modified/PCx/res.0.6/cluster_patient_distribution/stackplot_patient_contribution_each_cluster_proportion.pdf"
 pdf(myoutf,width=8,height=8)
 ggplot() + geom_bar(aes(y = value, x = cluster, fill = patient), data = data,
                    stat="identity")+theme(axis.text.x = element_text(angle = 60, hjust = 1,size=10),
                                           axis.title.x = element_blank())

 dev.off()



############Plot cluster proportions of total MNP##################################
tar<-unique(Myeloid_PW_Adult_only.integrated@meta.data$orig.ident)
proportion<-table(Myeloid_PW_Adult_only.integrated@meta.data$orig.ident,Myeloid_PW_Adult_only.integrated@meta.data$integrated_snn_res.0.6)
proportion<-as.matrix(proportion)
proportion<-as.data.frame(t(apply(proportion,1,function(x){x/sum(x)})))
proportion[,'GATA6']<-rep(0,nrow(proportion))

exp<-t(as.matrix(GetAssayData(Myeloid_PW_Adult_only.integrated,assay='RNA',slot='data')))

for (i in 1:length(tar))
{
  tag<-Myeloid_PW_Adult_only.integrated@meta.data$orig.ident==tar[i]
  cells<-row.names(Myeloid_PW_Adult_only.integrated@meta.data)[tag]
  exp_use<-exp[cells,]
  proportion[tar[i],'GATA6']<-sum(exp_use[,'GATA6']>0)/nrow(exp_use)
}


myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_only_Human_nCountRNA_modified/PCx/res.0.6/cluster_proportion_MNP/cluster_proportion_each_patient.xls"
write.table(proportion,myoutf,sep="\t",quote=F)


#################Perform GSEA########
myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_only_Human_nCountRNA_modified/Seruat_Myeloid_PW_Adult_only.integrated__finished.RDS"
Myeloid_PW_Adult_only.integrated<-readRDS(myinf)
tar<-"integrated_snn_res.0.6"
 Idents(object = Myeloid_PW_Adult_only.integrated) <- tar
 tar1<-as.character(unique(Myeloid_PW_Adult_only.integrated@meta.data$integrated_snn_res.0.6))
 tar1<-c('0','1','2','3','4','5','6','7','8','9','10')
for(i in 1:length(tar1))
{
  cluster.markers <- FindMarkers(Myeloid_PW_Adult_only.integrated, ident.1 =tar1[i],logfc.threshold = 0.001,min.pct = 0.001)
  myoutf<-paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_only_Human_nCountRNA_modified/PCx/res.0.6/GSEA/marker_integrated/Cluster_",tar1[i],"differentially_expressed_markers.xls")
  write.table(cluster.markers,myoutf,quote=F,sep="\t")
}


dir <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_only_Human_nCountRNA_modified/PCx/res.0.6/GSEA/marker/"
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
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_only_Human_nCountRNA_modified/PCx/res.0.6/GSEA/preranked_GSEA/",file_nam[i],".rnk")
  write.table(res,file=myoutf,quote=F,sep="\t",row.names=F)
}



##########Spliced versus unspliced GATA6########################
 myinf1<-paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_only_Human_nCountRNA_modified/Velocity/splice/spliced.csv")
 spliced<-read.csv(myinf1,sep='\t',header=TRUE)
 row.names(spliced)<-spliced$barcode
 spliced<-spliced[,2:ncol(spliced)]
 spliced_sle<-spliced[,c('SELP','GATA6','LYVE1','MRC1')]
 colnames(spliced_sle)<-paste0(colnames(spliced_sle),'_spliced')

 myinf2<-paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_only_Human_nCountRNA_modified/Velocity/splice/unspliced.csv")
 unspliced<-read.csv(myinf2,sep='\t',header=TRUE)
 row.names(unspliced)<-unspliced$barcode
 unspliced<-unspliced[,2:ncol(unspliced)]
 unspliced_sle<-unspliced[,c('SELP','GATA6','LYVE1','MRC1')]
 colnames(unspliced_sle)<-paste0(colnames(unspliced_sle),'_unspliced')


 transcript_all<-cbind(spliced_sle,unspliced_sle)
 
myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_only_Human_nCountRNA_modified/Seruat_Myeloid_PW_Adult_only.integrated__finished.RDS"
Myeloid_PW_Adult_only.integrated<-readRDS(myinf)

tag1<-GetAssayData(Myeloid_PW_Adult_only.integrated,assay='RNA',slot='data')['GATA6',]>0
cell_ID<-colnames(Myeloid_PW_Adult_only.integrated)[tag1]
transcript_all[,'GATA6_status']<-rep('Negative',nrow(transcript_all))
transcript_all[cell_ID,'GATA6_status']<-'Positive'
transcript_all[,'patient']<-Myeloid_PW_Adult_only.integrated@meta.data[row.names(transcript_all),'orig.ident']
transcript_all[,'cluster']<-Myeloid_PW_Adult_only.integrated@meta.data[row.names(transcript_all),'integrated_snn_res.0.6']

GATA6_splice_ratio_each_patient<-matrix(0,6,4)
row.names(GATA6_splice_ratio_each_patient)<-unique(transcript_all$patient)
colnames(GATA6_splice_ratio_each_patient)<-c('GATA6_pos_cell_splice','GATA6_pos_cell_unsplice','GATA6_neg_cell_splice','GATA6_neg_cell_unsplice')

tar<-unique(transcript_all$patient)
for (i in 1:length(tar))
{
  transcript_all_use<-transcript_all[transcript_all$patient==tar[i],]
  GATA6_pos<-transcript_all_use[transcript_all_use$GATA6_status=='Positive',]
  GATA6_neg<-transcript_all_use[transcript_all_use$GATA6_status=='Negative',]
  
  GATA6_pos_trans_all<-sum(GATA6_pos$GATA6_spliced[!is.na(GATA6_pos$GATA6_spliced)])+sum(GATA6_pos$GATA6_unspliced[!is.na(GATA6_pos$GATA6_unspliced)])
  GATA6_pos_splice_ratio<-sum(GATA6_pos$GATA6_spliced[!is.na(GATA6_pos$GATA6_spliced)])/GATA6_pos_trans_all

  GATA6_neg_trans_all<-sum(GATA6_neg$GATA6_spliced)+sum(GATA6_neg$GATA6_unspliced)
  GATA6_neg_splice_ratio<-sum(GATA6_neg$GATA6_spliced)/GATA6_neg_trans_all
  
  GATA6_splice_ratio_each_patient[tar[i],'GATA6_pos_cell_splice']<-GATA6_pos_splice_ratio
  GATA6_splice_ratio_each_patient[tar[i],'GATA6_pos_cell_unsplice']<-1-GATA6_pos_splice_ratio
  GATA6_splice_ratio_each_patient[tar[i],'GATA6_neg_cell_splice']<-GATA6_neg_splice_ratio
  GATA6_splice_ratio_each_patient[tar[i],'GATA6_neg_cell_unsplice']<-1-GATA6_neg_splice_ratio

}


myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revisiFon/Myeloid_PW_Adult_only_Human_nCountRNA_modified/Velocity/splice/GATA6_splice_ratio_each_patient.xls")
write.table(GATA6_splice_ratio_each_patient,myoutf,sep="\t",quote=F)













