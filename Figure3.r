#####Read all the data###########
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

xx[["percent.mt"]] <- PercentageFeatureSet(xx, pattern = "^MT-")
PW_Adult_kit.list<-list(PW_C,Kid1,Kid2,Kid4,Kid6,Kid7,Kid8,Kid9,Kid10,Kid11,PW_Q,PW_S,PW_U,PW_W,PW_Adu1,PW_Adu2)
PW_Adult_kit.list <- lapply(X = PW_Adult_kit.list, FUN = function(x) {
    x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")
    x <- subset(x, subset = nCount_RNA<40000 & nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})

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
 PW_Adult_kit.integrated<-FindNeighbors(PW_Adult_kit.integrated, dims = x)
 PW_Adult_kit.integrated <- FindClusters(PW_Adult_kit.integrated,resolution = seq(0.1,1.0,0.1))
 PW_Adult_kit.integrated <- RunUMAP(PW_Adult_kit.integrated, dims= x)

gender_info<-c('F','M','F','M','F','M','M','M','F','F','F','F','M','F','F','M')
age_info<-c('24','7','10','1.3','1.2','14','15','8','17','7','29','49','30','36','36','39')
tar<-unique(PW_Adult_kit.integrated@meta.data$orig.ident)
PW_Adult_kit.integrated@meta.data[,'gender']<-rep('Unknown',nrow(PW_Adult_kit.integrated@meta.data))
PW_Adult_kit.integrated@meta.data[,'age']<-rep('Unknown',nrow(PW_Adult_kit.integrated@meta.data))

PW_Adult_kit.integrated@meta.data[,'group']<-rep('Unassigned',nrow(PW_Adult_kit.integrated@meta.data))
tag1<-as.numeric(PW_Adult_kit.integrated@meta.data$age)>18&PW_Adult_kit.integrated@meta.data$gender=='M'
tag2<-as.numeric(PW_Adult_kit.integrated@meta.data$age)>18&PW_Adult_kit.integrated@meta.data$gender=='F'
tag3<-as.numeric(PW_Adult_kit.integrated@meta.data$age)<=18&PW_Adult_kit.integrated@meta.data$gender=='M'
tag4<-as.numeric(PW_Adult_kit.integrated@meta.data$age)<=18&PW_Adult_kit.integrated@meta.data$gender=='F'
PW_Adult_kit.integrated@meta.data$group[tag1]<-'Adult_Male'
PW_Adult_kit.integrated@meta.data$group[tag2]<-'Adult_Female'
PW_Adult_kit.integrated@meta.data$group[tag3]<-'Kid_Male'
PW_Adult_kit.integrated@meta.data$group[tag4]<-'Kid_Female'

for (i in 1:length(tar))
{
  tag<-PW_Adult_kit.integrated@meta.data$orig.ident==tar[i]
  PW_Adult_kit.integrated@meta.data[tag,'gender']<-gender_info[i]
  PW_Adult_kit.integrated@meta.data[tag,'age']<-age_info[i]
}

 myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/Seruat_PW_Adult_kit.integrated__finished.RDS"
 saveRDS(PW_Adult_kit.integrated, file = myoutf)

########Downstream Analysis#############
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
PW_Adult_kit.integrated@meta.data[,'GATA6']<-rep('0',nrow(PW_Adult_kit.integrated@meta.data))

tag1<-data['SELP',]>0
tag2<-data['TIMD4',]>0
tag3<-data['GATA6',]>0

name_1<-colnames(data)[tag1]
name_2<-colnames(data)[tag2]
name_3<-colnames(data)[tag3]

PW_Adult_kit.integrated@meta.data[name_1,'SELP']<-1
PW_Adult_kit.integrated@meta.data[name_2,'TIMD4']<-1
PW_Adult_kit.integrated@meta.data[name_3,'GATA6']<-1

SELP<-apply(table(PW_Adult_kit.integrated@meta.data$SELP,PW_Adult_kit.integrated@meta.data$orig.ident),2,function(x){x*100/sum(x)})
TIMD4<-apply(table(PW_Adult_kit.integrated@meta.data$TIMD4,PW_Adult_kit.integrated@meta.data$orig.ident),2,function(x){x*100/sum(x)})
GATA6<-apply(table(PW_Adult_kit.integrated@meta.data$GATA6,PW_Adult_kit.integrated@meta.data$orig.ident),2,function(x){x*100/sum(x)})

SELP<-table(PW_Adult_kit.integrated@meta.data$SELP,PW_Adult_kit.integrated@meta.data$orig.ident)
TIMD4<-table(PW_Adult_kit.integrated@meta.data$TIMD4,PW_Adult_kit.integrated@meta.data$orig.ident)
GATA6<-table(PW_Adult_kit.integrated@meta.data$GATA6,PW_Adult_kit.integrated@meta.data$orig.ident)

myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.6/cluster_patient_distribution/based_on_expression/SELP.xls"
 write.table(SELP,myoutf,sep="\t",quote=F)
 
myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.6/cluster_patient_distribution/based_on_expression/TIMD4.xls"
 write.table(TIMD4,myoutf,sep="\t",quote=F)
 
myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.6/cluster_patient_distribution/based_on_expression/GATA6.xls"
 write.table(GATA6,myoutf,sep="\t",quote=F)
 
########Feature plot#################
myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/Seruat_PW_Adult_kit.integrated__finished.RDS"
PW_Adult_kit.integrated<-readRDS(myinf)
DefaultAssay(PW_Adult_kit.integrated) <- "RNA"

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

#############Cluster contribution by each different patient##########
   
Idents(PW_Adult_kit.integrated)<-"integrated_snn_res.0.6"
 
 xx<-table(PW_Adult_kit.integrated@meta.data$integrated_snn_res.0.6,PW_Adult_kit.integrated@meta.data$orig.ident)
 yy<-t(apply(xx,2,function(x){x/sum(x)})   )
 myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.6/cluster_patient_distribution/cluster_proportion_CD45.xls"
 write.table(yy,myoutf,sep="\t",quote=F)
 
 myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.6/cluster_patient_distribution/cluster_counts_all.xls"
 write.table(xx,myoutf,sep="\t",quote=F)
 

myeloid<-c('0','1','2','4','7','11','15','17','18','19','21')
xx_myeloid<-xx[myeloid,]
yy_myeloid<-t(apply(xx_myeloid,2,function(x){x/sum(x)})   )
 myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.6/cluster_patient_distribution/cluster_proportion_myeloid.xls"
 write.table(yy_myeloid,myoutf,sep="\t",quote=F)
 
myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.6/cluster_patient_distribution/cluster_counts_myeloid.xls"
 write.table(xx_myeloid,myoutf,sep="\t",quote=F)
 
#################Output cell of each orig.ident and adult versus kid#############
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
 
  
###############DEG Cluster 1 vs Cluster 4 CD14+ vs CD14- CD1c+################
DefaultAssay(PW_Adult_kit.integrated) <- "integrated"
tar1 <- c("integrated_snn_res.0.6")
 
  Idents(object = PW_Adult_kit.integrated) <- tar1
  tmp <- PW_Adult_kit.integrated
  tmp.markers <- FindMarkers(object = tmp, only.pos = FALSE, ident.1='1',ident.2='4',min.pct = 0.1, logfc.threshold = 0.1)
  myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/PW_Adult_kid_normal_Human/PCx/res.0.6/cluster1_vs_cluster4.xls"
  write.table(tmp.markers,myoutf,sep="\t",quote=F)
