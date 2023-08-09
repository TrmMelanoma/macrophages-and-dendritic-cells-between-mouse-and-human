##Read all human adult data
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
######Read all B6 peritoneal wash data
myinf12 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/raw/Xin_Mu_PW_1/"
PW_Mu1.data <- Read10X(data.dir=myinf12)
PW_Mu1 <- CreateSeuratObject(counts = PW_Mu1.data, project = "Xin_B6_N1",min.cells=3,min.features=200)
counts <- GetAssayData(PW_Mu1, assay = "RNA")
Mu_genes<-row.names(counts)
Mouse_to_human <- function(x){

 human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
 mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

   genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T,verbose=F)

   humanx <- unique(genesV2[, 2])

   return(genesV2)
}
Hu_genes<-Mouse_to_human(Mu_genes)
tag1<-duplicated(Hu_genes[,1])
Hu_genes<-Hu_genes[!tag1,]
tag2<-duplicated(Hu_genes[,2])
Hu_genes<-Hu_genes[!tag2,]
  counts <- GetAssayData(PW_Mu1, assay = "RNA")
  counts <- counts[(which(rownames(counts) %in% Hu_genes$MGI.symbol)),]
  PW_Mu1 <- subset(PW_Mu1, features = rownames(counts))

for (i in 1:length(PW_Mu1@assays$RNA@counts@Dimnames[[1]]))
{
 cat("\r",i)
 mu_gene<-PW_Mu1@assays$RNA@counts@Dimnames[[1]][i] 
 hu_gene<-Hu_genes[Hu_genes$MGI.symbol==mu_gene,'HGNC.symbol']
 PW_Mu1@assays$RNA@counts@Dimnames[[1]][i]<-hu_gene
}
for (i in 1:length(PW_Mu1@assays$RNA@data@Dimnames[[1]]))
{
 cat("\r",i)
 mu_gene<-PW_Mu1@assays$RNA@data@Dimnames[[1]][i] 
 hu_gene<-Hu_genes[Hu_genes$MGI.symbol==mu_gene,'HGNC.symbol']
 PW_Mu1@assays$RNA@data@Dimnames[[1]][i]<-hu_gene
}

PW_Mu1<-RenameCells(object=PW_Mu1,add.cell.id=unique(PW_Mu1@meta.data$orig.ident))
PW_Mu1<-RenameCells(object=PW_Mu1,new.names=gsub('-1','',Cells(PW_Mu1)))


myinf13 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/raw/Xin_Mu_PW_3/"
PW_Mu3.data <- Read10X(data.dir=myinf13)
PW_Mu3 <- CreateSeuratObject(counts = PW_Mu3.data, project = "Xin_B6_N3",min.cells=3,min.features=200)
counts <- GetAssayData(PW_Mu3, assay = "RNA")
Mu_genes<-row.names(counts)
Mouse_to_human <- function(x){

 human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
 mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

   genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T,verbose=F)

   humanx <- unique(genesV2[, 2])

   return(genesV2)
}
Hu_genes<-Mouse_to_human(Mu_genes)
tag1<-duplicated(Hu_genes[,1])
Hu_genes<-Hu_genes[!tag1,]
tag2<-duplicated(Hu_genes[,2])
Hu_genes<-Hu_genes[!tag2,]

  counts <- GetAssayData(PW_Mu3, assay = "RNA")
  counts <- counts[(which(rownames(counts) %in% Hu_genes$MGI.symbol)),]
  PW_Mu3 <- subset(PW_Mu3, features = rownames(counts))

for (i in 1:length(PW_Mu3@assays$RNA@counts@Dimnames[[1]]))
{
 cat("\r",i)
 mu_gene<-PW_Mu3@assays$RNA@counts@Dimnames[[1]][i] 
 hu_gene<-Hu_genes[Hu_genes$MGI.symbol==mu_gene,'HGNC.symbol']
 PW_Mu3@assays$RNA@counts@Dimnames[[1]][i]<-hu_gene
}
for (i in 1:length(PW_Mu3@assays$RNA@data@Dimnames[[1]]))
{
 cat("\r",i)
 mu_gene<-PW_Mu3@assays$RNA@data@Dimnames[[1]][i] 
 hu_gene<-Hu_genes[Hu_genes$MGI.symbol==mu_gene,'HGNC.symbol']
 PW_Mu3@assays$RNA@data@Dimnames[[1]][i]<-hu_gene
}

PW_Mu3<-RenameCells(object=PW_Mu3,add.cell.id=unique(PW_Mu3@meta.data$orig.ident))
PW_Mu3<-RenameCells(object=PW_Mu3,new.names=gsub('-1','',Cells(PW_Mu3)))





myinf17 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/raw/Xin_Mu_PW_2/"
PW_Mu2.data <- Read10X(data.dir=myinf17)
PW_Mu2 <- CreateSeuratObject(counts = PW_Mu2.data, project = "Xin_B6_N2",min.cells=3,min.features=200)
counts <- GetAssayData(PW_Mu2, assay = "RNA")
Mu_genes<-row.names(counts)
Mouse_to_human <- function(x){

 human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
 mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

   genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T,verbose=F)

   humanx <- unique(genesV2[, 2])

   return(genesV2)
}
Hu_genes<-Mouse_to_human(Mu_genes)
tag1<-duplicated(Hu_genes[,1])
Hu_genes<-Hu_genes[!tag1,]
tag2<-duplicated(Hu_genes[,2])
Hu_genes<-Hu_genes[!tag2,]

  counts <- GetAssayData(PW_Mu2, assay = "RNA")
  counts <- counts[(which(rownames(counts) %in% Hu_genes$MGI.symbol)),]
  PW_Mu2 <- subset(PW_Mu2, features = rownames(counts))

for (i in 1:length(PW_Mu2@assays$RNA@counts@Dimnames[[1]]))
{
 cat("\r",i)
 mu_gene<-PW_Mu2@assays$RNA@counts@Dimnames[[1]][i] 
 hu_gene<-Hu_genes[Hu_genes$MGI.symbol==mu_gene,'HGNC.symbol']
 PW_Mu2@assays$RNA@counts@Dimnames[[1]][i]<-hu_gene
}
for (i in 1:length(PW_Mu2@assays$RNA@data@Dimnames[[1]]))
{
 cat("\r",i)
 mu_gene<-PW_Mu2@assays$RNA@data@Dimnames[[1]][i] 
 hu_gene<-Hu_genes[Hu_genes$MGI.symbol==mu_gene,'HGNC.symbol']
 PW_Mu2@assays$RNA@data@Dimnames[[1]][i]<-hu_gene
}

PW_Mu2<-RenameCells(object=PW_Mu2,add.cell.id=unique(PW_Mu2@meta.data$orig.ident))
PW_Mu2<-RenameCells(object=PW_Mu2,new.names=gsub('-1','',Cells(PW_Mu2)))


#####Integration
xx<- merge(PW_normal, y = c(PW_Mu1,PW_Mu2,PW_Mu3,PW_Q,PW_S,PW_U,PW_W,PW_Adu1,PW_Adu2), project = "original")
PW_Adult_all_b6_all.list<-SplitObject(xx, split.by = "orig.ident")
 xx[["percent.mt"]] <- PercentageFeatureSet(xx, pattern = "^MT-")
PW_Adult_all_b6_all.list <- lapply(X = PW_Adult_all_b6_all.list, FUN = function(x) {
    x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")
    x <- subset(x, subset = nCount_RNA<60000 & nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 10)
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

 features <- SelectIntegrationFeatures(object.list = PW_Adult_all_b6_all.list,nfeatures=2000)
 PW_Adult_all_b6_all.list <- lapply(X = PW_Adult_all_b6_all.list, FUN = function(x) {
 x <- ScaleData(x, features = features, verbose = FALSE)
 x <- RunPCA(x, features = features, verbose = FALSE)
 })

 immune.anchors <- FindIntegrationAnchors(object.list = PW_Adult_all_b6_all.list, anchor.features = features, reduction = "rpca")
 PW_Adult_all_b6_all.integrated <- IntegrateData(anchorset = immune.anchors)
 PW_Adult_all_b6_all.integrated <- ScaleData(PW_Adult_all_b6_all.integrated, verbose = FALSE)
 PW_Adult_all_b6_all.integrated <- RunPCA(PW_Adult_all_b6_all.integrated, npcs = 30, verbose = FALSE)
 x<-c(1:15)
 
 PW_Adult_all_b6_all.integrated<-FindNeighbors(PW_Adult_all_b6_all.integrated, dims = x)
 PW_Adult_all_b6_all.integrated <- FindClusters(PW_Adult_all_b6_all.integrated,resolution = seq(0.1,1.0,0.1))

 PW_Adult_all_b6_all.integrated <- RunUMAP(PW_Adult_all_b6_all.integrated, dims= x)

 myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/PW_Adult_all_b6_all.integrated_analysis/RPCA/TSNE_PW_Adult_all_b6_all.integrated_groupby_orig.ident.pdf"
  pdf(myoutf,width=8,height=5)
 DimPlot(object = PW_Adult_all_b6_all.integrated,group.by='orig.ident',reduction="umap",pt.size=0.1,raster=FALSE)
 dev.off()

 myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/PW_Adult_all_b6_all.integrated_analysis/RPCA/Seruat_PW_Adult_all_b6_all.integrated__finished.RDS"
 saveRDS(PW_Adult_all_b6_all.integrated, file = myoutf)

#######Vlnplot show expression of key genes to subset myeloid cell cluster##############

Immune_general<-c('CD3G','CD3E','CD8A','CD4','NKG7','KLRB1','ITGAM','ITGAX','MKI67','CD1C','CD14','CD19','GATA6','ADGRE1','LYVE1','MRC1','XCR1','TCF4','IRF7','IRF8')
DefaultAssay(PW_Adult_all_b6_all.integrated) <- "RNA"
Idents(PW_Adult_all_b6_all.integrated)<-"integrated_snn_res.0.7"
 
myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/PW_Adult_all_b6_all.integrated_analysis/RPCA/vlnplot_subset_myeloid_clusters_res.0.7.pdf"
pdf(myoutf,width=50,height=60)
VlnPlot(PW_Adult_all_b6_all.integrated, features=Immune_general,ncol=4)
dev.off()


##################Subset myeloid cells only and downsample for downstream analysis#########################
 tar<-"integrated_snn_res.0.7"
 Idents(object = PW_Adult_all_b6_all.integrated) <- tar
 subcluster<-c('0','1','3','4','5','6','8','9','11','13','14','16','17','18','19','20','25','26','27','29')
 tag<-PW_Adult_all_b6_all.integrated@meta.data[,tar]%in%subcluster
 tag1<-GetAssayData(PW_Adult_all_b6_all.integrated,assay='RNA', slot='data')['CD3E',]>0.05|GetAssayData(PW_Adult_all_b6_all.integrated,assay='RNA', slot='data')['CD19',]
 cell_ID<-row.names(PW_Adult_all_b6_all.integrated@meta.data)[tag&!tag1]
 
 myinf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/PW_Adult_all_b6_all.integrated_analysis/Seruat_PW_Adult_all_b6_all.integrated__raw.RDS"
 xx<-readRDS(myinf)
 Myeloid_PW_Adult_all_b6_all.integrated<-subset(x=xx,cells=cell_ID,nCount_RNA>3000)
 table(Myeloid_PW_Adult_all_b6_all.integrated@meta.data$orig.ident)
  #PW_Adu1   PW_Adu2 PW_normal      PW_Q      PW_S      PW_U      PW_W Xin_B6_N1 
  #   2847      2420      5167      2554      3583      1115      4068     16226 
  #Xin_B6_N2 Xin_B6_N3 
  #  10682     12202

cell_ID<-character()
tar<-c('Xin_B6_N1','Xin_B6_N2','Xin_B6_N3')
for (i in 1:length(tar))
{
  each_cell_ID<-colnames(Myeloid_PW_Adult_all_b6_all.integrated)[Myeloid_PW_Adult_all_b6_all.integrated$orig.ident==tar[i]]
  tt<-sample(each_cell_ID,round(0.4*length(each_cell_ID),digits=0),replace=FALSE)
  cell_ID<-c(cell_ID,tt)
}
Human_cell_ID<-colnames(Myeloid_PW_Adult_all_b6_all.integrated)[Myeloid_PW_Adult_all_b6_all.integrated$orig.ident%in%c('PW_Adu1','PW_Adu2','PW_normal','PW_Q','PW_S','PW_U','PW_W')]
cell_ID<-c(cell_ID,Human_cell_ID)
Myeloid_PW_Adult_all_b6_all.integrated_use<-subset(Myeloid_PW_Adult_all_b6_all.integrated,cells=cell_ID)
table(Myeloid_PW_Adult_all_b6_all.integrated_use@meta.data$orig.ident)
#  PW_Adu1   PW_Adu2 PW_normal      PW_Q      PW_S      PW_U      PW_W Xin_B6_N1 
#     2847      2420      5167      2554      3583      1115      4068      6490 
#  Xin_B6_N2 Xin_B6_N3 
#     4273      4881 

 Myeloid_PW_Adult_all_b6_all.integrated.list <- SplitObject(Myeloid_PW_Adult_all_b6_all.integrated_use, split.by = "orig.ident")
 Myeloid_PW_Adult_all_b6_all.integrated.list <- lapply(X = Myeloid_PW_Adult_all_b6_all.integrated.list, FUN = function(x) {
    x <- NormalizeData(x,assay='RNA')
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000,assay='RNA')
 })

 immune.anchors <- FindIntegrationAnchors(object.list = Myeloid_PW_Adult_all_b6_all.integrated.list, dims = 1:20,anchor.features = 3000)
 Myeloid_PW_Adult_all_b6_all.integrated_CCA<-IntegrateData(anchorset = immune.anchors, dims = 1:20)
 DefaultAssay(Myeloid_PW_Adult_all_b6_all.integrated_CCA) <- "integrated"
 Myeloid_PW_Adult_all_b6_all.integrated_CCA <- ScaleData(Myeloid_PW_Adult_all_b6_all.integrated_CCA, verbose = FALSE)
 Myeloid_PW_Adult_all_b6_all.integrated_CCA <- RunPCA(Myeloid_PW_Adult_all_b6_all.integrated_CCA, npcs = 30, verbose = FALSE)
 x<-c(1:15)
 Myeloid_PW_Adult_all_b6_all.integrated_CCA<-FindNeighbors(Myeloid_PW_Adult_all_b6_all.integrated_CCA, dims = x)
  Myeloid_PW_Adult_all_b6_all.integrated_CCA <- FindClusters(Myeloid_PW_Adult_all_b6_all.integrated_CCA,resolution = seq(0.1,1.0,0.1))
 tar1 <- c("integrated_snn_res.0.6")
 DefaultAssay(Myeloid_PW_Adult_all_b6_all.integrated_CCA) <- "RNA"
 for(i in 1 : length(tar1))
 {
  cat("\r",i)
  Idents(object = Myeloid_PW_Adult_all_b6_all.integrated_CCA) <- tar1[i]
  tmp <- Myeloid_PW_Adult_all_b6_all.integrated_CCA
  tmp.markers <- FindAllMarkers(object = tmp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_all_b6_all.integrated_analysis/downsample/CCA/PCx/cluster_markers_",tar1[i],"_ScRNA_RNA.xls")
  write.table(tmp.markers,myoutf,sep="\t",quote=F)
  myoutf1 <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_all_b6_all.integrated_analysis/downsample/CCA/PCx/heatmap_markers_",tar1[i],"_RNA.pdf")
  pdf(myoutf1,30,35)
  top20 <- tmp.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
  print(DoHeatmap(tmp, features = top20$gene) + NoLegend())
  dev.off()
 }
 Myeloid_PW_Adult_all_b6_all.integrated_CCA <- RunUMAP(Myeloid_PW_Adult_all_b6_all.integrated_CCA, dims= x)
 for(i in 1:length(tar))
 {
  Idents(object = Myeloid_PW_Adult_all_b6_all.integrated_CCA) <- tar[i]
  myoutf1 <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_all_b6_all.integrated_analysis/downsample/CCA/UMAP_SB6_",tar[i],"_undefined.pdf")
  pdf(myoutf1,7,5)
  print(DimPlot(Myeloid_PW_Adult_all_b6_all.integrated_CCA, reduction = "umap"))
  dev.off()
 }
moutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_all_b6_all.integrated_analysis/downsample/CCA/Seruat_PW_Adult_all_b6_all.integrated_CCA__finished.RDS"
 saveRDS(Myeloid_PW_Adult_all_b6_all.integrated_CCA, file = myoutf)

#############ClusterTree###################################

myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_all_b6_all.integrated_analysis/downsample/CCA/Seruat_PW_Adult_all_b6_all.integrated_CCA__finished.RDS"
Myeloid_PW_Adult_all_b6_all.integrated_CCA<-readRDS(myinf)

 myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_all_b6_all.integrated_analysis/downsample/CCA/PCx/clustertree/cluster_tree_res_0.1-1.pdf"
  pdf(myoutf,width=8,height=12)
  clustree(Myeloid_PW_Adult_all_b6_all.integrated_CCA, prefix = "integrated_snn_res.")
 dev.off()


##############The best resolution is 0.6 and recolor the clusters####################

Idents(object = Myeloid_PW_Adult_all_b6_all.integrated_CCA) <- 'integrated_snn_res.0.6'
levels(Myeloid_PW_Adult_all_b6_all.integrated_CCA@meta.data$integrated_snn_res.0.6)<-seq(0,17,1)
  myoutf1 <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_all_b6_all.integrated_analysis/downsample/CCA/PCx/res.0.6/UMAP_recolored.pdf"
  pdf(myoutf1,7,5)
  print(DimPlot(Myeloid_PW_Adult_all_b6_all.integrated_CCA, reduction = "umap",cols=c("#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231",
                                                                                      "#911eb4", "#46f0f0", "#f032e6", "#bcf60c", "#fabebe", "#008080",
                                                                                       "#e6beff", "#9a6324", "#fffac8", "#800000", "#aaffc3", "#808000", 
                                                                                       "#ffd8b1")))
  dev.off()


#######Calculate the DEGS between two CD1C clusters########

myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_all_b6_all.integrated_analysis/downsample/CCA/Seruat_PW_Adult_all_b6_all.integrated_CCA__finished.RDS"
Myeloid_PW_Adult_all_b6_all.integrated_CCA<-readRDS(myinf)
DefaultAssay(Myeloid_PW_Adult_all_b6_all.integrated_CCA) <- "RNA"
Idents(Myeloid_PW_Adult_all_b6_all.integrated_CCA)<-'integrated_snn_res.0.6'
 cluster.markers <- FindMarkers(Myeloid_PW_Adult_all_b6_all.integrated_CCA, ident.1 ='6',ident.2='8',logfc.threshold = 0.2,min.pct = 0.2)
  myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_all_b6_all.integrated_analysis/downsample/CCA/PCx/Cluster6_vs_9.xls"
  write.table(cluster.markers,myoutf,quote=F,sep="\t")

#########Perform heatmap for CD1C two clusters##################
exp<-Myeloid_PW_Adult_all_b6_all.integrated_CCA
DefaultAssay(exp) <- "RNA"
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
res <- as.data.frame(res)
data<-res
cluster.markers<-cluster.markers[order(cluster.markers$avg_log2FC),]
DEG_C1_C4<-c(row.names(cluster.markers)[1:50],tail(row.names(cluster.markers),50))
DEG_C1_C4<-c('CXCL9','TREM2','CX3CR1','CD84','CD163',DEG_C1_C4)
cells<-row.names(info)[info$integrated_snn_res.0.6%in%c('6','8')]
data_DEG_C1_C4<-data[DEG_C1_C4,cells]
data_DEG_C1_C4<<-na.omit(data_DEG_C1_C4)
info_heatmap<-info[cells,]
info_heatmap[,'heatmap_group']<-paste0(info_heatmap$integrated_snn_res.0.6,'_',info_heatmap$orig.ident)
tar <- unique(info_heatmap$heatmap_group)
#tar <- as.numeric(tar)
tar <- tar[order(tar)]
for(i in 1 : length(tar))
{
  res <- matrix(0, nrow(data_DEG_C1_C4),4)
  row.names(res) <- row.names(data_DEG_C1_C4)
  colnames(res) <- c("C1_avg_expression","C1_median_expression","C2_avg_expression","C2_median_expression")
  
  tag <- info_heatmap$heatmap_group == tar[i]
  
  sam1 <- row.names(info_heatmap)[tag==1]
  sam2 <- row.names(info_heatmap)[tag==0]
  
  for(k in 1 : nrow(data_DEG_C1_C4))
  {
    cat("\r",i,"----->",k)
    xx1 <- as.numeric(data_DEG_C1_C4[k, sam1])
    xx2 <- as.numeric(data_DEG_C1_C4[k, sam2])
    
    res[k,"C1_avg_expression"] <- mean(xx1)
    res[k,"C2_avg_expression"] <- mean(xx2)
    res[k,"C1_median_expression"]<-median(xx1)
    res[k,"C2_median_expression"] <- median(xx2)
  }
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_all_b6_all.integrated_analysis/downsample/CCA/PCx/res.0.6/pheatmap_CD1C/Cluster6_vs_8_Heatmap_",tar[i],"_vs_others_normalized_expression.xls")
  write.table(res,myoutf,sep="\t",quote=F)
}

res<-matrix(0,nrow(data_DEG_C1_C4),length(tar))
 row.names(res)<-row.names(data_DEG_C1_C4)
 colnames(res)<-tar[order(tar)]
for (i in 1:length(tar))
{
  myinf<-paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_all_b6_all.integrated_analysis/downsample/CCA/PCx/res.0.6/pheatmap_CD1C/Cluster6_vs_8_Heatmap_",tar[i],"_vs_others_normalized_expression.xls")
  xx<-read.table(myinf,sep='\t')
  res[,i]<-xx$C1_avg_expression
}
tag<-grep('PW',colnames(res))
res<-res[,tag]
res <- apply(res,1, function(arg) (arg-mean(arg))/sd(arg))
res <- t(res)
breaklist_1<-seq(-2.5,0,0.2)
breaklist_2<-seq(0,2.6,0.2)
breaklist<-unique(c(breaklist_1,breaklist_2))

myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_all_b6_all.integrated_analysis/downsample/CCA/PCx/res.0.6/pheatmap_CD1C/Cluster6_vs_cluster8_Heatmap_gene_each_cluster.pdf")
pdf(myoutf,width=20,height=30)
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


############Calculate cluster proportion of each sample and compare mouse and human##################
tar<-unique(Myeloid_PW_Adult_all_b6_all.integrated_CCA@meta.data$orig.ident)
proportion<-table(Myeloid_PW_Adult_all_b6_all.integrated_CCA@meta.data$orig.ident,Myeloid_PW_Adult_all_b6_all.integrated_CCA@meta.data$integrated_snn_res.0.6)
proportion<-as.matrix(proportion)
proportion<-as.data.frame(t(apply(proportion,1,function(x){x/sum(x)})))
proportion[,'GATA6']<-rep(0,nrow(proportion))
exp<-t(as.matrix(GetAssayData(Myeloid_PW_Adult_all_b6_all.integrated_CCA,assay='RNA',slot='data')))
for (i in 1:length(tar))
{
  tag<-Myeloid_PW_Adult_all_b6_all.integrated_CCA@meta.data$orig.ident==tar[i]
  cells<-row.names(Myeloid_PW_Adult_all_b6_all.integrated_CCA@meta.data)[tag]
  exp_use<-exp[cells,]
  proportion[tar[i],'GATA6']<-sum(exp_use[,'GATA6']>0)/nrow(exp_use)
}
myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_all_b6_all.integrated_analysis/downsample/CCA/PCx/res.0.6/cluster_proportion_MNP/cluster_proportion_each_patient.xls"
write.table(proportion,myoutf,sep="\t",quote=F)


######Featureplot#########
myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_all_b6_all.integrated_analysis/downsample/CCA/Seruat_Myeloid_PW_Adult_all_b6_all.integrated__finished.RDS"
Myeloid_PW_Adult_all_b6_all.integrated_CCA<-readRDS(myinf)
DefaultAssay(Myeloid_PW_Adult_all_b6_all.integrated_CCA) <- "RNA"
Figure2<-c('GATA6','SELP','TIMD4','CCR2','MRC1','LYVE1','CD14','CD1C','CD226')
  tag1<-grep('PW',Myeloid_PW_Adult_all_b6_all.integrated_CCA@meta.data$orig.ident)
  tag2<-grep('Xin',Myeloid_PW_Adult_all_b6_all.integrated_CCA@meta.data$orig.ident)
    cell_ID1<-row.names(Myeloid_PW_Adult_all_b6_all.integrated_CCA@meta.data)[tag1]
    cell_ID2<-row.names(Myeloid_PW_Adult_all_b6_all.integrated_CCA@meta.data)[tag2]
p<-FeaturePlot(object = Myeloid_PW_Adult_all_b6_all.integrated_CCA, cells=cell_ID1, features = Figure2,cols= c('#deebf7', "#bd0026"),
               reduction = "umap",slot = "data",max.cutoff=4,pt.size=0.2,combine=F)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + theme(plot.title = element_text(size=40),
                          axis.title = element_text(family="Helvetica",face = "bold", color = "black",size=20), 
                           axis.text=element_text(family="Helvetica",face = "bold", color = "black",size=16),
                           axis.line = element_line(colour = 'black', size = 1),
                           axis.ticks = element_line(colour = "black", size = 1),
                           axis.ticks.length=unit(.3, "cm"))
}
myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_all_b6_all.integrated_analysis/downsample/CCA/PCx/Figure2_marker_human.pdf"
pdf(myoutf,width=35,height=30)
print(cowplot::plot_grid(plotlist = p,ncol=3))
dev.off()

p<-FeaturePlot(object = Myeloid_PW_Adult_all_b6_all.integrated_CCA, cells=cell_ID2, features = Figure2,cols= c('#deebf7', "#bd0026"),
               reduction = "umap",slot = "data",max.cutoff=4,pt.size=0.2,combine=F)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + theme(plot.title = element_text(size=40),
                          axis.title = element_text(family="Helvetica",face = "bold", color = "black",size=20), 
                           axis.text=element_text(family="Helvetica",face = "bold", color = "black",size=16),
                           axis.line = element_line(colour = 'black', size = 1),
                           axis.ticks = element_line(colour = "black", size = 1),
                           axis.ticks.length=unit(.3, "cm"))
}
 
myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_all_b6_all.integrated_analysis/downsample/CCA/PCx/Figure2_marker_mouse.pdf"
pdf(myoutf,width=35,height=30)
print(cowplot::plot_grid(plotlist = p,ncol=3))
dev.off()

#################Output cell of each orig.ident#############
Idents(Myeloid_PW_Adult_all_b6_all.integrated_CCA)<-"integrated_snn_res.0.6"
 tar<-unique(Myeloid_PW_Adult_all_b6_all.integrated_CCA@meta.data$orig.ident)
 for (i in 1:length(tar))
 {
    tag<-Myeloid_PW_Adult_all_b6_all.integrated_CCA@meta.data$orig.ident==tar[i]
    cell_ID<-row.names(Myeloid_PW_Adult_all_b6_all.integrated_CCA@meta.data)[tag]
    myoutf1 <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_all_b6_all.integrated_analysis/downsample/CCA/PCx/res.0.6/each_sample/UMAP_SB6_",tar[i],"_undefined.pdf")
  pdf(myoutf1,7,5)
  print(DimPlot(Myeloid_PW_Adult_all_b6_all.integrated_CCA, reduction = "umap",pt.size=0.1,cells =cell_ID))
  dev.off()

 }
  
  tag1<-grep('PW',Myeloid_PW_Adult_all_b6_all.integrated_CCA@meta.data$orig.ident)
  tag2<-grep('Xin',Myeloid_PW_Adult_all_b6_all.integrated_CCA@meta.data$orig.ident)
  
    cell_ID1<-row.names(Myeloid_PW_Adult_all_b6_all.integrated_CCA@meta.data)[tag1]
    cell_ID2<-row.names(Myeloid_PW_Adult_all_b6_all.integrated_CCA@meta.data)[tag2]
  
  Myeloid_PW_Adult_all_b6_all.integrated_CCA@meta.data[,'species']<-rep('Human',nrow(Myeloid_PW_Adult_all_b6_all.integrated_CCA@meta.data))
  Myeloid_PW_Adult_all_b6_all.integrated_CCA@meta.data$species[tag2]<-'Mouse'

  xx_hu<-subset(Myeloid_PW_Adult_all_b6_all.integrated_CCA,subset = species == 'Human')

  myoutf1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_all_b6_all.integrated_analysis/downsample/CCA/PCx/res.0.6/each_sample/all_human.pdf"
  pdf(myoutf1,7,5)
  print(DimPlot(Myeloid_PW_Adult_all_b6_all.integrated_CCA,cells =cell_ID1, reduction = "umap",pt.size=0.1,cols=c("#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231",
                                                                                      "#911eb4", "#46f0f0", "#f032e6", "#bcf60c", "#fabebe", "#008080",
                                                                                       "#e6beff", "#9a6324", "#fffac8", "#800000", "#aaffc3", "#808000", 
                                                                                       "#ffd8b1")))
  dev.off()

myoutf1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_all_b6_all.integrated_analysis/downsample/CCA/PCx/res.0.6/each_sample/all_mouse.pdf"
  pdf(myoutf1,7,5)
  print(DimPlot(Myeloid_PW_Adult_all_b6_all.integrated_CCA, reduction = "umap",pt.size=0.1,cells =cell_ID2,cols=c("#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231",
                                                                                      "#911eb4", "#46f0f0", "#f032e6", "#bcf60c", "#fabebe", "#008080",
                                                                                       "#e6beff", "#9a6324", "#fffac8", "#800000", "#aaffc3", "#808000", 
                                                                                       "#ffd8b1")))
  dev.off()




#############Look at what are the clusters in adult only analysis are now in this new space##########
myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_only_Human_nCountRNA_modified/Seruat_Myeloid_PW_Adult_only.integrated__finished.RDS"
Myeloid_PW_Adult_only.integrated<-readRDS(myinf)

myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_all_b6_all.integrated_analysis/downsample/CCA/Seruat_PW_Adult_all_b6_all.integrated_CCA__finished.RDS"
Myeloid_PW_Adult_all_b6_all.integrated_CCA<-readRDS(myinf)

Myeloid_PW_Adult_all_b6_all.integrated_CCA@meta.data[,'previous_object']<-rep('NA',nrow(Myeloid_PW_Adult_all_b6_all.integrated_CCA@meta.data))
Myeloid_PW_Adult_all_b6_all.integrated_CCA@meta.data[row.names(Myeloid_PW_Adult_only.integrated@meta.data),'previous_object']<-as.character(Myeloid_PW_Adult_only.integrated@meta.data$integrated_snn_res.0.6)
Myeloid_PW_Adult_all_b6_all.integrated_CCA$previous_object<-factor(Myeloid_PW_Adult_all_b6_all.integrated_CCA$previous_object,levels=c(seq(0,15,1),'NA'))
Idents(Myeloid_PW_Adult_all_b6_all.integrated_CCA)='previous_object'
myoutf1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_all_b6_all.integrated_analysis/downsample/CCA/PCx/res.0.6/previous_adult_only_UMAP.pdf"
  pdf(myoutf1,7,5)
  print(DimPlot(Myeloid_PW_Adult_all_b6_all.integrated_CCA, reduction = "umap",pt.size=0.1,cols=c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
  "#8dd3c7", "#bebada", "#fb8072", "#80b1d3", "#fdb462",
  "#b3de69","#ececec")))
  dev.off()



#######Dotplot show GATA6+ clusters###########

myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_all_b6_all.integrated_analysis/downsample/CCA/Seruat_PW_Adult_all_b6_all.integrated_CCA__finished.RDS"
Myeloid_PW_Adult_all_b6_all.integrated_CCA<-readRDS(myinf)

LCM_markers<-c('GATA6','SELP','ADGRE1','LYVE1','MRC1','CCR2')

myoutf <-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_all_b6_all.integrated_analysis/downsample/CCA/PCx/res.0.6/Dot_plot_LCM_vs_CM_macrophage.pdf"
pdf(myoutf,width=4.5,height=4)
print(DotPlot(Myeloid_PW_Adult_all_b6_all.integrated_CCA, cols=c("#08519c","#ef3b2c"),group.by="integrated_snn_res.0.6",
              features = LCM_markers,assay='RNA',idents = c('0','1','2','3','4','5','7','9','10','11'))+
        theme(axis.text.x = element_text(angle = 45,size=10,hjust = 1)))

dev.off()



################Cell cycle scoring###################
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

myinf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_all_b6_all.integrated_analysis/downsample/CCA/Seruat_PW_Adult_all_b6_all.integrated_CCA__finished.RDS"
Myeloid_PW_Adult_all_b6_all.integrated_CCA<-readRDS(myinf)
DefaultAssay(Myeloid_PW_Adult_all_b6_all.integrated_CCA)<-'RNA'
Myeloid_PW_Adult_all_b6_all.integrated_CCA<-CellCycleScoring(Myeloid_PW_Adult_all_b6_all.integrated_CCA, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
cell_cycle_score<-t(Myeloid_PW_Adult_all_b6_all.integrated_CCA@meta.data[,c("S.Score",'G2M.Score')])
adt_assay <- CreateAssayObject(counts = cell_cycle_score)
Myeloid_PW_Adult_all_b6_all.integrated_CCA[['cell_cycle']]<-adt_assay
DefaultAssay(Myeloid_PW_Adult_all_b6_all.integrated_CCA) <- "cell_cycle"
  tag1<-grep('PW',Myeloid_PW_Adult_all_b6_all.integrated_CCA@meta.data$orig.ident)
  tag2<-grep('Xin',Myeloid_PW_Adult_all_b6_all.integrated_CCA@meta.data$orig.ident)

    cell_ID1<-row.names(Myeloid_PW_Adult_all_b6_all.integrated_CCA@meta.data)[tag1]
    cell_ID2<-row.names(Myeloid_PW_Adult_all_b6_all.integrated_CCA@meta.data)[tag2]
  

myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_all_b6_all.integrated_analysis/downsample/CCA/By_cell_cycle_Myeloid_PW_Adult_all_b6_all.integrated_CCA_Human.pdf"
  pdf(myoutf,width=6,height=5)
 DimPlot(object = Myeloid_PW_Adult_all_b6_all.integrated_CCA,cells=cell_ID1,group.by='Phase',reduction="umap",pt.size=0.3)
 dev.off()


myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_all_b6_all.integrated_analysis/downsample/CCA/By_cell_cycle_Myeloid_PW_Adult_all_b6_all.integrated_CCA_Mouse.pdf"
  pdf(myoutf,width=6,height=5)
 DimPlot(object = Myeloid_PW_Adult_all_b6_all.integrated_CCA,cells=cell_ID2,group.by='Phase',reduction="umap",pt.size=0.3)
 dev.off()

xx<-table(Myeloid_PW_Adult_all_b6_all.integrated_CCA$orig.ident,Myeloid_PW_Adult_all_b6_all.integrated_CCA$Phase)
xx<-apply(xx,1,function(x){x/sum(x)})
myoutf<-"/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_all_b6_all.integrated_analysis/downsample/CCA/PCx/res.0.6/cell_cycle_statistics.xls"
write.table(xx,myoutf,sep="\t",quote=F)

########Perform pre-rank GSEA###############################
myinf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_all_b6_all.integrated_analysis/downsample/CCA/Seruat_Myeloid_PW_Adult_all_b6_all.integrated_CCA__finished.RDS"
 Myeloid_PW_Adult_all_b6_all.integrated_CCA<-readRDS(myinf)
tar<-"integrated_snn_res.0.6"
 Idents(object = Myeloid_PW_Adult_all_b6_all.integrated_CCA) <- tar
 tar1<-as.character(unique(Myeloid_PW_Adult_all_b6_all.integrated_CCA@meta.data$integrated_snn_res.0.6))
for(i in 1:length(tar1))
{
  cluster.markers <- FindMarkers(Myeloid_PW_Adult_all_b6_all.integrated_CCA, ident.1 =tar1[i],logfc.threshold = 0.001,min.pct = 0.001)
  myoutf<-paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_all_b6_all.integrated_analysis/downsample/CCA/PCx/res.0.6/GSEA/marker_integrated/Cluster_",tar1[i],"differentially_expressed_markers.xls")
  write.table(cluster.markers,myoutf,quote=F,sep="\t")
}

dir <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_all_b6_all.integrated_analysis/downsample/CCA/PCx/res.0.6/GSEA/marker_integrated/"
files <- list.files(dir)
file_nam <- gsub(".xls","",files)
for(i in 1:length(files))
{
  cat("\r",i)
  tmpinf <- paste0(dir,files[i])
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
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_all_b6_all.integrated_analysis/downsample/CCA/PCx/res.0.6/GSEA/preranked_GSEA_integrated/",file_nam[i],".rnk")
  write.table(res,file=myoutf,quote=F,sep="\t",row.names=F)
}








