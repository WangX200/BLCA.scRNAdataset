###
### GSE129845
### 

## Load packages 
library(Seurat); library(dplyr); library(harmony); library(SeuratWrappers); 
library(patchwork); library(ggplot2)
## Load data
Sample1 <- Read10X('./Data.GSE129845/GSM3723357_Sample1')
Sample2 <- Read10X('./Data.GSE129845/GSM3723358_Sample2')
Sample3 <- Read10X('./Data.GSE129845/GSM3723359_Sample3')
Sample4 <- Read10X('./Data.GSE129845/GSM3723360_Sample4')  
Sample5 <- Read10X('./Data.GSE129845/GSM3723361_Sample5')  
#
Seurat1 <- CreateSeuratObject(counts = Sample1, project = "H1", min.cells = 1, min.features = 200)
Seurat2 <- CreateSeuratObject(counts = Sample2, project = "H2", min.cells = 1, min.features = 200)
Seurat3 <- CreateSeuratObject(counts = Sample3, project = "H3", min.cells = 1, min.features = 200)
Seurat4 <- CreateSeuratObject(counts = Sample4, project = "M1", min.cells = 1, min.features = 200)
Seurat5 <- CreateSeuratObject(counts = Sample5, project = "M2", min.cells = 1, min.features = 200)
#
Seurat1[["percent.mt"]] <- PercentageFeatureSet(Seurat1, pattern = "^MT-")
Seurat2[["percent.mt"]] <- PercentageFeatureSet(Seurat2, pattern = "^MT-")
Seurat3[["percent.mt"]] <- PercentageFeatureSet(Seurat3, pattern = "^MT-")
Seurat4[["percent.mt"]] <- PercentageFeatureSet(Seurat4, pattern = "^mt-")
Seurat5[["percent.mt"]] <- PercentageFeatureSet(Seurat5, pattern = "^mt-")
#
VlnPlot(Seurat1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
VlnPlot(Seurat2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
VlnPlot(Seurat3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
VlnPlot(Seurat4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
VlnPlot(Seurat5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 


## Normalizing and Dimensional Reduction
SimpleRunSeuratPipline <- function(xx, nfeatures=3000, ndims=1:50){
   require(Seurat); require(dplyr)
   xx <- NormalizeData(xx, normalization.method = "LogNormalize", scale.factor = 10000) 
   xx <- FindVariableFeatures(xx, selection.method = "vst", nfeatures = nfeatures)
   xx <- ScaleData(xx, features = rownames(xx))
   xx <- RunPCA(xx, features = VariableFeatures(object = xx), npcs=100)
   xx <- FindNeighbors(xx, dims = ndims) %>% FindClusters(resolution = c(1, 0.8, 0.6)) %>% RunUMAP(dims = ndims) %>% RunTSNE(dims = ndims)
   # Cell Cycle analysis
   xx <- CellCycleScoring(xx, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = F)
   xx
} 
Homo01 <- SimpleRunSeuratPipline(Seurat1, ndims=1:25)
Homo02 <- SimpleRunSeuratPipline(Seurat2, ndims=1:25)
Homo03 <- SimpleRunSeuratPipline(Seurat3, ndims=1:25)
Mice01 <- SimpleRunSeuratPipline(Seurat4, ndims=1:25)
Mice02 <- SimpleRunSeuratPipline(Seurat5, ndims=1:25)


## Identification Doublet Cells
DoubletFinderWrapper <- function(seu_obj, pcs=1:50, clusters='seurat_clusters', gt=FALSE){
   require(DoubletFinder) ; print('pK Identification')
   if (gt==FALSE){
   ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
     sweep.res.list <- paramSweep_v3(seu_obj, PCs = pcs, sct = FALSE)
     sweep.stats    <- summarizeSweep(sweep.res.list, GT = FALSE)
     bcmvn          <- find.pK(sweep.stats)
   }else{
   ## pK Identification (ground-truth) ------------------------------------------------------------------------------------------
     sweep.res.list <- paramSweep_v3(seu_obj, PCs = pcs, sct = FALSE)
     gt.calls    <- seu_obj@meta.data[rownames(sweep.res.list[[1]]), "GT"]  
     sweep.stats <- summarizeSweep(sweep.res.list, GT = TRUE, GT.calls = gt.calls)
     bcmvn       <- find.pK(sweep.stats)
   }
   print ('Homotypic Doublet Proportion Estimate')
   homotypic.prop <- modelHomotypic(seu_obj@meta.data[ ,clusters])  
   nExp_poi       <- round(0.075*nrow(seu_obj@meta.data))           
   nExp_poi.adj   <- round(nExp_poi*(1-homotypic.prop))
   
   print ('Run DoubletFinder with varying classification stringencies')
   seu_obj <- doubletFinder_v3(seu_obj, PCs = pcs, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
   tmp        <- grep('pANN', colnames(seu_obj@meta.data), value=T)
   seu_obj <- doubletFinder_v3(seu_obj, PCs = pcs, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = tmp, sct = FALSE)
   colnames(seu_obj@meta.data)[(ncol(seu_obj@meta.data)-2):ncol(seu_obj@meta.data)] <- c('pANN', 'DF1', 'DF2')
   seu_obj
}
#
Homo01 <- DoubletFinderWrapper(Homo01, pcs=1:25, clusters='seurat_clusters', gt=FALSE)
Homo02 <- DoubletFinderWrapper(Homo02, pcs=1:25, clusters='seurat_clusters', gt=FALSE)
Homo03 <- DoubletFinderWrapper(Homo03, pcs=1:25, clusters='seurat_clusters', gt=FALSE)
Mice01 <- DoubletFinderWrapper(Mice01, pcs=1:25, clusters='seurat_clusters', gt=FALSE)
Mice02 <- DoubletFinderWrapper(Mice02, pcs=1:25, clusters='seurat_clusters', gt=FALSE)
save(Homo01,Homo02,Homo03, Mice01,Mice02, file=file.path('RImage_GSE129845_01.RData'))

##
## QC
Homo01QC <- subset(Homo01, nFeature_RNA>=500 & nFeature_RNA<=4000 & percent.mt<=15 & DF2=='Singlet')
Homo02QC <- subset(Homo02, nFeature_RNA>=500 & nFeature_RNA<=4000 & percent.mt<=15 & DF2=='Singlet')
Homo03QC <- subset(Homo03, nFeature_RNA>=500 & nFeature_RNA<=4000 & percent.mt<=15 & DF2=='Singlet')
Mice01QC <- subset(Mice01, nFeature_RNA>=500 & nFeature_RNA<=4000 & percent.mt<=15 & DF2=='Singlet')
Mice02QC <- subset(Mice02, nFeature_RNA>=500 & nFeature_RNA<=4000 & percent.mt<=15 & DF2=='Singlet')
## Normalizing and Dimensional Reduction
Homo01QC <- SimpleRunSeuratPipline(Homo01QC, ndims=1:20)
Homo02QC <- SimpleRunSeuratPipline(Homo02QC, ndims=1:20)
Homo03QC <- SimpleRunSeuratPipline(Homo03QC, ndims=1:20)
Mice01QC <- SimpleRunSeuratPipline(Mice01QC, ndims=1:20)
Mice02QC <- SimpleRunSeuratPipline(Mice02QC, ndims=1:20)
#
DimPlot(Homo01QC, label=T) + VlnPlot(Homo01QC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
DimPlot(Homo02QC, label=T) + VlnPlot(Homo02QC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
DimPlot(Homo03QC, label=T) + VlnPlot(Homo03QC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
DimPlot(Mice01QC, label=T) + VlnPlot(Mice01QC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
DimPlot(Mice02QC, label=T) + VlnPlot(Mice02QC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
save(Homo01QC,Homo02QC,Homo03QC, Mice01QC,Mice02QC, file=file.path('RImage_GSE129845_02.RData'))

## Markers 
MarkerJASN <- data.frame(Genes=c("KRT18","KRT19","KRT5","KRT17","KRT13","KRT20","UPK1A","UPK1B","UPK3B","UPK3A","UPK2",
                                 "VIM","S100A4","COL3A1","COL1A1","COL1A2","ACTA2","TAGLN","DES","CNN1","ACTG2","TPM2",
								 "SELE","PECAM1","VCAM1","CDH5","LYZ","MS4A7","CD14","CD209","CD3D","CD3E","MZB1","CD79A","GPM6A"),
                       Markers=c('Epi','Epi','Bas','Bas','BasInter','Umbr','UmbrInter','UmbrInter','UmbrInter','Umbr','Umbr',
					             'Interstitial','Fibr','Fibr','Fibr','Fibr','myoFibr','Muscle','Muscle','Muscle','Muscle','Muscle',
								 'Endo','Endo','Endo','Endo','Mono','Mono','Mono','Dendritic','Tcell','Tcell','Bcell','Bcell','Neuron'))
DoHeatmap(Homo01QC,features=MarkerJASN$Genes)
DoHeatmap(Homo02QC,features=MarkerJASN$Genes)
DoHeatmap(Homo03QC,features=MarkerJASN$Genes)
DoHeatmap(Mice01QC,features=stringr::str_to_title(MarkerJASN$Genes))
DoHeatmap(Mice02QC,features=stringr::str_to_title(MarkerJASN$Genes))
## Assigning Cell Type
Homo01QC@meta.data$CC1 <- dplyr::recode(Homo01QC@meta.data$seurat_clusters, '0'='Fibr', '1'='TCells','2'='myoFibr', '3'='Endo', '4'='Mono','5'='Muscle','6'='Epi')
Homo02QC@meta.data$CC1 <- dplyr::recode(Homo02QC@meta.data$seurat_clusters, '0'='Fibr', '1'='Fibr','2'='myoFibr', '3'='Fibr', '4'='Bas','5'='Umbr','6'='Fibr','7'='Umbr', '8'='Fibr','9'='Mono', '10'='Mono','11'='TCells','12'='Endo')
Homo03QC@meta.data$CC1 <- dplyr::recode(Homo03QC@meta.data$seurat_clusters, '0'='Bas', '1'='Umbr','2'='UmbrInter', '3'='UmbrInter', '4'='UmbrInter',
                                            '5'='Muscle','6'='Fibr','7'='Fibr','8'='Endo','9'='TCells', '10'='Mono','11'='myoFibr','12'='Umbr','13'='BCell','14'='Fibr')
Mice01QC@meta.data$CC1 <- dplyr::recode(Mice01QC@meta.data$seurat_clusters, '0'='Bas', '1'='Umbr','2'='UmbrInter', '3'='Fibr', '4'='Bas','5'='Fibr','6'='UmbrInter','7'='myoFibr','8'='Umbr',
                                            '9'='Fibr', '10'='Mono','11'='unknown','12'='Neuron','13'='unknown','14'='Muscle','15'='Endo')
Mice02QC@meta.data$CC1 <- dplyr::recode(Mice02QC@meta.data$seurat_clusters, '0'='UmbrInter', '1'='UmbrInter','2'='Fibr', '3'='Bas', '4'='UmbrInter','5'='Bas',
                                            '6'='Fibr','7'='myoFibr','8'='Fibr','9'='Mono', '10'='Bas','11'='Fibr','12'='Umbr','13'='TCells','14'='Muscle','15'='Endo')


### 
### Merging Samples
HomoAll <- merge(Homo01QC, list(Homo02QC, Homo03QC))
MiceAll <- merge(Mice01QC, Mice02QC)
HomoAll <- SimpleRunSeuratPipline(HomoAll, ndims=1:20)
MiceAll <- SimpleRunSeuratPipline(MiceAll, ndims=1:20)
# 
DimPlot(HomoAll, label=T, group.by='orig.ident') 
DimPlot(HomoAll, label=T, group.by='CC1') 
DimPlot(HomoAll, label=T, group.by='seurat_clusters') 
VlnPlot(HomoAll, group.by='CC1', features=c('PIN1'))
# 
DimPlot(MiceAll, label=T, group.by='orig.ident') 
DimPlot(MiceAll, label=T, group.by='CC1') 
DimPlot(MiceAll, label=T, group.by='seurat_clusters') 
VlnPlot(MiceAll, group.by='CC1', features=c('Pin1'))


### HMY intergration 
{
   library(harmony); library(SeuratWrappers); library(patchwork)
   HomoAll_HMY <- RunHarmony(HomoAll, group.by.vars = c("orig.ident"), assays='RNA')
   HomoAll_HMY <- RunUMAP(HomoAll_HMY, reduction = "harmony", dims = 1:25)
   HomoAll_HMY <- RunTSNE(HomoAll_HMY, reduction = "harmony", dims = 1:25)
   HomoAll_HMY <- FindNeighbors(HomoAll_HMY, reduction = "harmony", dims = 1:25) %>% FindClusters(reduction = "harmony", resolution = c(1, 0.8, 0.6, 0.4))
 
   library(harmony); library(SeuratWrappers); library(patchwork)
   MiceAll_HMY <- RunHarmony(MiceAll, group.by.vars = c("orig.ident"), assays='RNA')
   MiceAll_HMY <- RunUMAP(MiceAll_HMY, reduction = "harmony", dims = 1:25)
   MiceAll_HMY <- RunTSNE(MiceAll_HMY, reduction = "harmony", dims = 1:25)
   MiceAll_HMY <- FindNeighbors(MiceAll_HMY, reduction = "harmony", dims = 1:25) %>% FindClusters(reduction = "harmony", resolution = c(1, 0.8, 0.6, 0.4))
}
  


### 
## Mouse Cell Clusters Annotation 
DoHeatmap(MiceAll_HMY,features=stringr::str_to_title(MarkerJASN$Genes))
DoHeatmap(MiceAll_HMY, group.by='RNA_snn_res.0.6',features=stringr::str_to_title(MarkerJASN$Genes))
DimPlot(MiceAll_HMY, reduction='tsne', label=T, group.by='RNA_snn_res.0.6') + 
DimPlot(MiceAll_HMY, reduction='tsne', label=T, group.by='CC1')
MiceAll_HMY@meta.data$CC2 <- dplyr::recode(MiceAll_HMY@meta.data$RNA_snn_res.0.6, '15'='TCells', '12'='Monocyte', '13'='DendriticCells', '14'='Neurone','17'='Endo','16'='Muscle', 
                                       '6'='Fibroblast','3'='Fibroblast', '9'='Fibroblast', '5'='myoFibroblast','11'='myoFibroblast', 
									   '4'='Basal','1'='Basal','2'='Basal','7'='Basal','0'='Luminal','8'='Luminal','10'='Luminal','18'='Luminal')
VlnPlot(MiceAll_HMY, group.by='CC2', features='Pin1')

## Human Cell Clusters Annotation 
DoHeatmap(HomoAll_HMY,features=stringr::str_to_title(MarkerJASN$Genes))
DoHeatmap(HomoAll_HMY, group.by='RNA_snn_res.0.6',features=stringr::str_to_title(MarkerJASN$Genes))
DimPlot(HomoAll_HMY, reduction='tsne', label=T   ) +
DimPlot(HomoAll_HMY, reduction='tsne', label=T, group.by='CC1')
HomoAll_HMY@meta.data$CC2 <- dplyr::recode(HomoAll_HMY@meta.data$seurat_clusters, '9'='Epi_TNNT1', '2'='Epi_Basal', '0'='Epi_Intermediate', '10'='Epi_Luminal','3'='Epi_Luminal','4'='Epi_Luminal', 
                                           '7'='Monocyte','11'='TCells',  '14'='Fibroblast','1'='Fibroblast','5'='myoFibroblast', '6'='Muscle','12'='Muscle', '8'='Endo', '13'='BCells')

##
MiceAll_HMY@meta.data$CC2 <- factor(MiceAll_HMY@meta.data$CC2, levels=c('Basal','Luminal','Bas.Lum.Mixed','Fibroblast','myoFibroblast','Muscle','Monocyte','DendriticCells','TCells','Endo','Neurone'))
HomoAll_HMY@meta.data$CC2 <- factor(HomoAll_HMY@meta.data$CC2, levels=c('Epi_Basal','Epi_Intermediate','Epi_Luminal','Epi_TNNT1','Fibroblast','myoFibroblast','Muscle','TCells','BCells','Monocyte','Endo'))
#
VlnPlot(MiceAll,     group.by='CC2', features='Pin1')
VlnPlot(HomoAll_HMY, group.by='CC2', features='PIN1')
FeaturePlot(HomoAll_HMY, reduction='tsne', feature='TNNT1')
DimPlot(MiceAll_HMY, reduction='tsne', label=T, group.by='CC2')
DimPlot(HomoAll_HMY, reduction='tsne', label=T, group.by='CC2')
save(HomoAll_HMY, MiceAll_HMY, file='RImage_GSE129845_03.RData')

## RunALRA
MiceAllAlra <- RunALRA(MiceAll_HMY)
HomoAllAlra <- RunALRA(HomoAll_HMY)
save(MiceAllAlra, HomoAllAlra, file='RImage_GSE129845_04.RData')

## 
library(ggsci); library(ggplot2)
CreatGGplotFromSeurat <- function(xx, assay=NA, features=NA, meta=NA){
   require(dplyr); DefaultAssay(xx) <- assay
   features <- features[features %in% rownames(xx)]
   if (length(features)<1){ stop('None of the features include in the assay') }
   print(features)
   tmpAssay <- GetAssay(xx[features, ], assay=assay)@data
   tmpAssay <- tmpAssay %>% as.data.frame() %>% t()
   if (is.na(meta)){meta = colnames(xx@meta.data)}
   x <- cbind(xx@meta.data[,meta], tmpAssay)
   colnames(x) <- gsub('-', '_', colnames(x))
   x
}

WGCNA::sizeGrWindow(7,5)
tmpMarker1 <- CreatGGplotFromSeurat(MiceAllAlra, assay='RNA', features=c('Pin1'))
tmpMarker2 <- reshape::melt(tmpMarker1[,c('CC2','Pin1')], id=c('CC2'))
LymCol2    <- c( pal_npg( "nrc")(9), '#DFC297', '#FEB24C', '#ACD851')  
ggplot(tmpMarker2, aes(CC2, value, fill=CC2)) + 
   geom_violin(position=position_dodge(0.8), scale='width') + #facet_wrap(~variable, scales='free') + 
   geom_jitter(width=0.2, size=0.2) + 
   geom_boxplot(outlier.alpha=0, width=0.2, position=position_dodge(0.8), color='white') +
   theme_classic() + theme(legend.position = "none", validate = TRUE) + 
   theme(axis.text.x=element_text(angle=45, hjust=1, size=10, color='black')) +
   xlab('') + ylab('Expression Level (logUMI)') + scale_color_manual(values=LymCol2) + scale_fill_manual(values=LymCol2)  -> Fig03
ggsave(Fig03, file='./GSE129845.Mice.PIN1.RawExpr.pdf')

WGCNA::sizeGrWindow(7,5)
tmpMarker1 <- CreatGGplotFromSeurat(MiceAllAlra, assay='alra', features=c('Pin1'))
tmpMarker2 <- reshape::melt(tmpMarker1[,c('CC2','Pin1')], id=c('CC2'))
LymCol2    <- c( pal_npg( "nrc")(9), '#DFC297', '#FEB24C', '#ACD851')  
ggplot(tmpMarker2, aes(CC2, value, fill=CC2)) + 
   geom_violin(position=position_dodge(0.8), scale='width') + #facet_wrap(~variable, scales='free') + 
   geom_jitter(width=0.2, size=0.2) + 
   geom_boxplot(outlier.alpha=0, width=0.2, position=position_dodge(0.8), color='white') +
   theme_classic() + theme(legend.position = "none", validate = TRUE) + 
   theme(axis.text.x=element_text(angle=45, hjust=1, size=10, color='black')) +
   xlab('') + ylab('Expression Level (logUMI)') + scale_color_manual(values=LymCol2) + scale_fill_manual(values=LymCol2)  -> Fig04
ggsave(Fig04, file='./GSE129845.Mice.PIN1.Alra.pdf')


WGCNA::sizeGrWindow(7,5)
tmpMarker1 <- CreatGGplotFromSeurat(HomoAllAlra, assay='RNA', features=c('PIN1'))
tmpMarker2 <- reshape::melt(tmpMarker1[,c('CC2','PIN1')], id=c('CC2'))
LymCol2    <- c( pal_npg( "nrc")(9), '#DFC297', '#FEB24C', '#ACD851')  
ggplot(tmpMarker2, aes(CC2, value, fill=CC2)) + 
   geom_violin(position=position_dodge(0.8), scale='width') + #facet_wrap(~variable, scales='free') + 
   geom_jitter(width=0.2, size=0.2) + 
   geom_boxplot(outlier.alpha=0, width=0.2, position=position_dodge(0.8), color='white') +
   theme_classic() + theme(legend.position = "none", validate = TRUE) + 
   theme(axis.text.x=element_text(angle=45, hjust=1, size=10, color='black')) +
   xlab('') + ylab('Expression Level (logUMI)') + scale_color_manual(values=LymCol2) + scale_fill_manual(values=LymCol2)  -> Fig05
ggsave(Fig05, file='./GSE129845.Homo.PIN1.RawExpr.pdf')

WGCNA::sizeGrWindow(7,5)
tmpMarker1 <- CreatGGplotFromSeurat(HomoAllAlra, assay='alra', features=c('PIN1'))
tmpMarker2 <- reshape::melt(tmpMarker1[,c('CC2','PIN1')], id=c('CC2'))
LymCol2    <- c( pal_npg( "nrc")(9), '#DFC297', '#FEB24C', '#ACD851')  
ggplot(tmpMarker2, aes(CC2, value, fill=CC2)) + 
   geom_violin(position=position_dodge(0.8), scale='width') + #facet_wrap(~variable, scales='free') + 
   geom_jitter(width=0.2, size=0.2) + 
   geom_boxplot(outlier.alpha=0, width=0.2, position=position_dodge(0.8), color='white') +
   theme_classic() + theme(legend.position = "none", validate = TRUE) + 
   theme(axis.text.x=element_text(angle=45, hjust=1, size=10, color='black')) +
   xlab('') + ylab('Expression Level (logUMI)') + scale_color_manual(values=LymCol2) + scale_fill_manual(values=LymCol2)  -> Fig06
ggsave(Fig06, file='./GSE129845.Homo.PIN1.Alra.pdf')


##
WGCNA::sizeGrWindow(7,5)
DimPlot(MiceAllAlra, reduction='tsne', label=F, group.by='CC2')  -> Fig11
DimPlot(HomoAllAlra, reduction='tsne', label=F, group.by='CC2')  -> Fig21
#
DefaultAssay(MiceAllAlra) <- 'RNA'  ; FeaturePlot(MiceAllAlra, reduction='tsne', order=T, pt.size=1, features='Pin1') + scale_color_viridis_c() -> Fig12
DefaultAssay(MiceAllAlra) <- 'alra' ; FeaturePlot(MiceAllAlra, reduction='tsne', order=T, pt.size=1, features='Pin1') + scale_color_viridis_c() -> Fig13

DefaultAssay(HomoAllAlra) <- 'RNA'  ; FeaturePlot(HomoAllAlra, reduction='tsne', order=T, pt.size=1, features='PIN1') + scale_color_viridis_c() -> Fig22
DefaultAssay(HomoAllAlra) <- 'alra' ; FeaturePlot(HomoAllAlra, reduction='tsne', order=T, pt.size=1, features='PIN1') + scale_color_viridis_c() -> Fig23
# 
ggsave(Fig11, width=7, height=5, file='./GSE129845.Mice.PIN1.TSNE.Clusters.pdf')
ggsave(Fig12, width=7, height=6, file='./GSE129845.Mice.PIN1.TSNE.RNA.pdf')
ggsave(Fig13, width=7, height=6, file='./GSE129845.Mice.PIN1.TSNE.Alra.pdf')
ggsave(Fig21, width=7, height=5, file='./GSE129845.Homo.PIN1.TSNE.Clusters.pdf')
ggsave(Fig22, width=7, height=6, file='./GSE129845.Homo.PIN1.TSNE.RNA.pdf')
ggsave(Fig23, width=7, height=6, file='./GSE129845.Homo.PIN1.TSNE.Alra.pdf')



### Packages ### 
sessionInfo()
R version 4.2.2 (2022-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 22631)

Matrix products: default

locale:
[1] LC_COLLATE=Chinese (Simplified)_China.utf8  LC_CTYPE=Chinese (Simplified)_China.utf8   
[3] LC_MONETARY=Chinese (Simplified)_China.utf8 LC_NUMERIC=C                               
[5] LC_TIME=Chinese (Simplified)_China.utf8    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] ggplot2_3.4.3        patchwork_1.1.3      SeuratWrappers_0.3.1 harmony_1.0.1        Rcpp_1.0.11          dplyr_1.1.3         
[7] SeuratObject_4.1.3   Seurat_4.3.0.1      

loaded via a namespace (and not attached):
  [1] backports_1.4.1        Hmisc_5.1-1            systemfonts_1.0.4      plyr_1.8.8             igraph_1.5.1          
  [6] lazyeval_0.2.2         sp_2.0-0               splines_4.2.2          listenv_0.9.0          scattermore_1.2       
 [11] GenomeInfoDb_1.34.9    digest_0.6.33          foreach_1.5.2          htmltools_0.5.6        GO.db_3.16.0          
 [16] fansi_1.0.4            checkmate_2.2.0        magrittr_2.0.3         memoise_2.0.1          tensor_1.5            
 [21] cluster_2.1.4          doParallel_1.0.17      ROCR_1.0-11            remotes_2.4.2.1        globals_0.16.2        
 [26] fastcluster_1.2.3      Biostrings_2.66.0      matrixStats_1.0.0      R.utils_2.12.2         spatstat.sparse_3.0-2 
 [31] colorspace_2.1-0       blob_1.2.4             ggrepel_0.9.3          textshaping_0.3.6      xfun_0.40             
 [36] crayon_1.5.2           RCurl_1.98-1.12        jsonlite_1.8.7         progressr_0.14.0       spatstat.data_3.0-1   
 [41] impute_1.72.3          survival_3.5-7         zoo_1.8-12             iterators_1.0.14       glue_1.6.2            
 [46] polyclip_1.10-4        gtable_0.3.4           zlibbioc_1.44.0        XVector_0.38.0         leiden_0.4.3          
 [51] future.apply_1.11.0    BiocGenerics_0.44.0    abind_1.4-5            scales_1.2.1           DBI_1.1.3             
 [56] spatstat.random_3.1-6  miniUI_0.1.1.1         htmlTable_2.4.1        viridisLite_0.4.2      xtable_1.8-4          
 [61] reticulate_1.32.0      foreign_0.8-85         bit_4.0.5              rsvd_1.0.5             preprocessCore_1.60.2 
 [66] Formula_1.2-5          stats4_4.2.2           htmlwidgets_1.6.2      httr_1.4.7             RColorBrewer_1.1-3    
 [71] ellipsis_0.3.2         ica_1.0-3              farver_2.1.1           reshape_0.8.9          pkgconfig_2.0.3       
 [76] R.methodsS3_1.8.2      nnet_7.3-19            uwot_0.1.16            deldir_1.0-9           utf8_1.2.3            
 [81] dynamicTreeCut_1.63-1  labeling_0.4.3         tidyselect_1.2.0       rlang_1.1.1            reshape2_1.4.4        
 [86] later_1.3.1            AnnotationDbi_1.60.2   munsell_0.5.0          tools_4.2.2            cachem_1.0.8          
 [91] cli_3.6.1              generics_0.1.3         RSQLite_2.3.1          ggridges_0.5.4         evaluate_0.21         
 [96] stringr_1.5.0          fastmap_1.1.1          ragg_1.2.5             goftest_1.2-3          knitr_1.44            
[101] bit64_4.0.5            fitdistrplus_1.1-11    purrr_1.0.2            RANN_2.6.1             KEGGREST_1.38.0       
[106] pbapply_1.7-2          future_1.33.0          nlme_3.1-163           mime_0.12              R.oo_1.25.0           
[111] rstudioapi_0.15.0      compiler_4.2.2         plotly_4.10.2          png_0.1-8              spatstat.utils_3.0-3  
[116] tibble_3.2.1           stringi_1.7.12         lattice_0.21-8         Matrix_1.6-1.1         vctrs_0.6.3           
[121] pillar_1.9.0           lifecycle_1.0.3        BiocManager_1.30.22    spatstat.geom_3.2-5    lmtest_0.9-40         
[126] RcppAnnoy_0.0.21       data.table_1.14.8      cowplot_1.1.1          bitops_1.0-7           irlba_2.3.5.1         
[131] httpuv_1.6.11          R6_2.5.1               promises_1.2.1         KernSmooth_2.23-22     gridExtra_2.3         
[136] IRanges_2.32.0         parallelly_1.36.0      codetools_0.2-19       MASS_7.3-60            withr_2.5.0           
[141] sctransform_0.4.0      GenomeInfoDbData_1.2.9 S4Vectors_0.36.2       parallel_4.2.2         rpart_4.1.19          
[146] grid_4.2.2             tidyr_1.3.0            rmarkdown_2.25         Rtsne_0.16             spatstat.explore_3.2-3
[151] base64enc_0.1-3        Biobase_2.58.0         shiny_1.7.5            WGCNA_1.72-1          
