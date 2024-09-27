### 
### GSE135337
### 

## Load packages 
library(Seurat); library(dplyr); library(harmony); library(SeuratWrappers); 
library(patchwork); library(ggplot2)
## Load data
BC1 <- data.table::fread("./Data.GSE135337/GSM4006644_BC1_gene_cell_exprs_table.txt.gz") %>% as.data.frame()
BC2 <- data.table::fread("./Data.GSE135337/GSM4006645_BC2_gene_cell_exprs_table.txt.gz") %>% as.data.frame()
BC3 <- data.table::fread("./Data.GSE135337/GSM4006646_BC3_gene_cell_exprs_table.txt.gz") %>% as.data.frame()
BC4 <- data.table::fread("./Data.GSE135337/GSM4006647_BC4_gene_cell_exprs_table.txt.gz") %>% as.data.frame()
BC5 <- data.table::fread("./Data.GSE135337/GSM4006648_BC5_gene_cell_exprs_table.txt.gz") %>% as.data.frame()
BC6 <- data.table::fread("./Data.GSE135337/GSM4751267_BC6_gene_cell_exprs_table.txt.gz") %>% as.data.frame()
BC7 <- data.table::fread("./Data.GSE135337/GSM4751268_BC7_gene_cell_exprs_table.txt.gz") %>% as.data.frame()
BCN <- data.table::fread("./Data.GSE135337/GSM5329919_BCN_gene_cell_exprs_table.xls.gz") %>% as.data.frame()
#
BC1 <- aggregate(BC1[ ,3:ncol(BC1)], by=list(BC1$Symbol), sum); rownames(BC1) <- BC1$Group.1; BC1 <- BC1[ ,-1]
BC2 <- aggregate(BC2[ ,3:ncol(BC2)], by=list(BC2$Symbol), sum); rownames(BC2) <- BC2$Group.1; BC2 <- BC2[ ,-1]
BC3 <- aggregate(BC3[ ,3:ncol(BC3)], by=list(BC3$Symbol), sum); rownames(BC3) <- BC3$Group.1; BC3 <- BC3[ ,-1]
BC4 <- aggregate(BC4[ ,3:ncol(BC4)], by=list(BC4$Symbol), sum); rownames(BC4) <- BC4$Group.1; BC4 <- BC4[ ,-1]
BC5 <- aggregate(BC5[ ,3:ncol(BC5)], by=list(BC5$Symbol), sum); rownames(BC5) <- BC5$Group.1; BC5 <- BC5[ ,-1]
BC6 <- aggregate(BC6[ ,3:ncol(BC6)], by=list(BC6$Symbol), sum); rownames(BC6) <- BC6$Group.1; BC6 <- BC6[ ,-1]
BC7 <- aggregate(BC7[ ,3:ncol(BC7)], by=list(BC7$Symbol), sum); rownames(BC7) <- BC7$Group.1; BC7 <- BC7[ ,-1]
BCN <- aggregate(BCN[ ,3:ncol(BCN)], by=list(BCN$Symbol), sum); rownames(BCN) <- BCN$Group.1; BCN <- BCN[ ,-1]
## 
BC1st <- CreateSeuratObject(counts = BC1, project = "BC1", min.cells = 1, min.features = 200)
BC2st <- CreateSeuratObject(counts = BC2, project = "BC2", min.cells = 1, min.features = 200)
BC3st <- CreateSeuratObject(counts = BC3, project = "BC3", min.cells = 1, min.features = 200)
BC4st <- CreateSeuratObject(counts = BC4, project = "BC4", min.cells = 1, min.features = 200)
BC5st <- CreateSeuratObject(counts = BC5, project = "BC5", min.cells = 1, min.features = 200)
BC6st <- CreateSeuratObject(counts = BC6, project = "BC6", min.cells = 1, min.features = 200)
BC7st <- CreateSeuratObject(counts = BC7, project = "BC7", min.cells = 1, min.features = 200)
BCNst <- CreateSeuratObject(counts = BCN, project = "BCN", min.cells = 1, min.features = 200)
#
BC1st[["percent.mt"]] <- PercentageFeatureSet(BC1st, pattern = "^MT-")
BC2st[["percent.mt"]] <- PercentageFeatureSet(BC2st, pattern = "^MT-")
BC3st[["percent.mt"]] <- PercentageFeatureSet(BC3st, pattern = "^MT-")
BC4st[["percent.mt"]] <- PercentageFeatureSet(BC4st, pattern = "^MT-")
BC5st[["percent.mt"]] <- PercentageFeatureSet(BC5st, pattern = "^MT-")
BC6st[["percent.mt"]] <- PercentageFeatureSet(BC6st, pattern = "^MT-")
BC7st[["percent.mt"]] <- PercentageFeatureSet(BC7st, pattern = "^MT-")
BCNst[["percent.mt"]] <- PercentageFeatureSet(BCNst, pattern = "^MT-")
#
AddMetaData <- function(xx, Patient=NA, Gender=NA, Age=NA, Stage=NA, Types=NA, Source=NA, addID){
  xx[['Source' ]] <- Source  #  
  xx[['Patient']] <- Patient #  
  xx[['Gender' ]] <- Gender  #  
  xx[['Age'    ]] <- Age     #  
  xx[['Stage'  ]] <- Stage   #  
  xx[['Types'  ]] <- Types   #  
  xx <- RenameCells(xx, add.cell.id = addID) # 
  xx
}
BC1st <- AddMetaData(BC1st, Patient='BC1', addID='BC1', Source='SH', Gender='M',  Age=NA, Stage='pTa', Types='T')
BC2st <- AddMetaData(BC2st, Patient='BC2', addID='BC2', Source='SH', Gender='M',  Age=NA, Stage='pT1', Types='T')
BC3st <- AddMetaData(BC3st, Patient='BC3', addID='BC3', Source='SH', Gender='M',  Age=NA, Stage='pT1', Types='T')
BC4st <- AddMetaData(BC4st, Patient='BC4', addID='BC4', Source='SH', Gender='M',  Age=NA, Stage='pT2', Types='T')
BC5st <- AddMetaData(BC5st, Patient='BC5', addID='BC5', Source='SH', Gender='M',  Age=NA, Stage='pT3', Types='T')
BC6st <- AddMetaData(BC6st, Patient='BC6', addID='BC6', Source='SH', Gender='M',  Age=NA, Stage='pTa', Types='T')
BC7st <- AddMetaData(BC7st, Patient='BC7', addID='BC7', Source='SH', Gender='M',  Age=NA, Stage='pTa', Types='T')
BCNst <- AddMetaData(BCNst, Patient='BCN', addID='BCN', Source='SH', Gender='M',  Age=NA, Stage='pT2', Types='N')

## Normalizing and Dimensional Reduction
SimpleRunSeuratPipline <- function(xx, nfeatures=3000, ndims=1:50){
   require(Seurat); require(dplyr)
   xx <- NormalizeData(xx, normalization.method = "LogNormalize", scale.factor = 10000) 
   xx <- FindVariableFeatures(xx, selection.method = "vst", nfeatures = nfeatures)
   xx <- ScaleData(xx, features = rownames(xx))
   xx <- RunPCA(xx, features = VariableFeatures(object = xx), npcs=100)
   xx <- FindNeighbors(xx, dims = ndims) %>% FindClusters(resolution = c(1, 0.8, 0.6)) %>% RunUMAP(dims = ndims) %>% RunTSNE(dims = ndims)
   ## Cell Cycle analysis
   xx <- CellCycleScoring(xx, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = F)
   xx
}
BC1st <- SimpleRunSeuratPipline(BC1st, ndims=1:25)
BC2st <- SimpleRunSeuratPipline(BC2st, ndims=1:25)
BC3st <- SimpleRunSeuratPipline(BC3st, ndims=1:25)
BC4st <- SimpleRunSeuratPipline(BC4st, ndims=1:25)
BC5st <- SimpleRunSeuratPipline(BC5st, ndims=1:25)
BC6st <- SimpleRunSeuratPipline(BC6st, ndims=1:25)
BC7st <- SimpleRunSeuratPipline(BC7st, ndims=1:25)
BCNst <- SimpleRunSeuratPipline(BCNst, ndims=1:25)
save(BC1st,BC2st,BC3st,BC4st,BC5st,BC6st,BC7st,BCNst, file=file.path('RImage_GSE135337_01.RData'))


DimPlot(BC1st, label=T); VlnPlot(BC1st, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
DimPlot(BC2st, label=T); VlnPlot(BC2st, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
DimPlot(BC3st, label=T); VlnPlot(BC3st, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
DimPlot(BC4st, label=T); VlnPlot(BC4st, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
DimPlot(BC5st, label=T); VlnPlot(BC5st, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
DimPlot(BC6st, label=T); VlnPlot(BC6st, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
DimPlot(BC7st, label=T); VlnPlot(BC7st, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
DimPlot(BCNst, label=T); VlnPlot(BCNst, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 

# Markers 
MarkerJASN <- data.frame(Genes=c('EPCAM',"KRT18","KRT19","KRT5","KRT17","KRT13","KRT20","UPK1A","UPK1B","UPK3B","UPK3A","UPK2",
                                 "VIM","S100A4","COL3A1","COL1A1","COL1A2","ACTA2","TAGLN","DES","CNN1","ACTG2","TPM2",
								 "SELE","PECAM1","VCAM1","CDH5","LYZ","MS4A7","CD14","CD209","CD3D","CD3E","MZB1","CD79A","GPM6A"),
                       Markers=c('Epi','Epi','Epi','Bas','Bas','BasInter','Umbr','UmbrInter','UmbrInter','UmbrInter','Umbr','Umbr',
					             'Interstitial','Fibr','Fibr','Fibr','Fibr','myoFibr','Muscle','Muscle','Muscle','Muscle','Muscle',
								 'Endo','Endo','Endo','Endo','Mono','Mono','Mono','Dendritic','Tcell','Tcell','Bcell','Bcell','Neuron'))

### Merging Samples
BLCA_IJC <- merge(BC1st, list(BC2st,BC3st,BC4st,BC5st,BC6st,BC7st,BCNst))
BLCA_IJC <- SimpleRunSeuratPipline(BLCA_IJC, ndims=1:25)
### HMY intergration 
BLCAIJC_HMY <- RunHarmony(BLCA_IJC, group.by.vars = c("orig.ident"), assays='RNA')
BLCAIJC_HMY <- RunUMAP(BLCAIJC_HMY, reduction = "harmony", dims = 1:25)
BLCAIJC_HMY <- RunTSNE(BLCAIJC_HMY, reduction = "harmony", dims = 1:25)
BLCAIJC_HMY <- FindNeighbors(BLCAIJC_HMY, reduction = "harmony", dims = 1:25) %>% FindClusters(reduction = "harmony", resolution = c(1, 0.8, 0.6, 0.4))
#
DimPlot(BLCAIJC_HMY, label=T, group.by='seurat_clusters') 
DimPlot(BLCAIJC_HMY, label=T, reduction='tsne', group.by='seurat_clusters') 
DimPlot(BLCAIJC_HMY, label=T, reduction='tsne', group.by='RNA_snn_res.0.6') 
DimPlot(BLCAIJC_HMY, label=T, reduction='tsne', group.by='Types') 
DimPlot(BLCAIJC_HMY, label=T, reduction='tsne', group.by='orig.ident') 
DimPlot(BLCAIJC_HMY, label=T, reduction='tsne', group.by='Phase') 
#
DoHeatmap(MiceAll_HMY,features=MarkerJASN$Genes)
BLCAIJC_HMY@meta.data$CC1 <- dplyr::recode(BLCAIJC_HMY@meta.data$seurat_clusters, '0'='Basal', '1'='Basal','3'='Basal_Cycling', '4'='Basal','8'='Basal', '2'='Luminal',  '5'='Luminal', '7'='Luminal_Cyclling',
														   '9'='Monocyte',  '6'='CAFs', '10'='TCells','11'='BCells', '12'='Muscle','13'='Endo' )
BLCAIJC_HMY@meta.data$CC1 <- factor(BLCAIJC_HMY@meta.data$CC1, levels=c("Basal", "Basal_Cycling", "Luminal", "Luminal_Cyclling", "CAFs", "Muscle", "Monocyte", "TCells","BCells", "Endo"))
#
BLCAIJC_HMY$CC2 <- paste(as.character(BLCAIJC_HMY$Types), as.character(BLCAIJC_HMY$CC1), sep='_')
BLCAIJC_HMY$CC2 <- dplyr::recode(BLCAIJC_HMY$CC2, "T_TCells"="TCells", "N_TCells"="TCells", "T_Monocyte"='Monocyte', "N_Monocyte"="Monocyte",  "T_Endo"='Endo', "N_Endo"="Endo", 
                        "T_BCells"='BCells', "N_BCells"="BCells", "T_Muscle"="Muscle","N_Muscle"="Muscle", "T_CAFs"='Fibrobast',"N_CAFs"="Fibrobast" )
BLCAIJC_HMY$CC2 <- factor(BLCAIJC_HMY$CC2, levels=c("N_Basal", "N_Luminal", "N_Luminal_Cyclling", "T_Basal", "T_Basal_Cycling", "T_Luminal", "T_Luminal_Cyclling",
                        "Fibrobast", "Muscle", "TCells", "BCells","Monocyte", "Endo"))
DimPlot(BLCAIJC_HMY, label=T, reduction='tsne', group.by='CC2') 

save(BLCAIJC_HMY, file='RImage_GSE135337_03_HMY.RData')


### RunAlra 
BLCAIJC_Alra <- RunALRA(BLCAIJC_HMY)
save(BLCAIJC_Alra, file='RImage_GSE135337_02_Alra.RData')



## Plotting ##
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
##
WGCNA::sizeGrWindow(7,5)
tmpMarker1 <- CreatGGplotFromSeurat(BLCAIJC_HMY, assay='RNA', features=c('PIN1'))
tmpMarker2 <- reshape::melt(tmpMarker1[,c('CC1','PIN1')], id=c('CC1'))
LymCol2    <- c( ggsci::pal_npg( "nrc")(9), '#DFC297', '#FEB24C', '#ACD851')  
ggplot(tmpMarker2, aes(CC1, value, fill=CC1)) + 
   geom_violin(position=position_dodge(0.8), scale='width') + #facet_wrap(~variable, scales='free') + 
   geom_jitter(width=0.2, size=0.2) + 
   geom_boxplot(outlier.alpha=0, width=0.2, position=position_dodge(0.8), color='white') +
   theme_classic() + theme(legend.position = "none", validate = TRUE) + 
   theme(axis.text.x=element_text(angle=45, hjust=1, size=10, color='black')) +
   xlab('') + ylab('Expression Level (logUMI)') + scale_color_manual(values=LymCol2) + scale_fill_manual(values=LymCol2)  -> Fig01
ggsave(Fig01, file='./GSE135337.PIN1.RawExpr.pdf')
#
WGCNA::sizeGrWindow(7,5)
tmpMarker1 <- CreatGGplotFromSeurat(BLCAIJC_Alra, assay='alra', features=c('PIN1'))
tmpMarker2 <- reshape::melt(tmpMarker1[,c('CC1','PIN1')], id=c('CC1'))
LymCol2    <- c( ggsci::pal_npg( "nrc")(9), '#DFC297', '#FEB24C', '#ACD851')  
ggplot(tmpMarker2, aes(CC1, value, fill=CC1)) + 
   geom_violin(position=position_dodge(0.8), scale='width') + geom_jitter(width=0.2, size=0.2) + 
   geom_boxplot(outlier.alpha=0, width=0.2, position=position_dodge(0.8), color='white') +
   theme_classic() + theme(legend.position = "none", validate = TRUE) + 
   theme(axis.text.x=element_text(angle=45, hjust=1, size=10, color='black')) +
   xlab('') + ylab('Expression Level (logUMI)') + scale_color_manual(values=LymCol2) + scale_fill_manual(values=LymCol2) -> Fig02  
ggsave(Fig02, file='./GSE135337.PIN1.Alra.pdf')
#
WGCNA::sizeGrWindow(7,5)
tmpMarker1 <- CreatGGplotFromSeurat(BLCAIJC_Alra, assay='alra', features=c('PIN1'))
tmpMarker2 <- reshape::melt(tmpMarker1[,c('CC2','PIN1')], id=c('CC2'))
LymCol2    <- c( ggsci::pal_npg( "nrc")(10), '#DFC297', '#FEB24C', '#ACD851')  
ggplot(tmpMarker2, aes(CC2, value, fill=CC2)) + 
   geom_violin(position=position_dodge(0.8), scale='width') + geom_jitter(width=0.2, size=0.2) + 
   geom_boxplot(outlier.alpha=0, width=0.2, position=position_dodge(0.8), color='white') +
   theme_classic() + theme(legend.position = "none", validate = TRUE) + 
   theme(axis.text.x=element_text(angle=45, hjust=1, size=10, color='black')) +
   xlab('') + ylab('Expression Level (logUMI)') + scale_color_manual(values=LymCol2) + scale_fill_manual(values=LymCol2) -> Fig02  
Fig02
ggsave(Fig02, file='./GSE135337.PIN1.Alra2.pdf')


# 
## Plotting ##
WGCNA::sizeGrWindow(7,5)
#DimPlot(BLCAIJC_Alra, reduction='umap', label=T, group.by='CC1') #+ NoLegend() # -> Fig11
DimPlot(subset(BLCAIJC_Alra, Types=='N'), reduction='tsne', label=F, group.by='CC1' ) -> Fig11
DimPlot(subset(BLCAIJC_Alra, Types=='T'), reduction='tsne', label=F, group.by='CC1' ) -> Fig12
#
WGCNA::sizeGrWindow(7,6)
DefaultAssay(BLCAIJC_Alra) <- 'RNA'  ; FeaturePlot(subset(BLCAIJC_Alra, Types=='N'), reduction='tsne', order=T, pt.size=1, features='PIN1') + scale_color_viridis_c() -> Fig13
DefaultAssay(BLCAIJC_Alra) <- 'RNA'  ; FeaturePlot(subset(BLCAIJC_Alra, Types=='T'), reduction='tsne', order=T, pt.size=1, features='PIN1') + scale_color_viridis_c() -> Fig14
DefaultAssay(BLCAIJC_Alra) <- 'alra' ; FeaturePlot(subset(BLCAIJC_Alra, Types=='N'), reduction='tsne', order=T, pt.size=1, features='PIN1') + scale_color_viridis_c() -> Fig15
DefaultAssay(BLCAIJC_Alra) <- 'alra' ; FeaturePlot(subset(BLCAIJC_Alra, Types=='T'), reduction='tsne', order=T, pt.size=1, features='PIN1') + scale_color_viridis_c() -> Fig16
#
ggsave(Fig11, width=7, height=5, file='./GSE135337.Homo.PIN1.TSNE.Clusters.Nor.pdf')
ggsave(Fig12, width=7, height=5, file='./GSE135337.Homo.PIN1.TSNE.Clusters.Tum.pdf')
ggsave(Fig13, width=7, height=6, file='./GSE135337.Homo.PIN1.TSNE.RNA.Nor.pdf')
ggsave(Fig14, width=7, height=6, file='./GSE135337.Homo.PIN1.TSNE.RNA.Tum.pdf')
ggsave(Fig15, width=7, height=6, file='./GSE135337.Homo.PIN1.TSNE.Alra.Nor.pdf')
ggsave(Fig16, width=7, height=6, file='./GSE135337.Homo.PIN1.TSNE.Alra.Tum.pdf')



### Packages ### 
sessionInfo()
R version 4.2.2 (2022-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 22631)

Matrix products: default

locale:
[1] LC_COLLATE=Chinese (Simplified)_China.utf8  LC_CTYPE=Chinese (Simplified)_China.utf8    LC_MONETARY=Chinese (Simplified)_China.utf8
[4] LC_NUMERIC=C                                LC_TIME=Chinese (Simplified)_China.utf8    

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
 [71] ellipsis_0.3.2         ica_1.0-3              farver_2.1.1           pkgconfig_2.0.3        R.methodsS3_1.8.2     
 [76] nnet_7.3-19            uwot_0.1.16            deldir_1.0-9           utf8_1.2.3             dynamicTreeCut_1.63-1 
 [81] labeling_0.4.3         tidyselect_1.2.0       rlang_1.1.1            reshape2_1.4.4         later_1.3.1           
 [86] AnnotationDbi_1.60.2   munsell_0.5.0          tools_4.2.2            cachem_1.0.8           cli_3.6.1             
 [91] generics_0.1.3         RSQLite_2.3.1          ggridges_0.5.4         evaluate_0.21          stringr_1.5.0         
 [96] fastmap_1.1.1          ragg_1.2.5             goftest_1.2-3          knitr_1.44             bit64_4.0.5           
[101] fitdistrplus_1.1-11    purrr_1.0.2            RANN_2.6.1             KEGGREST_1.38.0        pbapply_1.7-2         
[106] future_1.33.0          nlme_3.1-163           mime_0.12              R.oo_1.25.0            rstudioapi_0.15.0     
[111] compiler_4.2.2         plotly_4.10.2          png_0.1-8              spatstat.utils_3.0-3   tibble_3.2.1          
[116] stringi_1.7.12         lattice_0.21-8         Matrix_1.6-1.1         vctrs_0.6.3            pillar_1.9.0          
[121] lifecycle_1.0.3        BiocManager_1.30.22    spatstat.geom_3.2-5    lmtest_0.9-40          RcppAnnoy_0.0.21      
[126] data.table_1.14.8      cowplot_1.1.1          bitops_1.0-7           irlba_2.3.5.1          httpuv_1.6.11         
[131] R6_2.5.1               promises_1.2.1         KernSmooth_2.23-22     gridExtra_2.3          IRanges_2.32.0        
[136] parallelly_1.36.0      codetools_0.2-19       MASS_7.3-60            withr_2.5.0            sctransform_0.4.0     
[141] GenomeInfoDbData_1.2.9 S4Vectors_0.36.2       parallel_4.2.2         rpart_4.1.19           grid_4.2.2            
[146] tidyr_1.3.0            rmarkdown_2.25         Rtsne_0.16             spatstat.explore_3.2-3 base64enc_0.1-3       
[151] Biobase_2.58.0         shiny_1.7.5            WGCNA_1.72-1          
