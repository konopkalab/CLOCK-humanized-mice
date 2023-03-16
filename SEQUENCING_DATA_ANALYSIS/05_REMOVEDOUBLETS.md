# REMOVE DOUBLETS USING DOUBLETFINDER

## P07 WT
- python/3.7.x-anaconda
- gcc/8.3.0
- hdf5_18/1.8.17
- R/4.0.2-gccmkl
```{R}
##------------------------------------
## SEURAT CLUSTERING
##------------------------------------
## Libraries
rm(list = ls())
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(reshape2)
library(clustree)
set.seed(10)

##------------------------------------
## reference gene symbol list
ref <- read.table("gencode.vM17.protein_coding_gene_id_symbol.txt.gz", header = TRUE, sep = "\t")

##------------------------------------
## Read CellBender H5
p07wt02 <- Read10X_h5("CELLBENDER_OUTPUT/P07_WT/YL_P07_WT_02_CB/CellBender_Out_filtered.h5") # 21967  6990
colnames(p07wt02) <- paste("YL_P07_WT_02", colnames(p07wt02), sep = "_")

p07wt04 <- Read10X_h5("CELLBENDER_OUTPUT/P07_WT/YL_P07_WT_04_CB/CellBender_Out_filtered.h5") # 21967     8
colnames(p07wt04) <- paste("YL_P07_WT_04", colnames(p07wt04), sep = "_")

p07wt11 <- Read10X_h5("CELLBENDER_OUTPUT/P07_WT/YL_P07_WT_11_CB/CellBender_Out_filtered.h5") # 21967 11406
colnames(p07wt11) <- paste("YL_P07_WT_11", colnames(p07wt11), sep = "_")

p07wt16 <- Read10X_h5("CELLBENDER_OUTPUT/P07_WT/YL_P07_WT_16_CB/CellBender_Out_filtered.h5") # 21967  6405
colnames(p07wt16) <- paste("YL_P07_WT_16", colnames(p07wt16), sep = "_")

##------------------------------------
## fetch genes or rows corresponding to gencode protein coding gene symbols
p07wt02.temp <- p07wt02[row.names(p07wt02) %in% ref$GeneSymbol,]
p07wt04.temp <- p07wt04[row.names(p07wt04) %in% ref$GeneSymbol,]
p07wt11.temp <- p07wt11[row.names(p07wt11) %in% ref$GeneSymbol,]
p07wt16.temp <- p07wt16[row.names(p07wt16) %in% ref$GeneSymbol,]

##------------------------------------
## make a data from from matrix
YL_P07_WT_02 <- as.data.frame(as.matrix(p07wt02.temp))
YL_P07_WT_04 <- as.data.frame(as.matrix(p07wt04.temp))
YL_P07_WT_11 <- as.data.frame(as.matrix(p07wt11.temp))
YL_P07_WT_16 <- as.data.frame(as.matrix(p07wt16.temp))

##------------------------------------
## add gene symbol as a column
YL_P07_WT_02$Genes <- row.names(YL_P07_WT_02)
YL_P07_WT_04$Genes <- row.names(YL_P07_WT_04)
YL_P07_WT_11$Genes <- row.names(YL_P07_WT_11)
YL_P07_WT_16$Genes <- row.names(YL_P07_WT_16)

##------------------------------------
## combine individual tables into a giant data frame
## Sample "YL_P07_WT_04" is dropped being an outlier
dataCombined <- list("YL_P07_WT_02" = YL_P07_WT_02, "YL_P07_WT_11" = YL_P07_WT_11, "YL_P07_WT_16" = YL_P07_WT_16)

combinedData <- Reduce(function(x, y) { merge(x, y, all = TRUE, by = "Genes") } , dataCombined)
row.names(combinedData) <- combinedData$Genes
combinedData$Genes <- NULL

yl.data <- combinedData[row.names(combinedData) %in% ref$GeneSymbol,]

yl.data[is.na(yl.data)] <- 0
print(dim(yl.data)) ## 21967 45970

##------------------------------------
## Save Count Table
p07.wt.data <- yl.data
save(p07.wt.data, file = "YL_P07_WT_CB_COUNTS.RData")

##------------------------------------
## Initialize the Seurat object with the raw (non-normalized data).
seuObj <- CreateSeuratObject(counts = p07.wt.data, project = "YL_P07_WT_CB")
print(dim(seuObj)) ## 21967 24801

##------------------------------------
## Add and update meta data
metaTemp <- as.data.frame(seuObj@meta.data)
metaTemp2 <- as.data.frame(matrix(unlist(strsplit(row.names(metaTemp), "_")), ncol = 5, byrow = TRUE))
row.names(metaTemp2) <- row.names(metaTemp)
metaData <- merge(metaTemp, metaTemp2, by = "row.names")
row.names(metaData) <- metaData$Row.names
metaData$Row.names <- NULL
colnames(metaData) <- c("Ident", "nUMI", "nGenes", "Project", "Age", "Genotype", "Sample", "CellBarcode")

metaSample <- metaData$Sample
names(metaSample) <- row.names(metaData)

metaAge <- metaData$Age
names(metaAge) <- row.names(metaData)

metaGeno <- metaData$Genotype
names(metaGeno) <- row.names(metaData)

seuObj$Genotype <- metaGeno
seuObj$Sample <- metaSample
seuObj$Age <- metaAge
seuObj$GenoAgeSample <- paste(seuObj$Genotype, seuObj$Age, seuObj$Sample, sep = "_")

##------------------------------------
## Calculate Percent Mito
seuObj[["pMito_RNA"]] <- PercentageFeatureSet(seuObj, pattern = "^mt-")

## Set default identities to GenoAgeSample
Idents(seuObj) <- "GenoAgeSample"

table(seuObj@active.ident)
# WT_P07_02 WT_P07_11 WT_P07_16 
#      6990     11406      6405

##------------------------------------
## Visualize Data QC
p9 <- VlnPlot(seuObj, features = c("nFeature_RNA", "nCount_RNA", "pMito_RNA"), ncol = 3, pt.size = 0)
p2 <- FeatureScatter(seuObj, feature1 = "nCount_RNA", feature2 = "pMito_RNA")
p3 <- FeatureScatter(seuObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

ggsave(filename = "SEURAT_YL_QC_1.pdf", plot = p9, width = 8, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_YL_QC_2.pdf", plot = p2, width = 5, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_YL_QC_3.pdf", plot = p3, width = 5, height = 4, units = "in", dpi = 150)

##------------------------------------
## Data Normalization
seuObj <- NormalizeData(seuObj)

##------------------------------------
## Identify top 2000 highly variable genes
seuObj <- FindVariableFeatures(seuObj, selection.method = "vst", nfeatures = 2000)

## Identify the 10 most highly variable genes
top90 <- head(VariableFeatures(seuObj), 10)

## plot variable features with and without labels
p4 <- VariableFeaturePlot(seuObj)
p5 <- LabelPoints(plot = p4, points = top90, repel = TRUE)

ggsave(filename = "SEURAT_YL_VARGENES_1.pdf", plot = p4, width = 8, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_YL_VARGENES_2.pdf", plot = p5, width = 8, height = 4, units = "in", dpi = 150)

##------------------------------------
## Data Scaling
seuObj <- ScaleData(seuObj, vars.to.regress = c("nCount_RNA"))

##------------------------------------
## PCA & Jackstraw
## Compute PCA
seuObj <- RunPCA(seuObj, features = VariableFeatures(object = seuObj), npcs = 100)

## Compute and Score Jackstraw
seuObj <- JackStraw(seuObj, num.replicate = 100, dims = 100)
seuObj <- ScoreJackStraw(seuObj, dims = 1:100)

## Plot PCA loadings and Jackstraw scores
p6 <- ElbowPlot(seuObj, ndims = 100)
p7 <- JackStrawPlot(seuObj, dims = 1:100)

ggsave(filename = "SEURAT_YL_PCA_1.pdf", plot = p6, width = 6, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_YL_PCA_2.pdf", plot = p7, width = 12, height = 6, units = "in", dpi = 150)

##------------------------------------
## Data Clustering
## IDENTIFY PCS
pcScores <- seuObj@reductions$pca@jackstraw$overall.p.values
pcScoresNoSign <- pcScores[pcScores[,2] > 0.05,]
selpcs <- min(pcScoresNoSign[,1]) - 1
print(selpcs)
# 53

seuObj <- FindNeighbors(seuObj, dims = 1:selpcs)
seuObj <- FindClusters(seuObj, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 0.8, 1.4, 1.6, 1.8, 2.0))
seuObj <- RunUMAP(seuObj, dims = 1:selpcs)

##------------------------------------
## Save RData
save(seuObj, file = "YL_P07_WT_CB_CLUST.RData")


##------------------------------------
## DOUBLET FINDER
##------------------------------------
## Libraries
rm(list = ls())
library(dplyr)
library(patchwork)
library(ggplot2)
library(reshape2)
library(clustree)
library(Seurat)
library(DoubletFinder)
set.seed(10)

load("YL_P07_WT_CB_CLUST.RData")
# seuObj

pcScores <- seuObj@reductions$pca@jackstraw$overall.p.values
pcScoresNoSign <- pcScores[pcScores[,2] > 0.05,]
selpcs <- min(pcScoresNoSign[,1]) - 1
print(selpcs)
# 53

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_str <- paramSweep_v3(seuObj, PCs = 1:selpcs, sct = FALSE)
sweep.stats_str <- summarizeSweep(sweep.res.list_str, GT = FALSE)
bcmvn_str <- find.pK(sweep.stats_str)
selpk <- 0.005

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(seuObj@meta.data$RNA_snn_res.0.8) ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.1*nrow(seuObj@meta.data))  ## Assuming 10% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seuObjDF <- doubletFinder_v3(seuObj, PCs = 1:selpcs, pN = 0.25, pK = selpk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

seuObjDF <- doubletFinder_v3(seuObjDF, PCs = 1:selpcs, pN = 0.25, pK = selpk, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_2480", sct = FALSE)

table(seuObjDF$GenoAgeSample, seuObjDF$DF.classifications_0.25_0.005_2480)
#             Doublet Singlet
#   WT_P07_02     541    6449
#   WT_P07_11    1264   10142
#   WT_P07_16     675    5730

table(seuObjDF$GenoAgeSample, seuObjDF$DF.classifications_0.25_0.005_2308)
#             Doublet Singlet
#   WT_P07_02     498    6492
#   WT_P07_11    1199   10207
#   WT_P07_16     611    5794

seuMeta <- as.data.frame(seuObjDF@meta.data)
seuMetaFiltTemp <- seuMeta[seuMeta$DF.classifications_0.25_0.005_2480 == "Singlet",]
seuMetaFilt <- seuMetaFiltTemp[seuMetaFiltTemp$DF.classifications_0.25_0.005_2308 == "Singlet",]

## Before Doublets Removal
table(seuMeta$GenoAgeSample)
# WT_P07_02 WT_P07_11 WT_P07_16 
#      6990     11406      6405

## After Doublets Removal
table(seuMetaFilt$GenoAgeSample)
# WT_P07_02 WT_P07_11 WT_P07_16 
#      6449     10142      5730

seuMetaFilt$DoubletFinder <- rep("Retained", nrow(seuMetaFilt))

df.meta <- seuMetaFilt$DoubletFinder
names(df.meta) <- row.names(seuMetaFilt)

seuObjDF$DoubletFinder <- df.meta
seuObjDF$DoubletFinder[is.na(seuObjDF$DoubletFinder)] <- "Discarded"

##------------------------------------
## Save RData
save(seuObjDF, file = "YL_P07_WT_CB_DF.RData")
save(seuMeta, seuMetaFilt, file = "YL_P07_WT_CB_DF_FILT_META.RData")
```
<br></br>


## P07 KO
- python/3.7.x-anaconda
- gcc/8.3.0
- hdf5_18/1.8.17
- R/4.0.2-gccmkl
```{R}
##------------------------------------
## SEURAT CLUSTERING
##------------------------------------
## Libraries
rm(list = ls())
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(reshape2)
library(clustree)
set.seed(10)

##------------------------------------
## reference gene symbol list
ref <- read.table("gencode.vM17.protein_coding_gene_id_symbol.txt.gz", header = TRUE, sep = "\t")

##------------------------------------
## Read CellBender H5
p07ko01 <- Read10X_h5("CELLBENDER_OUTPUT/P07_KO/YL_P07_KO_01_CB/CellBender_Out_filtered.h5") # 21967  9512
colnames(p07ko01) <- paste("YL_P07_KO_01", colnames(p07ko01), sep = "_")

p07ko06 <- Read10X_h5("CELLBENDER_OUTPUT/P07_KO/YL_P07_KO_06_CB/CellBender_Out_filtered.h5") # 21967  16
colnames(p07ko06) <- paste("YL_P07_KO_06", colnames(p07ko06), sep = "_")

p07ko09 <- Read10X_h5("CELLBENDER_OUTPUT/P07_KO/YL_P07_KO_09_CB/CellBender_Out_filtered.h5") # 21967  10599
colnames(p07ko09) <- paste("YL_P07_KO_09", colnames(p07ko09), sep = "_")

p07ko18 <- Read10X_h5("CELLBENDER_OUTPUT/P07_KO/YL_P07_KO_18_CB/CellBender_Out_filtered.h5") # 21967  7934
colnames(p07ko18) <- paste("YL_P07_KO_18", colnames(p07ko18), sep = "_")

##------------------------------------
## fetch genes or rows corresponding to gencode protein coding gene symbols
p07ko01.temp <- p07ko01[row.names(p07ko01) %in% ref$GeneSymbol,]
p07ko06.temp <- p07ko06[row.names(p07ko06) %in% ref$GeneSymbol,]
p07ko09.temp <- p07ko09[row.names(p07ko09) %in% ref$GeneSymbol,]
p07ko18.temp <- p07ko18[row.names(p07ko18) %in% ref$GeneSymbol,]

##------------------------------------
## make a data from from matrix
YL_P07_KO_01 <- as.data.frame(as.matrix(p07ko01.temp))
YL_P07_KO_06 <- as.data.frame(as.matrix(p07ko06.temp))
YL_P07_KO_09 <- as.data.frame(as.matrix(p07ko09.temp))
YL_P07_KO_18 <- as.data.frame(as.matrix(p07ko18.temp))

##------------------------------------
## add gene symbol as a column
YL_P07_KO_01$Genes <- row.names(YL_P07_KO_01)
YL_P07_KO_06$Genes <- row.names(YL_P07_KO_06)
YL_P07_KO_09$Genes <- row.names(YL_P07_KO_09)
YL_P07_KO_18$Genes <- row.names(YL_P07_KO_18)

##------------------------------------
## combine individual tables into a giant data frame
## Sample "YL_P07_KO_06" is dropped being an outlier 
dataCombined <- list("YL_P07_KO_01" = YL_P07_KO_01, "YL_P07_KO_09" = YL_P07_KO_09, "YL_P07_KO_18" = YL_P07_KO_18)

combinedData <- Reduce(function(x, y) { merge(x, y, all = TRUE, by = "Genes") } , dataCombined)
row.names(combinedData) <- combinedData$Genes
combinedData$Genes <- NULL

yl.data <- combinedData[row.names(combinedData) %in% ref$GeneSymbol,]

yl.data[is.na(yl.data)] <- 0
print(dim(yl.data)) ## 21967 45970

##------------------------------------
## Save Count Table
p07.ko.data <- yl.data
save(p07.ko.data, file = "YL_P07_KO_CB_COUNTS.RData")

##------------------------------------
## Initialize the Seurat object with the raw (non-normalized data).
seuObj <- CreateSeuratObject(counts = p07.ko.data, project = "YL_P07_KO_CB")
print(dim(seuObj)) ## 21967 24801

##------------------------------------
## Add and update meta data
metaTemp <- as.data.frame(seuObj@meta.data)
metaTemp2 <- as.data.frame(matrix(unlist(strsplit(row.names(metaTemp), "_")), ncol = 5, byrow = TRUE))
row.names(metaTemp2) <- row.names(metaTemp)
metaData <- merge(metaTemp, metaTemp2, by = "row.names")
row.names(metaData) <- metaData$Row.names
metaData$Row.names <- NULL
colnames(metaData) <- c("Ident", "nUMI", "nGenes", "Project", "Age", "Genotype", "Sample", "CellBarcode")

metaSample <- metaData$Sample
names(metaSample) <- row.names(metaData)

metaAge <- metaData$Age
names(metaAge) <- row.names(metaData)

metaGeno <- metaData$Genotype
names(metaGeno) <- row.names(metaData)

seuObj$Genotype <- metaGeno
seuObj$Sample <- metaSample
seuObj$Age <- metaAge
seuObj$GenoAgeSample <- paste(seuObj$Genotype, seuObj$Age, seuObj$Sample, sep = "_")

##------------------------------------
## Calculate Percent Mito
seuObj[["pMito_RNA"]] <- PercentageFeatureSet(seuObj, pattern = "^mt-")

## Set default identities to GenoAgeSample
Idents(seuObj) <- "GenoAgeSample"

table(seuObj@active.ident)
# KO_P07_01 KO_P07_09 KO_P07_18 
#      9512     10599      7934

##------------------------------------
## Visualize Data QC
p9 <- VlnPlot(seuObj, features = c("nFeature_RNA", "nCount_RNA", "pMito_RNA"), ncol = 3, pt.size = 0)
p2 <- FeatureScatter(seuObj, feature1 = "nCount_RNA", feature2 = "pMito_RNA")
p3 <- FeatureScatter(seuObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

ggsave(filename = "SEURAT_YL_QC_1.pdf", plot = p9, width = 8, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_YL_QC_2.pdf", plot = p2, width = 5, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_YL_QC_3.pdf", plot = p3, width = 5, height = 4, units = "in", dpi = 150)

##------------------------------------
## Data Normalization
seuObj <- NormalizeData(seuObj)

##------------------------------------
## Identify top 2000 highly variable genes
seuObj <- FindVariableFeatures(seuObj, selection.method = "vst", nfeatures = 2000)

## Identify the 10 most highly variable genes
top90 <- head(VariableFeatures(seuObj), 10)

## plot variable features with and without labels
p4 <- VariableFeaturePlot(seuObj)
p5 <- LabelPoints(plot = p4, points = top90, repel = TRUE)

ggsave(filename = "SEURAT_YL_VARGENES_1.pdf", plot = p4, width = 8, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_YL_VARGENES_2.pdf", plot = p5, width = 8, height = 4, units = "in", dpi = 150)

##------------------------------------
## Data Scaling
# all.genes <- rownames(seuObj)
seuObj <- ScaleData(seuObj, vars.to.regress = c("nCount_RNA"))

##------------------------------------
## PCA & Jackstraw
## Compute PCA
seuObj <- RunPCA(seuObj, features = VariableFeatures(object = seuObj), npcs = 100)

## Compute and Score Jackstraw
seuObj <- JackStraw(seuObj, num.replicate = 100, dims = 100)
seuObj <- ScoreJackStraw(seuObj, dims = 1:100)

## Plot PCA loadings and Jackstraw scores
p6 <- ElbowPlot(seuObj, ndims = 100)
p7 <- JackStrawPlot(seuObj, dims = 1:100)

ggsave(filename = "SEURAT_YL_PCA_1.pdf", plot = p6, width = 6, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_YL_PCA_2.pdf", plot = p7, width = 12, height = 6, units = "in", dpi = 150)

##------------------------------------
## Data Clustering
## IDENTIFY PCS
pcScores <- seuObj@reductions$pca@jackstraw$overall.p.values
pcScoresNoSign <- pcScores[pcScores[,2] > 0.05,]
selpcs <- min(pcScoresNoSign[,1]) - 1
print(selpcs)
# 59

seuObj <- FindNeighbors(seuObj, dims = 1:selpcs)
seuObj <- FindClusters(seuObj, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 0.8, 1.4, 1.6, 1.8, 2.0))
seuObj <- RunUMAP(seuObj, dims = 1:selpcs)

##------------------------------------
## Save RData
save(seuObj, file = "YL_P07_WT_CB_CLUST.RData")


##------------------------------------
## DOUBLET FINDER
##------------------------------------
## Libraries
rm(list = ls())
library(dplyr)
library(patchwork)
library(ggplot2)
library(reshape2)
library(clustree)
library(Seurat)
library(DoubletFinder)
set.seed(10)

load("YL_P07_WT_CB_CLUST.RData")
# seuObj

pcScores <- seuObj@reductions$pca@jackstraw$overall.p.values
pcScoresNoSign <- pcScores[pcScores[,2] > 0.05,]
selpcs <- min(pcScoresNoSign[,1]) - 1
print(selpcs)
# 59

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_str <- paramSweep_v3(seuObj, PCs = 1:selpcs, sct = FALSE)
sweep.stats_str <- summarizeSweep(sweep.res.list_str, GT = FALSE)
bcmvn_str <- find.pK(sweep.stats_str)
selpk <- 0.005

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(seuObj@meta.data$RNA_snn_res.0.8) ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.1*nrow(seuObj@meta.data))  ## Assuming 10% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seuObjDF <- doubletFinder_v3(seuObj, PCs = 1:selpcs, pN = 0.25, pK = selpk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

seuObjDF <- doubletFinder_v3(seuObjDF, PCs = 1:selpcs, pN = 0.25, pK = selpk, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_2804", sct = FALSE)

table(seuObjDF$GenoAgeSample, seuObjDF$DF.classifications_0.25_0.005_2804)        
#             Doublet Singlet
#   KO_P07_01    1390    8122
#   KO_P07_09     646    9953
#   KO_P07_18     768    7166

table(seuObjDF$GenoAgeSample, seuObjDF$DF.classifications_0.25_0.005_2634)             
#             Doublet Singlet
#   KO_P07_01    1308    8204
#   KO_P07_09     611    9988
#   KO_P07_18     715    7219

seuMeta <- as.data.frame(seuObjDF@meta.data)
seuMetaFiltTemp <- seuMeta[seuMeta$DF.classifications_0.25_0.005_2804 == "Singlet",]
seuMetaFilt <- seuMetaFiltTemp[seuMetaFiltTemp$DF.classifications_0.25_0.005_2634 == "Singlet",]

## Before Doublets Removal
table(seuMeta$GenoAgeSample)
# KO_P07_01 KO_P07_09 KO_P07_18 
#      9512     10599      7934

## After Doublets Removal
table(seuMetaFilt$GenoAgeSample)
# KO_P07_01 KO_P07_09 KO_P07_18 
#      8122      9953      7166

seuMetaFilt$DoubletFinder <- rep("Retained", nrow(seuMetaFilt))

df.meta <- seuMetaFilt$DoubletFinder
names(df.meta) <- row.names(seuMetaFilt)

seuObjDF$DoubletFinder <- df.meta
seuObjDF$DoubletFinder[is.na(seuObjDF$DoubletFinder)] <- "Discarded"

##------------------------------------
## Save RData
save(seuObjDF, file = "YL_P07_KO_CB_DF.RData")
save(seuMeta, seuMetaFilt, file = "YL_P07_KO_CB_DF_FILT_META.RData")
```
<br></br>


## P07 HU
- python/3.7.x-anaconda
- gcc/8.3.0
- hdf5_18/1.8.17
- R/4.0.2-gccmkl
```{R}
##------------------------------------
## SEURAT CLUSTERING
##------------------------------------
## Libraries
rm(list = ls())
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(reshape2)
library(clustree)
set.seed(10)

##------------------------------------
## reference gene symbol list
ref <- read.table("gencode.vM17.protein_coding_gene_id_symbol.txt.gz", header = TRUE, sep = "\t")

##------------------------------------
## Read CellBender H5
p07hu03 <- Read10X_h5("CELLBENDER_OUTPUT/P07_HU/YL_P07_HU_03_CB/CellBender_Out_filtered.h5") # 21967  10538
colnames(p07hu03) <- paste("YL_P07_HU_03", colnames(p07hu03), sep = "_")

p07hu05 <- Read10X_h5("CELLBENDER_OUTPUT/P07_HU/YL_P07_HU_05_CB/CellBender_Out_filtered.h5") # 21967  8668
colnames(p07hu05) <- paste("YL_P07_HU_05", colnames(p07hu05), sep = "_")

p07hu12 <- Read10X_h5("CELLBENDER_OUTPUT/P07_HU/YL_P07_HU_12_CB/CellBender_Out_filtered.h5") # 21967  11589
colnames(p07hu12) <- paste("YL_P07_HU_12", colnames(p07hu12), sep = "_")

p07hu17 <- Read10X_h5("CELLBENDER_OUTPUT/P07_HU/YL_P07_HU_17_CB/CellBender_Out_filtered.h5") # 21967  13
colnames(p07hu17) <- paste("YL_P07_HU_17", colnames(p07hu17), sep = "_")

##------------------------------------
## fetch genes or rows corresponding to gencode protein coding gene symbols
p07hu03.temp <- p07hu03[row.names(p07hu03) %in% ref$GeneSymbol,]
p07hu05.temp <- p07hu05[row.names(p07hu05) %in% ref$GeneSymbol,]
p07hu12.temp <- p07hu12[row.names(p07hu12) %in% ref$GeneSymbol,]
p07hu17.temp <- p07hu17[row.names(p07hu17) %in% ref$GeneSymbol,]

##------------------------------------
## make a data from from matrix
YL_P07_HU_03 <- as.data.frame(as.matrix(p07hu03.temp))
YL_P07_HU_05 <- as.data.frame(as.matrix(p07hu05.temp))
YL_P07_HU_12 <- as.data.frame(as.matrix(p07hu12.temp))
YL_P07_HU_17 <- as.data.frame(as.matrix(p07hu17.temp))

##------------------------------------
## add gene symbol as a column
YL_P07_HU_03$Genes <- row.names(YL_P07_HU_03)
YL_P07_HU_05$Genes <- row.names(YL_P07_HU_05)
YL_P07_HU_12$Genes <- row.names(YL_P07_HU_12)
YL_P07_HU_17$Genes <- row.names(YL_P07_HU_17)

##------------------------------------
## combine individual tables into a giant data frame
## Sample "YL_P07_HU_17" is dropped being an outlier 
dataCombined <- list("YL_P07_HU_03" = YL_P07_HU_03, "YL_P07_HU_05" = YL_P07_HU_05, "YL_P07_HU_12" = YL_P07_HU_12)

combinedData <- Reduce(function(x, y) { merge(x, y, all = TRUE, by = "Genes") } , dataCombined)
row.names(combinedData) <- combinedData$Genes
combinedData$Genes <- NULL

yl.data <- combinedData[row.names(combinedData) %in% ref$GeneSymbol,]

yl.data[is.na(yl.data)] <- 0
print(dim(yl.data)) ## 21967 30808

##------------------------------------
## Save Count Table
p07.hu.data <- yl.data
save(p07.hu.data, file = "YL_P07_HU_CB_COUNTS.RData")

##------------------------------------
## Initialize the Seurat object with the raw (non-normalized data).
seuObj <- CreateSeuratObject(counts = p07.hu.data, project = "YL_P07_HU_CB")
print(dim(seuObj)) ## 21967 30808

##------------------------------------
## Add and update meta data
metaTemp <- as.data.frame(seuObj@meta.data)
metaTemp2 <- as.data.frame(matrix(unlist(strsplit(row.names(metaTemp), "_")), ncol = 5, byrow = TRUE))
row.names(metaTemp2) <- row.names(metaTemp)
metaData <- merge(metaTemp, metaTemp2, by = "row.names")
row.names(metaData) <- metaData$Row.names
metaData$Row.names <- NULL
colnames(metaData) <- c("Ident", "nUMI", "nGenes", "Project", "Age", "Genotype", "Sample", "CellBarcode")

metaSample <- metaData$Sample
names(metaSample) <- row.names(metaData)

metaAge <- metaData$Age
names(metaAge) <- row.names(metaData)

metaGeno <- metaData$Genotype
names(metaGeno) <- row.names(metaData)

seuObj$Genotype <- metaGeno
seuObj$Sample <- metaSample
seuObj$Age <- metaAge
seuObj$GenoAgeSample <- paste(seuObj$Genotype, seuObj$Age, seuObj$Sample, sep = "_")

##------------------------------------
## Calculate Percent Mito
seuObj[["pMito_RNA"]] <- PercentageFeatureSet(seuObj, pattern = "^mt-")

## Set default identities to GenoAgeSample
Idents(seuObj) <- "GenoAgeSample"

table(seuObj@active.ident)
# HU_P07_03 HU_P07_05 HU_P07_12
#     10538      8668     11589

##------------------------------------
## Visualize Data QC
p9 <- VlnPlot(seuObj, features = c("nFeature_RNA", "nCount_RNA", "pMito_RNA"), ncol = 3, pt.size = 0)
p2 <- FeatureScatter(seuObj, feature1 = "nCount_RNA", feature2 = "pMito_RNA")
p3 <- FeatureScatter(seuObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

ggsave(filename = "SEURAT_YL_QC_1.pdf", plot = p9, width = 8, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_YL_QC_2.pdf", plot = p2, width = 5, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_YL_QC_3.pdf", plot = p3, width = 5, height = 4, units = "in", dpi = 150)

##------------------------------------
## Data Normalization
seuObj <- NormalizeData(seuObj)

##------------------------------------
## Identify top 2000 highly variable genes
seuObj <- FindVariableFeatures(seuObj, selection.method = "vst", nfeatures = 2000)

## Identify the 10 most highly variable genes
top90 <- head(VariableFeatures(seuObj), 10)

## plot variable features with and without labels
p4 <- VariableFeaturePlot(seuObj)
p5 <- LabelPoints(plot = p4, points = top90, repel = TRUE)

ggsave(filename = "SEURAT_YL_VARGENES_1.pdf", plot = p4, width = 8, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_YL_VARGENES_2.pdf", plot = p5, width = 8, height = 4, units = "in", dpi = 150)

##------------------------------------
## Data Scaling
# all.genes <- rownames(seuObj)
seuObj <- ScaleData(seuObj, vars.to.regress = c("nCount_RNA"))

##------------------------------------
## PCA & Jackstraw
## Compute PCA
seuObj <- RunPCA(seuObj, features = VariableFeatures(object = seuObj), npcs = 100)

## Compute and Score Jackstraw
seuObj <- JackStraw(seuObj, num.replicate = 100, dims = 100)
seuObj <- ScoreJackStraw(seuObj, dims = 1:100)

## Plot PCA loadings and Jackstraw scores
p6 <- ElbowPlot(seuObj, ndims = 100)
p7 <- JackStrawPlot(seuObj, dims = 1:100)

ggsave(filename = "SEURAT_YL_PCA_1.pdf", plot = p6, width = 6, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_YL_PCA_2.pdf", plot = p7, width = 12, height = 6, units = "in", dpi = 150)

##------------------------------------
## Data Clustering
## IDENTIFY PCS
pcScores <- seuObj@reductions$pca@jackstraw$overall.p.values
pcScoresNoSign <- pcScores[pcScores[,2] > 0.05,]
selpcs <- min(pcScoresNoSign[,1]) - 1
print(selpcs)
# 56

seuObj <- FindNeighbors(seuObj, dims = 1:selpcs)
seuObj <- FindClusters(seuObj, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 0.8, 1.4, 1.6, 1.8, 2.0))
seuObj <- RunUMAP(seuObj, dims = 1:selpcs)

##------------------------------------
## Save RData
save(seuObj, file = "YL_P07_HU_CB_CLUST.RData")



##------------------------------------
## DOUBLET FINDER
##------------------------------------
## Libraries
rm(list = ls())
library(dplyr)
library(patchwork)
library(ggplot2)
library(reshape2)
library(clustree)
library(Seurat)
library(DoubletFinder)
set.seed(10)

load("YL_P07_HU_CB_CLUST.RData")
# seuObj

pcScores <- seuObj@reductions$pca@jackstraw$overall.p.values
pcScoresNoSign <- pcScores[pcScores[,2] > 0.05,]
selpcs <- min(pcScoresNoSign[,1]) - 1
print(selpcs)
# 56

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_str <- paramSweep_v3(seuObj, PCs = 1:selpcs, sct = FALSE)
sweep.stats_str <- summarizeSweep(sweep.res.list_str, GT = FALSE)
bcmvn_str <- find.pK(sweep.stats_str)
selpk <- 0.005

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(seuObj@meta.data$RNA_snn_res.0.8) ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.1*nrow(seuObj@meta.data))  ## Assuming 10% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seuObjDF <- doubletFinder_v3(seuObj, PCs = 1:selpcs, pN = 0.25, pK = selpk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
seuObjDF <- doubletFinder_v3(seuObjDF, PCs = 1:selpcs, pN = 0.25, pK = selpk, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_3080", sct = FALSE)

table(seuObjDF$GenoAgeSample, seuObjDF$DF.classifications_0.25_0.005_3080)        
#             Doublet Singlet
#   HU_P07_03     725    9813
#   HU_P07_05    1118    7550
#   HU_P07_12    1237   10352

table(seuObjDF$GenoAgeSample, seuObjDF$DF.classifications_0.25_0.005_2886)             
#             Doublet Singlet
#   HU_P07_03     671    9867
#   HU_P07_05    1075    7593
#   HU_P07_12    1140   10449

seuMeta <- as.data.frame(seuObjDF@meta.data)
seuMetaFiltTemp <- seuMeta[seuMeta$DF.classifications_0.25_0.005_3080 == "Singlet",]
seuMetaFilt <- seuMetaFiltTemp[seuMetaFiltTemp$DF.classifications_0.25_0.005_2886 == "Singlet",]

## Before Doublets Removal
table(seuMeta$GenoAgeSample)
# HU_P07_03 HU_P07_05 HU_P07_12 
#     10538      8668     11589

## After Doublets Removal
table(seuMetaFilt$GenoAgeSample)
# HU_P07_03 HU_P07_05 HU_P07_12 
#      9813      7550     10352

seuMetaFilt$DoubletFinder <- rep("Retained", nrow(seuMetaFilt))

df.meta <- seuMetaFilt$DoubletFinder
names(df.meta) <- row.names(seuMetaFilt)

seuObjDF$DoubletFinder <- df.meta
seuObjDF$DoubletFinder[is.na(seuObjDF$DoubletFinder)] <- "Discarded"

##------------------------------------
## Save RData
save(seuObjDF, file = "YL_P07_HU_CB_DF.RData")
save(seuMeta, seuMetaFilt, file = "YL_P07_HU_CB_DF_FILT_META.RData")
```
<br></br>



## P56 WT
- python/3.7.x-anaconda
- gcc/8.3.0
- hdf5_18/1.8.17
- R/4.0.2-gccmkl
```{R}
##------------------------------------
## SEURAT CLUSTERING
##------------------------------------
## Libraries
rm(list = ls())
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(reshape2)
library(clustree)
set.seed(10)

##------------------------------------
## reference gene symbol list
ref <- read.table("gencode.vM17.protein_coding_gene_id_symbol.txt.gz", header = TRUE, sep = "\t")

##------------------------------------
## Read CellBender H5
p56wt01 <- Read10X_h5("/work/Neuroinformatics_Core/akulk1/MANUAL/YL_DATA/CELLBENDER_MANUAL/P56_WT/YL_P56_WT_01_CB/CellBender_Out_filtered.h5") # 21967  10029
colnames(p56wt01) <- paste("YL_P56_WT_01", colnames(p56wt01), sep = "_")

p56wt08 <- Read10X_h5("/work/Neuroinformatics_Core/akulk1/MANUAL/YL_DATA/CELLBENDER_MANUAL/P56_WT/YL_P56_WT_08_CB/CellBender_Out_filtered.h5") # 21967  12924
colnames(p56wt08) <- paste("YL_P56_WT_08", colnames(p56wt08), sep = "_")

p56wt10 <- Read10X_h5("/work/Neuroinformatics_Core/akulk1/MANUAL/YL_DATA/CELLBENDER_MANUAL/P56_WT/YL_P56_WT_10_CB/CellBender_Out_filtered.h5") # 21967  12592
colnames(p56wt10) <- paste("YL_P56_WT_10", colnames(p56wt10), sep = "_")

p56wt12 <- Read10X_h5("/work/Neuroinformatics_Core/akulk1/MANUAL/YL_DATA/CELLBENDER_MANUAL/P56_WT/YL_P56_WT_12_CB/CellBender_Out_filtered.h5") #  21967  8882
colnames(p56wt12) <- paste("YL_P56_WT_12", colnames(p56wt12), sep = "_")

##------------------------------------
## fetch genes or rows corresponding to gencode protein coding gene symbols
p56wt01.temp <- p56wt01[row.names(p56wt01) %in% ref$GeneSymbol,]
p56wt08.temp <- p56wt08[row.names(p56wt08) %in% ref$GeneSymbol,]
p56wt10.temp <- p56wt10[row.names(p56wt10) %in% ref$GeneSymbol,]
p56wt12.temp <- p56wt12[row.names(p56wt12) %in% ref$GeneSymbol,]

##------------------------------------
## make a data from from matrix
YL_P56_WT_01 <- as.data.frame(as.matrix(p56wt01.temp))
YL_P56_WT_08 <- as.data.frame(as.matrix(p56wt08.temp))
YL_P56_WT_10 <- as.data.frame(as.matrix(p56wt10.temp))
YL_P56_WT_12 <- as.data.frame(as.matrix(p56wt12.temp))

##------------------------------------
## add gene symbol as a column
YL_P56_WT_01$Genes <- row.names(YL_P56_WT_01)
YL_P56_WT_08$Genes <- row.names(YL_P56_WT_08)
YL_P56_WT_10$Genes <- row.names(YL_P56_WT_10)
YL_P56_WT_12$Genes <- row.names(YL_P56_WT_12)

##------------------------------------
## combine individual tables into a giant data frame
dataCombined <- list("YL_P56_WT_01" = YL_P56_WT_01, "YL_P56_WT_08" = YL_P56_WT_08, "YL_P56_WT_10" = YL_P56_WT_10, "YL_P56_WT_12" = YL_P56_WT_12)

combinedData <- Reduce(function(x, y) { merge(x, y, all = TRUE, by = "Genes") } , dataCombined)
row.names(combinedData) <- combinedData$Genes
combinedData$Genes <- NULL

yl.data <- combinedData[row.names(combinedData) %in% ref$GeneSymbol,]

yl.data[is.na(yl.data)] <- 0
print(dim(yl.data)) ## 21967 44427

##------------------------------------
## Save Count Table
p56.wt.data <- yl.data
save(p56.wt.data, file = "YL_P56_WT_CB_COUNTS.RData")

##------------------------------------
## Initialize the Seurat object with the raw (non-normalized data).
seuObj <- CreateSeuratObject(counts = p56.wt.data, project = "YL_P56_WT_CB")
print(dim(seuObj)) ## 21967 44427

##------------------------------------
## Add and update meta data
metaTemp <- as.data.frame(seuObj@meta.data)
metaTemp2 <- as.data.frame(matrix(unlist(strsplit(row.names(metaTemp), "_")), ncol = 5, byrow = TRUE))
row.names(metaTemp2) <- row.names(metaTemp)
metaData <- merge(metaTemp, metaTemp2, by = "row.names")
row.names(metaData) <- metaData$Row.names
metaData$Row.names <- NULL
colnames(metaData) <- c("Ident", "nUMI", "nGenes", "Project", "Age", "Genotype", "Sample", "CellBarcode")

metaSample <- metaData$Sample
names(metaSample) <- row.names(metaData)

metaAge <- metaData$Age
names(metaAge) <- row.names(metaData)

metaGeno <- metaData$Genotype
names(metaGeno) <- row.names(metaData)

seuObj$Genotype <- metaGeno
seuObj$Sample <- metaSample
seuObj$Age <- metaAge
seuObj$GenoAgeSample <- paste(seuObj$Genotype, seuObj$Age, seuObj$Sample, sep = "_")

##------------------------------------
## Calculate Percent Mito
seuObj[["pMito_RNA"]] <- PercentageFeatureSet(seuObj, pattern = "^mt-")

## Set default identities to GenoAgeSample
Idents(seuObj) <- "GenoAgeSample"

table(seuObj@active.ident)
# WT_P56_01 WT_P56_08 WT_P56_10 WT_P56_12 
#     10029     12924     12592      8882

##------------------------------------
## Visualize Data QC
p9 <- VlnPlot(seuObj, features = c("nFeature_RNA", "nCount_RNA", "pMito_RNA"), ncol = 3, pt.size = 0)
p2 <- FeatureScatter(seuObj, feature1 = "nCount_RNA", feature2 = "pMito_RNA")
p3 <- FeatureScatter(seuObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

ggsave(filename = "SEURAT_YL_QC_1.pdf", plot = p9, width = 8, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_YL_QC_2.pdf", plot = p2, width = 5, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_YL_QC_3.pdf", plot = p3, width = 5, height = 4, units = "in", dpi = 150)

##------------------------------------
## Data Normalization
seuObj <- NormalizeData(seuObj)

##------------------------------------
## Identify top 2000 highly variable genes
seuObj <- FindVariableFeatures(seuObj, selection.method = "vst", nfeatures = 2000)

## Identify the 10 most highly variable genes
top90 <- head(VariableFeatures(seuObj), 10)

## plot variable features with and without labels
p4 <- VariableFeaturePlot(seuObj)
p5 <- LabelPoints(plot = p4, points = top90, repel = TRUE)

ggsave(filename = "SEURAT_YL_VARGENES_1.pdf", plot = p4, width = 8, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_YL_VARGENES_2.pdf", plot = p5, width = 8, height = 4, units = "in", dpi = 150)

##------------------------------------
## Data Scaling
# all.genes <- rownames(seuObj)
seuObj <- ScaleData(seuObj, vars.to.regress = c("nCount_RNA"))

##------------------------------------
## PCA & Jackstraw
## Compute PCA
seuObj <- RunPCA(seuObj, features = VariableFeatures(object = seuObj), npcs = 100)

## Compute and Score Jackstraw
seuObj <- JackStraw(seuObj, num.replicate = 100, dims = 100)
seuObj <- ScoreJackStraw(seuObj, dims = 1:100)

## Plot PCA loadings and Jackstraw scores
p6 <- ElbowPlot(seuObj, ndims = 100)
p7 <- JackStrawPlot(seuObj, dims = 1:100)

ggsave(filename = "SEURAT_YL_PCA_1.pdf", plot = p6, width = 6, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_YL_PCA_2.pdf", plot = p7, width = 12, height = 6, units = "in", dpi = 150)

##------------------------------------
## Data Clustering
## IDENTIFY PCS
pcScores <- seuObj@reductions$pca@jackstraw$overall.p.values
pcScoresNoSign <- pcScores[pcScores[,2] > 0.05,]
selpcs <- min(pcScoresNoSign[,1]) - 1
print(selpcs)
# 67

seuObj <- FindNeighbors(seuObj, dims = 1:selpcs)
seuObj <- FindClusters(seuObj, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 0.8, 1.4, 1.6, 1.8, 2.0))
seuObj <- RunUMAP(seuObj, dims = 1:selpcs)

##------------------------------------
## Save RData
save(seuObj, file = "YL_P56_WT_CB_CLUST.RData")



##------------------------------------
## DOUBLET FINDER
##------------------------------------
## Libraries
rm(list = ls())
library(dplyr)
library(patchwork)
library(ggplot2)
library(reshape2)
library(clustree)
library(Seurat)
library(DoubletFinder)
set.seed(10)

load("YL_P56_WT_CB_CLUST.RData")
# seuObj

pcScores <- seuObj@reductions$pca@jackstraw$overall.p.values
pcScoresNoSign <- pcScores[pcScores[,2] > 0.05,]
selpcs <- min(pcScoresNoSign[,1]) - 1
print(selpcs)
# 67

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_str <- paramSweep_v3(seuObj, PCs = 1:selpcs, sct = FALSE)
sweep.stats_str <- summarizeSweep(sweep.res.list_str, GT = FALSE)
bcmvn_str <- find.pK(sweep.stats_str)
selpk <- 0.005

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(seuObj@meta.data$RNA_snn_res.0.8) ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.1*nrow(seuObj@meta.data))  ## Assuming 10% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seuObjDF <- doubletFinder_v3(seuObj, PCs = 1:selpcs, pN = 0.25, pK = selpk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

seuObjDF <- doubletFinder_v3(seuObjDF, PCs = 1:selpcs, pN = 0.25, pK = selpk, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_4443", sct = FALSE)

table(seuObjDF$GenoAgeSample, seuObjDF$DF.classifications_0.25_0.005_4443)        
#             Doublet Singlet
#   WT_P56_01    1110    8919
#   WT_P56_08    1167   11757
#   WT_P56_10    1153   11439
#   WT_P56_12    1013    7869

table(seuObjDF$GenoAgeSample, seuObjDF$DF.classifications_0.25_0.005_4142)             
#             Doublet Singlet
#   WT_P56_01    1054    8975
#   WT_P56_08    1100   11824
#   WT_P56_10    1062   11530
#   WT_P56_12     926    7956

seuMeta <- as.data.frame(seuObjDF@meta.data)
seuMetaFiltTemp <- seuMeta[seuMeta$DF.classifications_0.25_0.005_4443 == "Singlet",]
seuMetaFilt <- seuMetaFiltTemp[seuMetaFiltTemp$DF.classifications_0.25_0.005_4142 == "Singlet",]

## Before Doublets Removal
table(seuMeta$GenoAgeSample)
# WT_P56_01 WT_P56_08 WT_P56_10 WT_P56_12 
#     10029     12924     12592      8882

## After Doublets Removal
table(seuMetaFilt$GenoAgeSample)
# WT_P56_01 WT_P56_08 WT_P56_10 WT_P56_12 
#      8919     11757     11439      7869

seuMetaFilt$DoubletFinder <- rep("Retained", nrow(seuMetaFilt))

df.meta <- seuMetaFilt$DoubletFinder
names(df.meta) <- row.names(seuMetaFilt)

seuObjDF$DoubletFinder <- df.meta
seuObjDF$DoubletFinder[is.na(seuObjDF$DoubletFinder)] <- "Discarded"

##------------------------------------
## Save RData
save(seuObjDF, file = "YL_P56_WT_CB_DF.RData")
save(seuMeta, seuMetaFilt, file = "YL_P56_WT_CB_DF_FILT_META.RData")
```
<br></br>


## P56 KO
- python/3.7.x-anaconda
- gcc/8.3.0
- hdf5_18/1.8.17
- R/4.0.2-gccmkl
```{R}
##------------------------------------
## SEURAT CLUSTERING
##------------------------------------
## Libraries
rm(list = ls())
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(reshape2)
library(clustree)
set.seed(10)

##------------------------------------
## reference gene symbol list
ref <- read.table("gencode.vM17.protein_coding_gene_id_symbol.txt.gz", header = TRUE, sep = "\t")

##------------------------------------
## Read CellBender H5

#### P56 KO
p56ko03 <- Read10X_h5("CELLBENDER_OUTPUT/P56_KO/YL_P56_KO_03_CB/CellBender_Out_filtered.h5") #  21967  9304
colnames(p56ko03) <- paste("YL_P56_KO_03", colnames(p56ko03), sep = "_")

p56ko04 <- Read10X_h5("CELLBENDER_OUTPUT/P56_KO/YL_P56_KO_04_CB/CellBender_Out_filtered.h5") # 21967  10521
colnames(p56ko04) <- paste("YL_P56_KO_04", colnames(p56ko04), sep = "_")

p56ko07 <- Read10X_h5("CELLBENDER_OUTPUT/P56_KO/YL_P56_KO_07_CB/CellBender_Out_filtered.h5") # 21967  11679
colnames(p56ko07) <- paste("YL_P56_KO_07", colnames(p56ko07), sep = "_")

p56ko13 <- Read10X_h5("CELLBENDER_OUTPUT/P56_KO/YL_P56_KO_13_CB/CellBender_Out_filtered.h5") #  21967  8596
colnames(p56ko13) <- paste("YL_P56_KO_13", colnames(p56ko13), sep = "_")

##------------------------------------
## fetch genes or rows corresponding to gencode protein coding gene symbols
p56ko03.temp <- p56ko03[row.names(p56ko03) %in% ref$GeneSymbol,]
p56ko04.temp <- p56ko04[row.names(p56ko04) %in% ref$GeneSymbol,]
p56ko07.temp <- p56ko07[row.names(p56ko07) %in% ref$GeneSymbol,]
p56ko13.temp <- p56ko13[row.names(p56ko13) %in% ref$GeneSymbol,]

##------------------------------------
## make a data from from matrix
YL_P56_KO_03 <- as.data.frame(as.matrix(p56ko03.temp))
YL_P56_KO_04 <- as.data.frame(as.matrix(p56ko04.temp))
YL_P56_KO_07 <- as.data.frame(as.matrix(p56ko07.temp))
YL_P56_KO_13 <- as.data.frame(as.matrix(p56ko13.temp))

##------------------------------------
## add gene symbol as a column
YL_P56_KO_03$Genes <- row.names(YL_P56_KO_03)
YL_P56_KO_04$Genes <- row.names(YL_P56_KO_04)
YL_P56_KO_07$Genes <- row.names(YL_P56_KO_07)
YL_P56_KO_13$Genes <- row.names(YL_P56_KO_13)

##------------------------------------
## combine individual tables into a giant data frame
dataCombined <- list("YL_P56_KO_03" = YL_P56_KO_03, "YL_P56_KO_04" = YL_P56_KO_04, "YL_P56_KO_07" = YL_P56_KO_07, "YL_P56_KO_13" = YL_P56_KO_13)

combinedData <- Reduce(function(x, y) { merge(x, y, all = TRUE, by = "Genes") } , dataCombined)
row.names(combinedData) <- combinedData$Genes
combinedData$Genes <- NULL

yl.data <- combinedData[row.names(combinedData) %in% ref$GeneSymbol,]

yl.data[is.na(yl.data)] <- 0
print(dim(yl.data)) ## 21967 40100

##------------------------------------
## Save Count Table
p56.ko.data <- yl.data
save(p56.ko.data, file = "YL_P56_KO_CB_COUNTS.RData")

##------------------------------------
## Initialize the Seurat object with the raw (non-normalized data).
seuObj <- CreateSeuratObject(counts = p56.ko.data, project = "YL_P56_KO_CB")
print(dim(seuObj)) ## 21967 40100

##------------------------------------
## Add and update meta data
metaTemp <- as.data.frame(seuObj@meta.data)
metaTemp2 <- as.data.frame(matrix(unlist(strsplit(row.names(metaTemp), "_")), ncol = 5, byrow = TRUE))
row.names(metaTemp2) <- row.names(metaTemp)
metaData <- merge(metaTemp, metaTemp2, by = "row.names")
row.names(metaData) <- metaData$Row.names
metaData$Row.names <- NULL
colnames(metaData) <- c("Ident", "nUMI", "nGenes", "Project", "Age", "Genotype", "Sample", "CellBarcode")

metaSample <- metaData$Sample
names(metaSample) <- row.names(metaData)

metaAge <- metaData$Age
names(metaAge) <- row.names(metaData)

metaGeno <- metaData$Genotype
names(metaGeno) <- row.names(metaData)

seuObj$Genotype <- metaGeno
seuObj$Sample <- metaSample
seuObj$Age <- metaAge
seuObj$GenoAgeSample <- paste(seuObj$Genotype, seuObj$Age, seuObj$Sample, sep = "_")

##------------------------------------
## Calculate Percent Mito
seuObj[["pMito_RNA"]] <- PercentageFeatureSet(seuObj, pattern = "^mt-")

## Set default identities to GenoAgeSample
Idents(seuObj) <- "GenoAgeSample"

table(seuObj@active.ident)
# KO_P56_03 KO_P56_04 KO_P56_07 KO_P56_13 
#      9304     10521     11679      8596

## Visualize Data QC
p9 <- VlnPlot(seuObj, features = c("nFeature_RNA", "nCount_RNA", "pMito_RNA"), ncol = 3, pt.size = 0)
p2 <- FeatureScatter(seuObj, feature1 = "nCount_RNA", feature2 = "pMito_RNA")
p3 <- FeatureScatter(seuObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

ggsave(filename = "SEURAT_YL_QC_1.pdf", plot = p9, width = 8, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_YL_QC_2.pdf", plot = p2, width = 5, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_YL_QC_3.pdf", plot = p3, width = 5, height = 4, units = "in", dpi = 150)

##------------------------------------
## Data Normalization
seuObj <- NormalizeData(seuObj)

##------------------------------------
## Identify top 2000 highly variable genes
seuObj <- FindVariableFeatures(seuObj, selection.method = "vst", nfeatures = 2000)

## Identify the 10 most highly variable genes
top90 <- head(VariableFeatures(seuObj), 10)

## plot variable features with and without labels
p4 <- VariableFeaturePlot(seuObj)
p5 <- LabelPoints(plot = p4, points = top90, repel = TRUE)

ggsave(filename = "SEURAT_YL_VARGENES_1.pdf", plot = p4, width = 8, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_YL_VARGENES_2.pdf", plot = p5, width = 8, height = 4, units = "in", dpi = 150)

##------------------------------------
## Data Scaling
# all.genes <- rownames(seuObj)
seuObj <- ScaleData(seuObj, vars.to.regress = c("nCount_RNA"))

##------------------------------------
## PCA & Jackstraw
## Compute PCA
seuObj <- RunPCA(seuObj, features = VariableFeatures(object = seuObj), npcs = 100)

## Compute and Score Jackstraw
seuObj <- JackStraw(seuObj, num.replicate = 100, dims = 100)
seuObj <- ScoreJackStraw(seuObj, dims = 1:100)

## Plot PCA loadings and Jackstraw scores
p6 <- ElbowPlot(seuObj, ndims = 100)
p7 <- JackStrawPlot(seuObj, dims = 1:100)

ggsave(filename = "SEURAT_YL_PCA_1.pdf", plot = p6, width = 6, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_YL_PCA_2.pdf", plot = p7, width = 12, height = 6, units = "in", dpi = 150)

##------------------------------------
## Data Clustering
## IDENTIFY PCS
pcScores <- seuObj@reductions$pca@jackstraw$overall.p.values
pcScoresNoSign <- pcScores[pcScores[,2] > 0.05,]
selpcs <- min(pcScoresNoSign[,1]) - 1
print(selpcs)
# 64

seuObj <- FindNeighbors(seuObj, dims = 1:selpcs)
seuObj <- FindClusters(seuObj, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 0.8, 1.4, 1.6, 1.8, 2.0))
seuObj <- RunUMAP(seuObj, dims = 1:selpcs)

##------------------------------------
## Save RData
save(seuObj, file = "YL_P56_KO_CB_CLUST.RData")



##------------------------------------
## DOUBLET FINDER
##------------------------------------
## Libraries
rm(list = ls())
library(dplyr)
library(patchwork)
library(ggplot2)
library(reshape2)
library(clustree)
library(Seurat)
library(DoubletFinder)
set.seed(10)

load("YL_P56_KO_CB_CLUST.RData")
# seuObj

pcScores <- seuObj@reductions$pca@jackstraw$overall.p.values
pcScoresNoSign <- pcScores[pcScores[,2] > 0.05,]
selpcs <- min(pcScoresNoSign[,1]) - 1
print(selpcs)
# 58

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_str <- paramSweep_v3(seuObj, PCs = 1:selpcs, sct = FALSE)
sweep.stats_str <- summarizeSweep(sweep.res.list_str, GT = FALSE)
bcmvn_str <- find.pK(sweep.stats_str)
selpk <- 0.005

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(seuObj@meta.data$RNA_snn_res.0.8) ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.1*nrow(seuObj@meta.data))  ## Assuming 10% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seuObjDF <- doubletFinder_v3(seuObj, PCs = 1:selpcs, pN = 0.25, pK = selpk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

seuObjDF <- doubletFinder_v3(seuObjDF, PCs = 1:selpcs, pN = 0.25, pK = selpk, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_4010", sct = FALSE)

table(seuObjDF$GenoAgeSample, seuObjDF$DF.classifications_0.25_0.005_4010)        
#             Doublet Singlet
#   KO_P56_03     789    8515
#   KO_P56_04    1400    9121
#   KO_P56_07     722   10957
#   KO_P56_13    1099    7497

table(seuObjDF$GenoAgeSample, seuObjDF$DF.classifications_0.25_0.005_3743)             
#             Doublet Singlet
#   KO_P56_03     740    8564
#   KO_P56_04    1309    9212
#   KO_P56_07     665   11014
#   KO_P56_13    1029    7567

seuMeta <- as.data.frame(seuObjDF@meta.data)
seuMetaFiltTemp <- seuMeta[seuMeta$DF.classifications_0.25_0.005_4010 == "Singlet",]
seuMetaFilt <- seuMetaFiltTemp[seuMetaFiltTemp$DF.classifications_0.25_0.005_3743 == "Singlet",]

## Before Doublets Removal
table(seuMeta$GenoAgeSample)
# KO_P56_03 KO_P56_04 KO_P56_07 KO_P56_13 
#      9304     10521     11679      8596

## After Doublets Removal
table(seuMetaFilt$GenoAgeSample)
# KO_P56_03 KO_P56_04 KO_P56_07 KO_P56_13 
#      8515      9121     10957      7497

seuMetaFilt$DoubletFinder <- rep("Retained", nrow(seuMetaFilt))

df.meta <- seuMetaFilt$DoubletFinder
names(df.meta) <- row.names(seuMetaFilt)

seuObjDF$DoubletFinder <- df.meta
seuObjDF$DoubletFinder[is.na(seuObjDF$DoubletFinder)] <- "Discarded"

##------------------------------------
## Save RData
save(seuObjDF, file = "YL_P56_KO_CB_DF.RData")
save(seuMeta, seuMetaFilt, file = "YL_P56_KO_CB_DF_FILT_META.RData")
```
<br></br>


## P56 HU
- python/3.7.x-anaconda
- gcc/8.3.0
- hdf5_18/1.8.17
- R/4.0.2-gccmkl
```{R}
##------------------------------------
## SEURAT CLUSTERING
##------------------------------------
## Libraries
rm(list = ls())
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(reshape2)
library(clustree)
set.seed(10)

##------------------------------------
## reference gene symbol list
ref <- read.table("gencode.vM17.protein_coding_gene_id_symbol.txt.gz", header = TRUE, sep = "\t")

##------------------------------------
## Read CellBender H5
p56hu02 <- Read10X_h5("CELLBENDER_OUTPUT/P56_HU/YL_P56_HU_02_CB/CellBender_Out_filtered.h5") #  21967  9904
colnames(p56hu02) <- paste("YL_P56_HU_02", colnames(p56hu02), sep = "_")

p56hu05 <- Read10X_h5("CELLBENDER_OUTPUT/P56_HU/YL_P56_HU_05_CB/CellBender_Out_filtered.h5") #  21967  7608
colnames(p56hu05) <- paste("YL_P56_HU_05", colnames(p56hu05), sep = "_")

p56hu09 <- Read10X_h5("CELLBENDER_OUTPUT/P56_HU/YL_P56_HU_09_CB/CellBender_Out_filtered.h5") # 21967  14343
colnames(p56hu09) <- paste("YL_P56_HU_09", colnames(p56hu09), sep = "_")

p56hu11 <- Read10X_h5("CELLBENDER_OUTPUT/P56_HU/YL_P56_HU_11_CB/CellBender_Out_filtered.h5") #  21967  8256
colnames(p56hu11) <- paste("YL_P56_HU_11", colnames(p56hu11), sep = "_")

##------------------------------------
## fetch genes or rows corresponding to gencode protein coding gene symbols
p56hu02.temp <- p56hu02[row.names(p56hu02) %in% ref$GeneSymbol,]
p56hu05.temp <- p56hu05[row.names(p56hu05) %in% ref$GeneSymbol,]
p56hu09.temp <- p56hu09[row.names(p56hu09) %in% ref$GeneSymbol,]
p56hu11.temp <- p56hu11[row.names(p56hu11) %in% ref$GeneSymbol,]

##------------------------------------
## make a data from from matrix
YL_P56_HU_02 <- as.data.frame(as.matrix(p56hu02.temp))
YL_P56_HU_05 <- as.data.frame(as.matrix(p56hu05.temp))
YL_P56_HU_09 <- as.data.frame(as.matrix(p56hu09.temp))
YL_P56_HU_11 <- as.data.frame(as.matrix(p56hu11.temp))

##------------------------------------
## add gene symbol as a column
YL_P56_HU_02$Genes <- row.names(YL_P56_HU_02)
YL_P56_HU_05$Genes <- row.names(YL_P56_HU_05)
YL_P56_HU_09$Genes <- row.names(YL_P56_HU_09)
YL_P56_HU_11$Genes <- row.names(YL_P56_HU_11)

##------------------------------------
## combine individual tables into a giant data frame
dataCombined <- list("YL_P56_HU_02" = YL_P56_HU_02, "YL_P56_HU_05" = YL_P56_HU_05, "YL_P56_HU_09" = YL_P56_HU_09, "YL_P56_HU_11" = YL_P56_HU_11)

combinedData <- Reduce(function(x, y) { merge(x, y, all = TRUE, by = "Genes") } , dataCombined)
row.names(combinedData) <- combinedData$Genes
combinedData$Genes <- NULL

yl.data <- combinedData[row.names(combinedData) %in% ref$GeneSymbol,]

yl.data[is.na(yl.data)] <- 0
print(dim(yl.data)) ## 21967 40111

##------------------------------------
## Save Count Table
p56.hu.data <- yl.data
save(p56.hu.data, file = "YL_P56_HU_CB_COUNTS.RData")

##------------------------------------
## Initialize the Seurat object with the raw (non-normalized data).
seuObj <- CreateSeuratObject(counts = p56.hu.data, project = "YL_P56_HU_CB")
print(dim(seuObj)) ## 21967 40111

##------------------------------------
## Add and update meta data
metaTemp <- as.data.frame(seuObj@meta.data)
metaTemp2 <- as.data.frame(matrix(unlist(strsplit(row.names(metaTemp), "_")), ncol = 5, byrow = TRUE))
row.names(metaTemp2) <- row.names(metaTemp)
metaData <- merge(metaTemp, metaTemp2, by = "row.names")
row.names(metaData) <- metaData$Row.names
metaData$Row.names <- NULL
colnames(metaData) <- c("Ident", "nUMI", "nGenes", "Project", "Age", "Genotype", "Sample", "CellBarcode")

metaSample <- metaData$Sample
names(metaSample) <- row.names(metaData)

metaAge <- metaData$Age
names(metaAge) <- row.names(metaData)

metaGeno <- metaData$Genotype
names(metaGeno) <- row.names(metaData)

seuObj$Genotype <- metaGeno
seuObj$Sample <- metaSample
seuObj$Age <- metaAge
seuObj$GenoAgeSample <- paste(seuObj$Genotype, seuObj$Age, seuObj$Sample, sep = "_")

##------------------------------------
## Calculate Percent Mito
seuObj[["pMito_RNA"]] <- PercentageFeatureSet(seuObj, pattern = "^mt-")

## Set default identities to GenoAgeSample
Idents(seuObj) <- "GenoAgeSample"

table(seuObj@active.ident)
# HU_P56_02 HU_P56_05 HU_P56_09 HU_P56_11 
#      9904      7608     14343      8256

##------------------------------------
## Visualize Data QC
p9 <- VlnPlot(seuObj, features = c("nFeature_RNA", "nCount_RNA", "pMito_RNA"), ncol = 3, pt.size = 0)
p2 <- FeatureScatter(seuObj, feature1 = "nCount_RNA", feature2 = "pMito_RNA")
p3 <- FeatureScatter(seuObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

ggsave(filename = "SEURAT_YL_QC_1.pdf", plot = p9, width = 8, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_YL_QC_2.pdf", plot = p2, width = 5, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_YL_QC_3.pdf", plot = p3, width = 5, height = 4, units = "in", dpi = 150)

##------------------------------------
## Data Normalization
seuObj <- NormalizeData(seuObj)

##------------------------------------
## Identify top 2000 highly variable genes
seuObj <- FindVariableFeatures(seuObj, selection.method = "vst", nfeatures = 2000)

## Identify the 10 most highly variable genes
top90 <- head(VariableFeatures(seuObj), 10)

## plot variable features with and without labels
p4 <- VariableFeaturePlot(seuObj)
p5 <- LabelPoints(plot = p4, points = top90, repel = TRUE)

ggsave(filename = "SEURAT_YL_VARGENES_1.pdf", plot = p4, width = 8, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_YL_VARGENES_2.pdf", plot = p5, width = 8, height = 4, units = "in", dpi = 150)

##------------------------------------
## Data Scaling
# all.genes <- rownames(seuObj)
seuObj <- ScaleData(seuObj, vars.to.regress = c("nCount_RNA"))

##------------------------------------
## PCA & Jackstraw
## Compute PCA
seuObj <- RunPCA(seuObj, features = VariableFeatures(object = seuObj), npcs = 100)

## Compute and Score Jackstraw
seuObj <- JackStraw(seuObj, num.replicate = 100, dims = 100)
seuObj <- ScoreJackStraw(seuObj, dims = 1:100)

## Plot PCA loadings and Jackstraw scores
p6 <- ElbowPlot(seuObj, ndims = 100)
p7 <- JackStrawPlot(seuObj, dims = 1:100)

ggsave(filename = "SEURAT_YL_PCA_1.pdf", plot = p6, width = 6, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_YL_PCA_2.pdf", plot = p7, width = 12, height = 6, units = "in", dpi = 150)

##------------------------------------
## Data Clustering
## IDENTIFY PCS
pcScores <- seuObj@reductions$pca@jackstraw$overall.p.values
pcScoresNoSign <- pcScores[pcScores[,2] > 0.05,]
selpcs <- min(pcScoresNoSign[,1]) - 1
print(selpcs)
# 64

seuObj <- FindNeighbors(seuObj, dims = 1:selpcs)
seuObj <- FindClusters(seuObj, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 0.8, 1.4, 1.6, 1.8, 2.0))
seuObj <- RunUMAP(seuObj, dims = 1:selpcs)

##------------------------------------
## Save RData
save(seuObj, file = "YL_P56_HU_CB_CLUST.RData")



##------------------------------------
## DOUBLET FINDER
##------------------------------------
## Libraries
rm(list = ls())
library(dplyr)
library(patchwork)
library(ggplot2)
library(reshape2)
library(clustree)
library(Seurat)
library(DoubletFinder)
set.seed(10)

load("YL_P56_HU_CB_CLUST.RData")
# seuObj

pcScores <- seuObj@reductions$pca@jackstraw$overall.p.values
pcScoresNoSign <- pcScores[pcScores[,2] > 0.05,]
selpcs <- min(pcScoresNoSign[,1]) - 1
print(selpcs)
# 64

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_str <- paramSweep_v3(seuObj, PCs = 1:selpcs, sct = FALSE)
sweep.stats_str <- summarizeSweep(sweep.res.list_str, GT = FALSE)
bcmvn_str <- find.pK(sweep.stats_str)
selpk <- 0.005

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(seuObj@meta.data$RNA_snn_res.0.8) ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.1*nrow(seuObj@meta.data))  ## Assuming 10% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seuObjDF <- doubletFinder_v3(seuObj, PCs = 1:selpcs, pN = 0.25, pK = selpk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

seuObjDF <- doubletFinder_v3(seuObjDF, PCs = 1:selpcs, pN = 0.25, pK = selpk, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_4011", sct = FALSE)

table(seuObjDF$GenoAgeSample, seuObjDF$DF.classifications_0.25_0.005_4011)        
#             Doublet Singlet
#   HU_P56_02     754    9150
#   HU_P56_05    1419    6189
#   HU_P56_09     874   13469
#   HU_P56_11     964    7292

table(seuObjDF$GenoAgeSample, seuObjDF$DF.classifications_0.25_0.005_3738)             
#             Doublet Singlet
#   HU_P56_02     684    9220
#   HU_P56_05    1346    6262
#   HU_P56_09     807   13536
#   HU_P56_11     901    7355

seuMeta <- as.data.frame(seuObjDF@meta.data)
seuMetaFiltTemp <- seuMeta[seuMeta$DF.classifications_0.25_0.005_4011 == "Singlet",]
seuMetaFilt <- seuMetaFiltTemp[seuMetaFiltTemp$DF.classifications_0.25_0.005_3738 == "Singlet",]

## Before Doublets Removal
table(seuMeta$GenoAgeSample)
# HU_P56_02 HU_P56_05 HU_P56_09 HU_P56_11 
#      9904      7608     14343      8256

## After Doublets Removal
table(seuMetaFilt$GenoAgeSample)
# HU_P56_02 HU_P56_05 HU_P56_09 HU_P56_11 
#      9150      6189     13469      7292

seuMetaFilt$DoubletFinder <- rep("Retained", nrow(seuMetaFilt))

df.meta <- seuMetaFilt$DoubletFinder
names(df.meta) <- row.names(seuMetaFilt)

seuObjDF$DoubletFinder <- df.meta
seuObjDF$DoubletFinder[is.na(seuObjDF$DoubletFinder)] <- "Discarded"

##------------------------------------
## Save RData
save(seuObjDF, file = "YL_P56_HU_CB_DF.RData")
save(seuMeta, seuMetaFilt, file = "YL_P56_HU_CB_DF_FILT_META.RData")
```


-----
