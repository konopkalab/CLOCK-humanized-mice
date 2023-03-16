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
NK_P07_WT_02 <- as.data.frame(as.matrix(p07wt02.temp))
NK_P07_WT_04 <- as.data.frame(as.matrix(p07wt04.temp))
NK_P07_WT_11 <- as.data.frame(as.matrix(p07wt11.temp))
NK_P07_WT_16 <- as.data.frame(as.matrix(p07wt16.temp))

##------------------------------------
## add gene symbol as a column
NK_P07_WT_02$Genes <- row.names(NK_P07_WT_02)
NK_P07_WT_04$Genes <- row.names(NK_P07_WT_04)
NK_P07_WT_11$Genes <- row.names(NK_P07_WT_11)
NK_P07_WT_16$Genes <- row.names(NK_P07_WT_16)

##------------------------------------
## combine individual tables into a giant data frame
## Sample "NK_P07_WT_04" is dropped being an outlier
dataCombined <- list("NK_P07_WT_02" = NK_P07_WT_02, "NK_P07_WT_11" = NK_P07_WT_11, "NK_P07_WT_16" = NK_P07_WT_16)

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

ggsave(filename = "SEURAT_NK_QC_1.pdf", plot = p9, width = 8, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_NK_QC_2.pdf", plot = p2, width = 5, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_NK_QC_3.pdf", plot = p3, width = 5, height = 4, units = "in", dpi = 150)


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

ggsave(filename = "SEURAT_NK_VARGENES_1.pdf", plot = p4, width = 8, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_NK_VARGENES_2.pdf", plot = p5, width = 8, height = 4, units = "in", dpi = 150)


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

ggsave(filename = "SEURAT_NK_PCA_1.pdf", plot = p6, width = 6, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_NK_PCA_2.pdf", plot = p7, width = 12, height = 6, units = "in", dpi = 150)


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

