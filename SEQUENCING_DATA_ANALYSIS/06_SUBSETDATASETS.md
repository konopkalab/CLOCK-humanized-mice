# SUBSET DATASETS


## P07 WT
- python/3.7.x-anaconda
- gcc/8.3.0
- hdf5_18/1.8.17
- R/4.0.2-gccmkl
```{R}
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
## CellBender and DoubletFinder Cleaned Seurat Object
load("DOUBLETFINDER/B_USING_CELLBENDER_P07_WT/YL_P07_WT_CB_DF.RData")
# seuObjDF

table(seuObjDF@meta.data$DoubletFinder)
# Discarded  Retained 
#      2480     22321

cells2keep <- row.names(seuObjDF@meta.data[seuObjDF@meta.data$DoubletFinder == "Retained",])

##------------------------------------
## reference gene symbol list
ref <- read.table("gencode.vM17.protein_coding_gene_id_symbol.txt.gz", header = TRUE, sep = "\t")

##------------------------------------
## CELLBENDER MANUAL
#### P07 WT
p07wt02 <- Read10X_h5("/work/Neuroinformatics_Core/akulk1/MANUAL/YL_DATA/CELLBENDER_MANUAL/P07_WT/YL_P07_WT_02_CB/CellBender_Out_filtered.h5") # 21967  6990
colnames(p07wt02) <- paste("YL_P07_WT_02", colnames(p07wt02), sep = "_")

p07wt04 <- Read10X_h5("/work/Neuroinformatics_Core/akulk1/MANUAL/YL_DATA/CELLBENDER_MANUAL/CELLBENDER_AGAIN/YL_P07_WT_04_CB/CellBender_Out_filtered.h5") # 21967     8
colnames(p07wt04) <- paste("YL_P07_WT_04", colnames(p07wt04), sep = "_")

p07wt11 <- Read10X_h5("/work/Neuroinformatics_Core/akulk1/MANUAL/YL_DATA/CELLBENDER_MANUAL/P07_WT/YL_P07_WT_11_CB/CellBender_Out_filtered.h5") # 21967 11406
colnames(p07wt11) <- paste("YL_P07_WT_11", colnames(p07wt11), sep = "_")

p07wt16 <- Read10X_h5("/work/Neuroinformatics_Core/akulk1/MANUAL/YL_DATA/CELLBENDER_MANUAL/P07_WT/YL_P07_WT_16_CB/CellBender_Out_filtered.h5") # 21967  6405
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
dataCombined <- list("YL_P07_WT_02" = YL_P07_WT_02, "YL_P07_WT_11" = YL_P07_WT_11, "YL_P07_WT_16" = YL_P07_WT_16)

combinedData <- Reduce(function(x, y) { merge(x, y, all = TRUE, by = "Genes") } , dataCombined)
row.names(combinedData) <- combinedData$Genes
combinedData$Genes <- NULL

yl.data <- combinedData[row.names(combinedData) %in% ref$GeneSymbol,]

yl.data[is.na(yl.data)] <- 0
print(dim(yl.data)) ## 21967 24801

p07.wt.data <- yl.data[,colnames(yl.data) %in% cells2keep]

##------------------------------------
## Save count table
save(p07.wt.data, file = "YL_P07_WT_CB_DF_COUNTS.RData")

##------------------------------------
## Initialize the Seurat object with the raw (non-normalized data).
seuObj <- CreateSeuratObject(counts = p07.wt.data, project = "YL_P07_WT_CB_DF")
print(dim(seuObj)) ## 21967 22321

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

print(table(seuObj@active.ident))

##------------------------------------
## Visualize Data QC
p9 <- VlnPlot(seuObj, features = c("nFeature_RNA", "nCount_RNA", "pMito_RNA"), ncol = 3, pt.size = 0)
p2 <- FeatureScatter(seuObj, feature1 = "nCount_RNA", feature2 = "pMito_RNA")
p3 <- FeatureScatter(seuObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

ggsave(filename = "SEURAT_YL_QC_1.pdf", plot = p9, width = 8, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_YL_QC_2.pdf", plot = p2, width = 5, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_YL_QC_3.pdf", plot = p3, width = 5, height = 4, units = "in", dpi = 150)

##------------------------------------
## Data Filtering
seuObj <- subset(seuObj, subset = nCount_RNA < 10000 & pMito_RNA < 0.5)
print(dim(seuObj)) ## 21967 40128

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
all.genes <- rownames(seuObj)
seuObj <- ScaleData(seuObj, features = all.genes, vars.to.regress = c("nCount_RNA", "pMito_RNA"))

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
# 46

seuObj <- FindNeighbors(seuObj, dims = 1:selpcs)
seuObj <- FindClusters(seuObj, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 0.8, 1.4, 1.6, 1.8, 2.0))
seuObj <- RunUMAP(seuObj, dims = 1:selpcs)

##------------------------------------
## Save RData
save(seuObj, file = "YL_P07_WT_CB_DF_CLUST.RData")

##------------------------------------
## Subset
seuObj_ClockYes_ArntlYes <- subset(x = seuObj, subset = Clock > 0 & Arntl > 0) # 1152
seuObj_ClockNo_ArntlNo <- subset(x = seuObj, subset = Clock <= 0 & Arntl <= 0) # 14517
seuObj_ClockYes_ArntlNo <- subset(x = seuObj, subset = Clock > 0 & Arntl <= 0) # 3798
seuObj_ClockNo_ArntlYes <- subset(x = seuObj, subset = Arntl > 0 & Clock <= 0) # 2675

save(seuObj_ClockYes_ArntlYes, file = "YL_P07_WT_CB_DF_ClockYes_ArntlYes.RData")
save(seuObj_ClockNo_ArntlNo, file = "YL_P07_WT_CB_DF_ClockNo_ArntlNo.RData")
save(seuObj_ClockYes_ArntlNo, file = "YL_P07_WT_CB_DF_ClockYes_ArntlNo.RData")
save(seuObj_ClockNo_ArntlYes, file = "YL_P07_WT_CB_DF_ClockNo_ArntlYes.RData")

```
<br></br>


## P07 KO
- python/3.7.x-anaconda
- gcc/8.3.0
- hdf5_18/1.8.17
- R/4.0.2-gccmkl
```{R}
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
## CellBender and DoubletFinder Cleaned Seurat Object
load("DOUBLETFINDER/B_USING_CELLBENDER_P07_KO/YL_P07_KO_CB_DF.RData")
# seuObjDF

table(seuObjDF@meta.data$DoubletFinder)
# Discarded  Retained 
#      2804     25241

cells2keep <- row.names(seuObjDF@meta.data[seuObjDF@meta.data$DoubletFinder == "Retained",])

##------------------------------------
## reference gene symbol list
ref <- read.table("gencode.vM17.protein_coding_gene_id_symbol.txt.gz", header = TRUE, sep = "\t")

##------------------------------------
## CELLBENDER MANUAL
#### P07 KO
p07ko01 <- Read10X_h5("/work/Neuroinformatics_Core/akulk1/MANUAL/YL_DATA/CELLBENDER_MANUAL/P07_KO/YL_P07_KO_01_CB/CellBender_Out_filtered.h5") # 21967  9512
colnames(p07ko01) <- paste("YL_P07_KO_01", colnames(p07ko01), sep = "_")

p07ko06 <- Read10X_h5("/work/Neuroinformatics_Core/akulk1/MANUAL/YL_DATA/CELLBENDER_MANUAL/P07_KO/YL_P07_KO_06_CB/CellBender_Out_filtered.h5") # 21967  16
colnames(p07ko06) <- paste("YL_P07_KO_06", colnames(p07ko06), sep = "_")

p07ko09 <- Read10X_h5("/work/Neuroinformatics_Core/akulk1/MANUAL/YL_DATA/CELLBENDER_MANUAL/P07_KO/YL_P07_KO_09_CB/CellBender_Out_filtered.h5") # 21967  10599
colnames(p07ko09) <- paste("YL_P07_KO_09", colnames(p07ko09), sep = "_")

p07ko18 <- Read10X_h5("/work/Neuroinformatics_Core/akulk1/MANUAL/YL_DATA/CELLBENDER_MANUAL/P07_KO/YL_P07_KO_18_CB/CellBender_Out_filtered.h5") # 21967  7934
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
dataCombined <- list("YL_P07_KO_01" = YL_P07_KO_01, "YL_P07_KO_09" = YL_P07_KO_09, "YL_P07_KO_18" = YL_P07_KO_18)

combinedData <- Reduce(function(x, y) { merge(x, y, all = TRUE, by = "Genes") } , dataCombined)
row.names(combinedData) <- combinedData$Genes
combinedData$Genes <- NULL

yl.data <- combinedData[row.names(combinedData) %in% ref$GeneSymbol,]

yl.data[is.na(yl.data)] <- 0
print(dim(yl.data)) ## 21967 28045

p07.ko.data <- yl.data[,colnames(yl.data) %in% cells2keep]
print(dim(p07.ko.data))
# 21967 25241

##------------------------------------
## Save count table
save(p07.ko.data, file = "YL_P07_KO_CB_DF_COUNTS.RData")

##------------------------------------
## Initialize the Seurat object with the raw (non-normalized data).
seuObj <- CreateSeuratObject(counts = p07.ko.data, project = "YL_P07_KO_CB_DF")
print(dim(seuObj)) 

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

print(table(seuObj@active.ident))

## Visualize Data QC
p9 <- VlnPlot(seuObj, features = c("nFeature_RNA", "nCount_RNA", "pMito_RNA"), ncol = 3, pt.size = 0)
p2 <- FeatureScatter(seuObj, feature1 = "nCount_RNA", feature2 = "pMito_RNA")
p3 <- FeatureScatter(seuObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

ggsave(filename = "SEURAT_YL_QC_1.pdf", plot = p9, width = 8, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_YL_QC_2.pdf", plot = p2, width = 5, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_YL_QC_3.pdf", plot = p3, width = 5, height = 4, units = "in", dpi = 150)

##------------------------------------
## Data Filtering
seuObj <- subset(seuObj, subset = nCount_RNA < 10000 & pMito_RNA < 0.5)
print(dim(seuObj)) ## 21967 40128

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
all.genes <- rownames(seuObj)
seuObj <- ScaleData(seuObj, features = all.genes, vars.to.regress = c("nCount_RNA", "pMito_RNA"))

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

seuObj <- FindNeighbors(seuObj, dims = 1:selpcs)
seuObj <- FindClusters(seuObj, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 0.8, 1.4, 1.6, 1.8, 2.0))
seuObj <- RunUMAP(seuObj, dims = 1:selpcs)

##------------------------------------
## Save RData
save(seuObj, file = "YL_P07_KO_CB_DF_CLUST.RData")

##------------------------------------
## Subset
seuObj_ClockYes_ArntlYes <- subset(x = seuObj, subset = Clock > 0 & Arntl > 0) # 1951
seuObj_ClockNo_ArntlNo <- subset(x = seuObj, subset = Clock <= 0 & Arntl <= 0) # 14750
seuObj_ClockYes_ArntlNo <- subset(x = seuObj, subset = Clock > 0 & Arntl <= 0) # 4837
seuObj_ClockNo_ArntlYes <- subset(x = seuObj, subset = Arntl > 0 & Clock <= 0) # 3564

save(seuObj_ClockYes_ArntlYes, file = "YL_P07_KO_CB_DF_ClockYes_ArntlYes.RData")
save(seuObj_ClockNo_ArntlNo, file = "YL_P07_KO_CB_DF_ClockNo_ArntlNo.RData")
save(seuObj_ClockYes_ArntlNo, file = "YL_P07_KO_CB_DF_ClockYes_ArntlNo.RData")
save(seuObj_ClockNo_ArntlYes, file = "YL_P07_KO_CB_DF_ClockNo_ArntlYes.RData")

```
<br></br>



## P07 HU
- python/3.7.x-anaconda
- gcc/8.3.0
- hdf5_18/1.8.17
- R/4.0.2-gccmkl
```{R}
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
## COUNT TABLE AFTER CB and DF
p07hu03_m <- read.table("COUNTS_M/YL_P07_HU_03_Counts.tsv.gz", sep = "\t", header = TRUE, row.names = 1)
colnames(p07hu03_m) <- paste("YL_P07_HU_03", colnames(p07hu03_m), sep = "_")
p07hu03_h <- read.table("SEURAT_AFTER_CB_DF/P07_HU/COUNTS_H/YL_P07_HU_03_Counts.tsv.gz", sep = "\t", header = TRUE, row.names = 1)
colnames(p07hu03_h) <- paste("YL_P07_HU_03", colnames(p07hu03_h), sep = "_")

p07hu05_m <- read.table("SEURAT_AFTER_CB_DF/P07_HU/COUNTS_M/YL_P07_HU_05_Counts.tsv.gz", sep = "\t", header = TRUE, row.names = 1)
colnames(p07hu05_m) <- paste("YL_P07_HU_05", colnames(p07hu05_m), sep = "_")
p07hu05_h <- read.table("SEURAT_AFTER_CB_DF/P07_HU/COUNTS_H/YL_P07_HU_05_Counts.tsv.gz", sep = "\t", header = TRUE, row.names = 1)
colnames(p07hu05_h) <- paste("YL_P07_HU_05", colnames(p07hu05_h), sep = "_")

p07hu12_m <- read.table("SEURAT_AFTER_CB_DF/P07_HU/COUNTS_M/YL_P07_HU_12_Counts.tsv.gz", sep = "\t", header = TRUE, row.names = 1)
colnames(p07hu12_m) <- paste("YL_P07_HU_12", colnames(p07hu12_m), sep = "_")
p07hu12_h <- read.table("SEURAT_AFTER_CB_DF/P07_HU/COUNTS_H/YL_P07_HU_12_Counts.tsv.gz", sep = "\t", header = TRUE, row.names = 1)
colnames(p07hu12_h) <- paste("YL_P07_HU_12", colnames(p07hu12_h), sep = "_")

p07hu17_m <- read.table("SEURAT_AFTER_CB_DF/P07_HU/COUNTS_M/YL_P07_HU_17_Counts.tsv.gz", sep = "\t", header = TRUE, row.names = 1)
colnames(p07hu17_m) <- paste("YL_P07_HU_17", colnames(p07hu17_m), sep = "_")
p07hu17_h <- read.table("SEURAT_AFTER_CB_DF/P07_HU/COUNTS_H/YL_P07_HU_17_Counts.tsv.gz", sep = "\t", header = TRUE, row.names = 1)
colnames(p07hu17_h) <- paste("YL_P07_HU_17", colnames(p07hu17_h), sep = "_")

##------------------------------------
## CellBender and DoubletFinder Cleaned Seurat Object
load("DOUBLETFINDER/B_USING_CELLBENDER_P07_HU/YL_P07_HU_CB_DF.RData")
# seuObjDF

table(seuObjDF@meta.data$DoubletFinder)
# Discarded  Retained 
#      3080     27715

cells2keep <- row.names(seuObjDF@meta.data[seuObjDF@meta.data$DoubletFinder == "Retained",])

##------------------------------------
## CELLBENDER MANUAL
#### P07 HU
p07hu03 <- Read10X_h5("CELLBENDER_MANUAL/P07_HU/YL_P07_HU_03_CB/CellBender_Out_filtered.h5") # 21967  10538
colnames(p07hu03) <- paste("YL_P07_HU_03", colnames(p07hu03), sep = "_")

p07hu05 <- Read10X_h5("CELLBENDER_MANUAL/P07_HU/YL_P07_HU_05_CB/CellBender_Out_filtered.h5") # 21967  8668
colnames(p07hu05) <- paste("YL_P07_HU_05", colnames(p07hu05), sep = "_")

p07hu12 <- Read10X_h5("CELLBENDER_MANUAL/P07_HU/YL_P07_HU_12_CB/CellBender_Out_filtered.h5") # 21967  11589
colnames(p07hu12) <- paste("YL_P07_HU_12", colnames(p07hu12), sep = "_")

p07hu17 <- Read10X_h5("CELLBENDER_MANUAL/P07_HU/YL_P07_HU_17_CB/CellBender_Out_filtered.h5") # 21967  13
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
## Sample "YL_P07_HU_17" dropped as an outlier
dataCombined <- list("YL_P07_HU_03" = YL_P07_HU_03, "YL_P07_HU_05" = YL_P07_HU_05, "YL_P07_HU_12" = YL_P07_HU_12)

combinedData <- Reduce(function(x, y) { merge(x, y, all = TRUE, by = "Genes") } , dataCombined)
row.names(combinedData) <- combinedData$Genes
combinedData$Genes <- NULL

yl.data <- combinedData[row.names(combinedData) %in% ref$GeneSymbol,]

yl.data[is.na(yl.data)] <- 0
print(dim(yl.data)) ## 21967 30795

p07.hu.data <- yl.data[,colnames(yl.data) %in% cells2keep]
print(dim(p07.hu.data))
# 21967 27715

##------------------------------------
## Save count table
save(p07.hu.data, file = "YL_P07_HU_CB_DF_COUNTS.RData")

##------------------------------------
## organize for human and mouse clock
p7_hu_03 <- p07.hu.data[,grepl("^YL_P07_HU_03", colnames(p07.hu.data))]
p7_hu_05 <- p07.hu.data[,grepl("^YL_P07_HU_05", colnames(p07.hu.data))]
p7_hu_12 <- p07.hu.data[,grepl("^YL_P07_HU_12", colnames(p07.hu.data))]

p7_hu_03_tempm <- as.data.frame(t(p07hu03_m[row.names(p07hu03_m) %in% c("Clock"), colnames(p07hu03_m) %in% colnames(p7_hu_03)]))
p7_hu_03_temph <- as.data.frame(t(p07hu03_h[row.names(p07hu03_h) %in% c("CLOCK"), colnames(p07hu03_h) %in% colnames(p7_hu_03)]))
p7_hu_03_w <- merge(p7_hu_03_tempm, p7_hu_03_temph, by = "row.names", all.x = TRUE, all.y = TRUE)
row.names(p7_hu_03_w) <- p7_hu_03_w$Row.names
p7_hu_03_w$Row.names <- NULL
colnames(p7_hu_03_w) <- c("Clock", "Clock_human")

p7_hu_05_tempm <- as.data.frame(t(p07hu05_m[row.names(p07hu05_m) %in% c("Clock"), colnames(p07hu05_m) %in% colnames(p7_hu_05)]))
p7_hu_05_temph <- as.data.frame(t(p07hu05_h[row.names(p07hu05_h) %in% c("CLOCK"), colnames(p07hu05_h) %in% colnames(p7_hu_05)]))
p7_hu_05_w <- merge(p7_hu_05_tempm, p7_hu_05_temph, by = "row.names", all.x = TRUE, all.y = TRUE)
row.names(p7_hu_05_w) <- p7_hu_05_w$Row.names
p7_hu_05_w$Row.names <- NULL
colnames(p7_hu_05_w) <- c("Clock", "Clock_human")

p7_hu_12_tempm <- as.data.frame(t(p07hu12_m[row.names(p07hu12_m) %in% c("Clock"), colnames(p07hu12_m) %in% colnames(p7_hu_12)]))
p7_hu_12_temph <- as.data.frame(t(p07hu12_h[row.names(p07hu12_h) %in% c("CLOCK"), colnames(p07hu12_h) %in% colnames(p7_hu_12)]))
p7_hu_12_w <- merge(p7_hu_12_tempm, p7_hu_12_temph, by = "row.names", all.x = TRUE, all.y = TRUE)
row.names(p7_hu_12_w) <- p7_hu_12_w$Row.names
p7_hu_12_w$Row.names <- NULL
colnames(p7_hu_12_w) <- c("Clock", "Clock_human")

p7_hu_03_wo <- as.data.frame(t(p7_hu_03[!row.names(p7_hu_03) %in% c("Clock"),]))
p7_hu_05_wo <- as.data.frame(t(p7_hu_05[!row.names(p7_hu_05) %in% c("Clock"),]))
p7_hu_12_wo <- as.data.frame(t(p7_hu_12[!row.names(p7_hu_12) %in% c("Clock"),]))

p7_hu_03_updated <- merge(p7_hu_03_wo, p7_hu_03_w, by = "row.names", all.x = TRUE)
row.names(p7_hu_03_updated) <- p7_hu_03_updated$Row.names
p7_hu_03_updated$Row.names <- NULL
p7_hu_03_final <- as.data.frame(t(p7_hu_03_updated))
p7_hu_03_final <- p7_hu_03_final[order(row.names(p7_hu_03_final)),order(colnames(p7_hu_03_final))]
p7_hu_03_final[is.na(p7_hu_03_final)] <- 0

p7_hu_05_updated <- merge(p7_hu_05_wo, p7_hu_05_w, by = "row.names", all.x = TRUE)
row.names(p7_hu_05_updated) <- p7_hu_05_updated$Row.names
p7_hu_05_updated$Row.names <- NULL
p7_hu_05_final <- as.data.frame(t(p7_hu_05_updated))
p7_hu_05_final <- p7_hu_05_final[order(row.names(p7_hu_05_final)),order(colnames(p7_hu_05_final))]
p7_hu_05_final[is.na(p7_hu_05_final)] <- 0

p7_hu_12_updated <- merge(p7_hu_12_wo, p7_hu_12_w, by = "row.names", all.x = TRUE)
row.names(p7_hu_12_updated) <- p7_hu_12_updated$Row.names
p7_hu_12_updated$Row.names <- NULL
p7_hu_12_final <- as.data.frame(t(p7_hu_12_updated))
p7_hu_12_final <- p7_hu_12_final[order(row.names(p7_hu_12_final)),order(colnames(p7_hu_12_final))]
p7_hu_12_final[is.na(p7_hu_12_final)] <- 0

##------------------------------------
## add gene symbol as a column
p7_hu_03_final$Genes <- row.names(p7_hu_03_final)
p7_hu_05_final$Genes <- row.names(p7_hu_05_final)
p7_hu_12_final$Genes <- row.names(p7_hu_12_final)

##------------------------------------
## combine individual tables into a giant data frame
dataCombined <- list("YL_P07_HU_03" = p7_hu_03_final, "YL_P07_HU_05" = p7_hu_05_final, "YL_P07_HU_12" = p7_hu_12_final)

combinedData <- Reduce(function(x, y) { merge(x, y, all = TRUE, by = "Genes") } , dataCombined)
row.names(combinedData) <- combinedData$Genes
combinedData$Genes <- NULL

yl.data <- combinedData #[row.names(combinedData) %in% ref$GeneSymbol,]

yl.data[is.na(yl.data)] <- 0
print(dim(yl.data)) ## 21968 27715

p07.hu.data <- yl.data #[,colnames(yl.data) %in% cells2keep]
print(dim(p07.hu.data))
# 21968 27715

##------------------------------------
## Save count table
save(p07.hu.data, file = "YL_P07_HU_CB_DF_COUNTS_UPDATED.RData")

##------------------------------------
## Initialize the Seurat object with the raw (non-normalized data).
seuObj <- CreateSeuratObject(counts = p07.hu.data, project = "YL_P07_HU_CB_DF")
print(dim(seuObj)) ## 21967 22321

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

print(table(seuObj@active.ident))
# HU_P07_03 HU_P07_05 HU_P07_12 
#      9813      7550     10352

##------------------------------------
## Visualize Data QC
p9 <- VlnPlot(seuObj, features = c("nFeature_RNA", "nCount_RNA", "pMito_RNA"), ncol = 3, pt.size = 0)
p2 <- FeatureScatter(seuObj, feature1 = "nCount_RNA", feature2 = "pMito_RNA")
p3 <- FeatureScatter(seuObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

ggsave(filename = "SEURAT_YL_QC_1.pdf", plot = p9, width = 8, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_YL_QC_2.pdf", plot = p2, width = 5, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_YL_QC_3.pdf", plot = p3, width = 5, height = 4, units = "in", dpi = 150)

##------------------------------------
## Data Filtering
seuObj <- subset(seuObj, subset = nCount_RNA < 10000 & pMito_RNA < 0.5)
print(dim(seuObj)) ## 21968 27628

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
all.genes <- rownames(seuObj)
seuObj <- ScaleData(seuObj, features = all.genes, vars.to.regress = c("nCount_RNA", "pMito_RNA"))

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

seuObj <- FindNeighbors(seuObj, dims = 1:selpcs)
seuObj <- FindClusters(seuObj, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 0.8, 1.4, 1.6, 1.8, 2.0))
seuObj <- RunUMAP(seuObj, dims = 1:selpcs)

##------------------------------------
## Save RData
save(seuObj, file = "YL_P07_HU_CB_DF_CLUST.RData")

##------------------------------------
## Subset
seuObj_replace_plus <- subset(x = seuObj, subset = Clock > 0 & `Clock-human` > 0 & Arntl > 0) # 838
seuObj_replace_minus <- subset(x = seuObj, subset = Clock > 0 & `Clock-human` > 0 & Arntl <= 0) # 3039

seuObj_mis_plus <- subset(x = seuObj, subset = Clock <= 0 & `Clock-human` > 0 & Arntl > 0) # 2516
seuObj_mis_minus <- subset(x = seuObj, subset = Clock <= 0 & `Clock-human` > 0 & Arntl <= 0) # 11299

seuObj_ko_plus <- subset(x = seuObj, subset = Clock > 0 & `Clock-human` <= 0 & Arntl > 0) # 297
seuObj_ko_minus <- subset(x = seuObj, subset = Clock > 0 & `Clock-human` <= 0 & Arntl <= 0) # 1481

seuObj_no_plus <- subset(x = seuObj, subset = Clock <= 0 & `Clock-human` <= 0 & Arntl > 0) # 1122
seuObj_no_minus <- subset(x = seuObj, subset = Clock <= 0 & `Clock-human` <= 0 & Arntl <= 0) # 7036


save(seuObj_replace_plus, file = "YL_P07_HU_CB_DF_replace_plus.RData")
save(seuObj_replace_minus, file = "YL_P07_HU_CB_DF_replace_minus.RData")
save(seuObj_mis_plus, file = "YL_P07_HU_CB_DF_mis_plus.RData")
save(seuObj_mis_minus, file = "YL_P07_HU_CB_DF_mis_minus.RData")
save(seuObj_ko_plus, file = "YL_P07_HU_CB_DF_ko_plus.RData")
save(seuObj_ko_minus, file = "YL_P07_HU_CB_DF_ko_minus.RData")
save(seuObj_no_plus, file = "YL_P07_HU_CB_DF_no_plus.RData")
save(seuObj_no_minus, file = "YL_P07_HU_CB_DF_no_minus.RData")

```
<br></br>



## P56 WT
- python/3.7.x-anaconda
- gcc/8.3.0
- hdf5_18/1.8.17
- R/4.0.2-gccmkl
```{R}
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
## CellBender and DoubletFinder Cleaned Seurat Object
load("DOUBLETFINDER/B_USING_CELLBENDER_P56_WT/YL_P56_WT_CB_DF.RData")
# seuObjDF

table(seuObjDF@meta.data$DoubletFinder)
# Discarded  Retained 
#      4443     39984

cells2keep <- row.names(seuObjDF@meta.data[seuObjDF@meta.data$DoubletFinder == "Retained",])

##------------------------------------
## reference gene symbol list
ref <- read.table("gencode.vM17.protein_coding_gene_id_symbol.txt.gz", header = TRUE, sep = "\t")

##------------------------------------
## CELLBENDER MANUAL
#### P56 WT
p56wt01 <- Read10X_h5("CELLBENDER_MANUAL/P56_WT/YL_P56_WT_01_CB/CellBender_Out_filtered.h5") # 21967  10029
colnames(p56wt01) <- paste("YL_P56_WT_01", colnames(p56wt01), sep = "_")

p56wt08 <- Read10X_h5("CELLBENDER_MANUAL/P56_WT/YL_P56_WT_08_CB/CellBender_Out_filtered.h5") # 21967  12924
colnames(p56wt08) <- paste("YL_P56_WT_08", colnames(p56wt08), sep = "_")

p56wt10 <- Read10X_h5("CELLBENDER_MANUAL/P56_WT/YL_P56_WT_10_CB/CellBender_Out_filtered.h5") # 21967  12592
colnames(p56wt10) <- paste("YL_P56_WT_10", colnames(p56wt10), sep = "_")

p56wt12 <- Read10X_h5("CELLBENDER_MANUAL/P56_WT/YL_P56_WT_12_CB/CellBender_Out_filtered.h5") #  21967  8882
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
print(dim(yl.data)) ## 21967 24801

p56.wt.data <- yl.data[,colnames(yl.data) %in% cells2keep]

##------------------------------------
## Save count table
save(p56.wt.data, file = "YL_P56_WT_CB_DF_COUNTS.RData")

##------------------------------------
## Initialize the Seurat object with the raw (non-normalized data).
seuObj <- CreateSeuratObject(counts = p56.wt.data, project = "YL_P56_WT_CB_DF")
print(dim(seuObj)) ## 21967 39984

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

print(table(seuObj@active.ident))

##------------------------------------
## Visualize Data QC
p9 <- VlnPlot(seuObj, features = c("nFeature_RNA", "nCount_RNA", "pMito_RNA"), ncol = 3, pt.size = 0)
p2 <- FeatureScatter(seuObj, feature1 = "nCount_RNA", feature2 = "pMito_RNA")
p3 <- FeatureScatter(seuObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

ggsave(filename = "SEURAT_YL_QC_1.pdf", plot = p9, width = 8, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_YL_QC_2.pdf", plot = p2, width = 5, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_YL_QC_3.pdf", plot = p3, width = 5, height = 4, units = "in", dpi = 150)

##------------------------------------
## Data Filtering
seuObj <- subset(seuObj, subset = nCount_RNA < 10000 & pMito_RNA < 0.5)
print(dim(seuObj)) ## 21967 40128

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
all.genes <- rownames(seuObj)
seuObj <- ScaleData(seuObj, features = all.genes, vars.to.regress = c("nCount_RNA", "pMito_RNA"))

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
# 46

seuObj <- FindNeighbors(seuObj, dims = 1:selpcs)
seuObj <- FindClusters(seuObj, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 0.8, 1.4, 1.6, 1.8, 2.0))
seuObj <- RunUMAP(seuObj, dims = 1:selpcs)

##------------------------------------
## Save RData
save(seuObj, file = "YL_P56_WT_CB_DF_CLUST.RData")

##------------------------------------
## Subset
seuObj_ClockYes_ArntlYes <- subset(x = seuObj, subset = Clock > 0 & Arntl > 0) # 
seuObj_ClockNo_ArntlNo <- subset(x = seuObj, subset = Clock <= 0 & Arntl <= 0) # 
seuObj_ClockYes_ArntlNo <- subset(x = seuObj, subset = Clock > 0 & Arntl <= 0) # 
seuObj_ClockNo_ArntlYes <- subset(x = seuObj, subset = Arntl > 0 & Clock <= 0) # 

save(seuObj_ClockYes_ArntlYes, file = "YL_P56_WT_CB_DF_ClockYes_ArntlYes.RData")
save(seuObj_ClockNo_ArntlNo, file = "YL_P56_WT_CB_DF_ClockNo_ArntlNo.RData")
save(seuObj_ClockYes_ArntlNo, file = "YL_P56_WT_CB_DF_ClockYes_ArntlNo.RData")
save(seuObj_ClockNo_ArntlYes, file = "YL_P56_WT_CB_DF_ClockNo_ArntlYes.RData")

```
<br></br>




## P56 KO
- python/3.7.x-anaconda
- gcc/8.3.0
- hdf5_18/1.8.17
- R/4.0.2-gccmkl
```{R}
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
## CellBender and DoubletFinder Cleaned Seurat Object
load("DOUBLETFINDER/B_USING_CELLBENDER_P56_KO/YL_P56_KO_CB_DF.RData")
# seuObjDF

table(seuObjDF@meta.data$DoubletFinder)
# Discarded  Retained 
#      4010     36090

cells2keep <- row.names(seuObjDF@meta.data[seuObjDF@meta.data$DoubletFinder == "Retained",])

##------------------------------------
## reference gene symbol list
ref <- read.table("gencode.vM17.protein_coding_gene_id_symbol.txt.gz", header = TRUE, sep = "\t")

##------------------------------------
## CELLBENDER MANUAL
#### P56 KO
p56ko03 <- Read10X_h5("CELLBENDER_MANUAL/P56_KO/YL_P56_KO_03_CB/CellBender_Out_filtered.h5") #  21967  9304
colnames(p56ko03) <- paste("YL_P56_KO_03", colnames(p56ko03), sep = "_")

p56ko04 <- Read10X_h5("CELLBENDER_MANUAL/P56_KO/YL_P56_KO_04_CB/CellBender_Out_filtered.h5") # 21967  10521
colnames(p56ko04) <- paste("YL_P56_KO_04", colnames(p56ko04), sep = "_")

p56ko07 <- Read10X_h5("CELLBENDER_MANUAL/P56_KO/YL_P56_KO_07_CB/CellBender_Out_filtered.h5") # 21967  11679
colnames(p56ko07) <- paste("YL_P56_KO_07", colnames(p56ko07), sep = "_")

p56ko13 <- Read10X_h5("CELLBENDER_MANUAL/P56_KO/YL_P56_KO_13_CB/CellBender_Out_filtered.h5") #  21967  8596
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

p56.ko.data <- yl.data[,colnames(yl.data) %in% cells2keep]

##------------------------------------
## Save count table
save(p56.ko.data, file = "YL_P56_KO_CB_DF_COUNTS.RData")

##------------------------------------
## Initialize the Seurat object with the raw (non-normalized data).
seuObj <- CreateSeuratObject(counts = p56.ko.data, project = "YL_P56_KO_CB_DF")
print(dim(seuObj)) ## 21967 36090

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

print(table(seuObj@active.ident))

##------------------------------------
## Visualize Data QC
p9 <- VlnPlot(seuObj, features = c("nFeature_RNA", "nCount_RNA", "pMito_RNA"), ncol = 3, pt.size = 0)
p2 <- FeatureScatter(seuObj, feature1 = "nCount_RNA", feature2 = "pMito_RNA")
p3 <- FeatureScatter(seuObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

ggsave(filename = "SEURAT_YL_QC_1.pdf", plot = p9, width = 8, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_YL_QC_2.pdf", plot = p2, width = 5, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_YL_QC_3.pdf", plot = p3, width = 5, height = 4, units = "in", dpi = 150)

##------------------------------------
## Data Filtering
seuObj <- subset(seuObj, subset = nCount_RNA < 10000 & pMito_RNA < 0.5)
print(dim(seuObj)) ## 21967 40128

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
all.genes <- rownames(seuObj)
seuObj <- ScaleData(seuObj, features = all.genes, vars.to.regress = c("nCount_RNA", "pMito_RNA"))

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
# 46

seuObj <- FindNeighbors(seuObj, dims = 1:selpcs)
seuObj <- FindClusters(seuObj, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 0.8, 1.4, 1.6, 1.8, 2.0))
seuObj <- RunUMAP(seuObj, dims = 1:selpcs)

##------------------------------------
## Save RData
save(seuObj, file = "YL_P56_KO_CB_DF_CLUST.RData")

##------------------------------------
## Subset
seuObj_ClockYes_ArntlYes <- subset(x = seuObj, subset = Clock > 0 & Arntl > 0) # 4426
seuObj_ClockNo_ArntlNo <- subset(x = seuObj, subset = Clock <= 0 & Arntl <= 0) # 17954
seuObj_ClockYes_ArntlNo <- subset(x = seuObj, subset = Clock > 0 & Arntl <= 0) # 6232
seuObj_ClockNo_ArntlYes <- subset(x = seuObj, subset = Arntl > 0 & Clock <= 0) # 5046

save(seuObj_ClockYes_ArntlYes, file = "YL_P56_KO_CB_DF_ClockYes_ArntlYes.RData")
save(seuObj_ClockNo_ArntlNo, file = "YL_P56_KO_CB_DF_ClockNo_ArntlNo.RData")
save(seuObj_ClockYes_ArntlNo, file = "YL_P56_KO_CB_DF_ClockYes_ArntlNo.RData")
save(seuObj_ClockNo_ArntlYes, file = "YL_P56_KO_CB_DF_ClockNo_ArntlYes.RData")

```
<br></br>



## P56 HU
- python/3.7.x-anaconda
- gcc/8.3.0
- hdf5_18/1.8.17
- R/4.0.2-gccmkl
```{R}
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
## CellBender and DoubletFinder Cleaned Seurat Object
load("DOUBLETFINDER/B_USING_CELLBENDER_P56_HU/YL_P56_HU_CB_DF.RData")
# seuObjDF

table(seuObjDF@meta.data$DoubletFinder)
# Discarded  Retained 
#      3080     27715

cells2keep <- row.names(seuObjDF@meta.data[seuObjDF@meta.data$DoubletFinder == "Retained",])

##------------------------------------
## reference gene symbol list
ref <- read.table("gencode.vM17.protein_coding_gene_id_symbol.txt.gz", header = TRUE, sep = "\t")

## CELLBENDER MANUAL
#### P56 HU
p56hu02 <- Read10X_h5("CELLBENDER_MANUAL/P56_HU/YL_P56_HU_02_CB/CellBender_Out_filtered.h5") # 21967  9904
colnames(p56hu02) <- paste("YL_P56_HU_02", colnames(p56hu02), sep = "_")

p56hu05 <- Read10X_h5("CELLBENDER_MANUAL/P56_HU/YL_P56_HU_05_CB/CellBender_Out_filtered.h5") # 21967  7608
colnames(p56hu05) <- paste("YL_P56_HU_05", colnames(p56hu05), sep = "_")

p56hu09 <- Read10X_h5("CELLBENDER_MANUAL/P56_HU/YL_P56_HU_09_CB/CellBender_Out_filtered.h5") # 21967  14343
colnames(p56hu09) <- paste("YL_P56_HU_09", colnames(p56hu09), sep = "_")

p56hu11 <- Read10X_h5("CELLBENDER_MANUAL/P56_HU/YL_P56_HU_11_CB/CellBender_Out_filtered.h5") # 21967  8256
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

p56.hu.data <- yl.data[,colnames(yl.data) %in% cells2keep]
print(dim(p56.hu.data))
# 21967 36100

##------------------------------------
## Raw count tables
p56hu02_m <- read.table("SEURAT_AFTER_CB_DF/P56_HU/COUNTS_M/YL_P56_HU_02_Counts.tsv.gz", sep = "\t", header = TRUE, row.names = 1)
colnames(p56hu02_m) <- paste("YL_P56_HU_02", colnames(p56hu02_m), sep = "_")
p56hu02_h <- read.table("SEURAT_AFTER_CB_DF/P56_HU/COUNTS_H/YL_P56_HU_02_Counts.tsv.gz", sep = "\t", header = TRUE, row.names = 1)
colnames(p56hu02_h) <- paste("YL_P56_HU_02", colnames(p56hu02_h), sep = "_")

p56hu05_m <- read.table("SEURAT_AFTER_CB_DF/P56_HU/COUNTS_M/YL_P56_HU_05_Counts.tsv.gz", sep = "\t", header = TRUE, row.names = 1)
colnames(p56hu05_m) <- paste("YL_P56_HU_05", colnames(p56hu05_m), sep = "_")
p56hu05_h <- read.table("SEURAT_AFTER_CB_DF/P56_HU/COUNTS_H/YL_P56_HU_05_Counts.tsv.gz", sep = "\t", header = TRUE, row.names = 1)
colnames(p56hu05_h) <- paste("YL_P56_HU_05", colnames(p56hu05_h), sep = "_")

p56hu09_m <- read.table("SEURAT_AFTER_CB_DF/P56_HU/COUNTS_M/YL_P56_HU_09_Counts.tsv.gz", sep = "\t", header = TRUE, row.names = 1)
colnames(p56hu09_m) <- paste("YL_P56_HU_09", colnames(p56hu09_m), sep = "_")
p56hu09_h <- read.table("SEURAT_AFTER_CB_DF/P56_HU/COUNTS_H/YL_P56_HU_09_Counts.tsv.gz", sep = "\t", header = TRUE, row.names = 1)
colnames(p56hu09_h) <- paste("YL_P56_HU_09", colnames(p56hu09_h), sep = "_")

p56hu11_m <- read.table("SEURAT_AFTER_CB_DF/P56_HU/COUNTS_M/YL_P56_HU_11_Counts.tsv.gz", sep = "\t", header = TRUE, row.names = 1)
colnames(p56hu11_m) <- paste("YL_P56_HU_11", colnames(p56hu11_m), sep = "_")
p56hu11_h <- read.table("SEURAT_AFTER_CB_DF/P56_HU/COUNTS_H/YL_P56_HU_11_Counts.tsv.gz", sep = "\t", header = TRUE, row.names = 1)
colnames(p56hu11_h) <- paste("YL_P56_HU_11", colnames(p56hu11_h), sep = "_")

##------------------------------------
## organize for human and mouse clock

p56_hu_02 <- p56.hu.data[,grepl("^YL_P56_HU_02", colnames(p56.hu.data))]
p56_hu_05 <- p56.hu.data[,grepl("^YL_P56_HU_05", colnames(p56.hu.data))]
p56_hu_09 <- p56.hu.data[,grepl("^YL_P56_HU_09", colnames(p56.hu.data))]
p56_hu_11 <- p56.hu.data[,grepl("^YL_P56_HU_11", colnames(p56.hu.data))]

p56_hu_02_tempm <- as.data.frame(t(p56hu02_m[row.names(p56hu02_m) %in% c("Clock"), colnames(p56hu02_m) %in% colnames(p56_hu_02)]))
p56_hu_02_temph <- as.data.frame(t(p56hu02_h[row.names(p56hu02_h) %in% c("CLOCK"), colnames(p56hu02_h) %in% colnames(p56_hu_02)]))
p56_hu_02_w <- merge(p56_hu_02_tempm, p56_hu_02_temph, by = "row.names", all.x = TRUE, all.y = TRUE)
row.names(p56_hu_02_w) <- p56_hu_02_w$Row.names
p56_hu_02_w$Row.names <- NULL
colnames(p56_hu_02_w) <- c("Clock", "Clock_human")

p56_hu_05_tempm <- as.data.frame(t(p56hu05_m[row.names(p56hu05_m) %in% c("Clock"), colnames(p56hu05_m) %in% colnames(p56_hu_05)]))
p56_hu_05_temph <- as.data.frame(t(p56hu05_h[row.names(p56hu05_h) %in% c("CLOCK"), colnames(p56hu05_h) %in% colnames(p56_hu_05)]))
p56_hu_05_w <- merge(p56_hu_05_tempm, p56_hu_05_temph, by = "row.names", all.x = TRUE, all.y = TRUE)
row.names(p56_hu_05_w) <- p56_hu_05_w$Row.names
p56_hu_05_w$Row.names <- NULL
colnames(p56_hu_05_w) <- c("Clock", "Clock_human")

p56_hu_09_tempm <- as.data.frame(t(p56hu09_m[row.names(p56hu09_m) %in% c("Clock"), colnames(p56hu09_m) %in% colnames(p56_hu_09)]))
p56_hu_09_temph <- as.data.frame(t(p56hu09_h[row.names(p56hu09_h) %in% c("CLOCK"), colnames(p56hu09_h) %in% colnames(p56_hu_09)]))
p56_hu_09_w <- merge(p56_hu_09_tempm, p56_hu_09_temph, by = "row.names", all.x = TRUE, all.y = TRUE)
row.names(p56_hu_09_w) <- p56_hu_09_w$Row.names
p56_hu_09_w$Row.names <- NULL
colnames(p56_hu_09_w) <- c("Clock", "Clock_human")

p56_hu_11_tempm <- as.data.frame(t(p56hu11_m[row.names(p56hu11_m) %in% c("Clock"), colnames(p56hu11_m) %in% colnames(p56_hu_11)]))
p56_hu_11_temph <- as.data.frame(t(p56hu11_h[row.names(p56hu11_h) %in% c("CLOCK"), colnames(p56hu11_h) %in% colnames(p56_hu_11)]))
p56_hu_11_w <- merge(p56_hu_11_tempm, p56_hu_11_temph, by = "row.names", all.x = TRUE, all.y = TRUE)
row.names(p56_hu_11_w) <- p56_hu_09_w$Row.names
p56_hu_11_w$Row.names <- NULL
colnames(p56_hu_11_w) <- c("Clock", "Clock_human")

p56_hu_02_wo <- as.data.frame(t(p56_hu_02[!row.names(p56_hu_02) %in% c("Clock"),]))
p56_hu_05_wo <- as.data.frame(t(p56_hu_05[!row.names(p56_hu_05) %in% c("Clock"),]))
p56_hu_09_wo <- as.data.frame(t(p56_hu_09[!row.names(p56_hu_09) %in% c("Clock"),]))
p56_hu_11_wo <- as.data.frame(t(p56_hu_11[!row.names(p56_hu_11) %in% c("Clock"),]))

p56_hu_02_updated <- merge(p56_hu_02_wo, p56_hu_02_w, by = "row.names", all.x = TRUE)
row.names(p56_hu_02_updated) <- p56_hu_02_updated$Row.names
p56_hu_02_updated$Row.names <- NULL
p56_hu_02_final <- as.data.frame(t(p56_hu_02_updated))
p56_hu_02_final <- p56_hu_02_final[order(row.names(p56_hu_02_final)),order(colnames(p56_hu_02_final))]
p56_hu_02_final[is.na(p56_hu_02_final)] <- 0

p56_hu_05_updated <- merge(p56_hu_05_wo, p56_hu_05_w, by = "row.names", all.x = TRUE)
row.names(p56_hu_05_updated) <- p56_hu_05_updated$Row.names
p56_hu_05_updated$Row.names <- NULL
p56_hu_05_final <- as.data.frame(t(p56_hu_05_updated))
p56_hu_05_final <- p56_hu_05_final[order(row.names(p56_hu_05_final)),order(colnames(p56_hu_05_final))]
p56_hu_05_final[is.na(p56_hu_05_final)] <- 0

p56_hu_09_updated <- merge(p56_hu_09_wo, p56_hu_09_w, by = "row.names", all.x = TRUE)
row.names(p56_hu_09_updated) <- p56_hu_09_updated$Row.names
p56_hu_09_updated$Row.names <- NULL
p56_hu_09_final <- as.data.frame(t(p56_hu_09_updated))
p56_hu_09_final <- p56_hu_09_final[order(row.names(p56_hu_09_final)),order(colnames(p56_hu_09_final))]
p56_hu_09_final[is.na(p56_hu_09_final)] <- 0

p56_hu_11_updated <- merge(p56_hu_11_wo, p56_hu_11_w, by = "row.names", all.x = TRUE)
row.names(p56_hu_11_updated) <- p56_hu_11_updated$Row.names
p56_hu_11_updated$Row.names <- NULL
p56_hu_11_final <- as.data.frame(t(p56_hu_11_updated))
p56_hu_11_final <- p56_hu_11_final[order(row.names(p56_hu_11_final)),order(colnames(p56_hu_11_final))]
p56_hu_11_final[is.na(p56_hu_11_final)] <- 0

##------------------------------------
## add gene symbol as a column
p56_hu_02_final$Genes <- row.names(p56_hu_02_final)
p56_hu_05_final$Genes <- row.names(p56_hu_05_final)
p56_hu_09_final$Genes <- row.names(p56_hu_09_final)
p56_hu_11_final$Genes <- row.names(p56_hu_11_final)

##------------------------------------
## combine individual tables into a giant data frame
dataCombined <- list("YL_P56_HU_02" = p56_hu_02_final, "YL_P56_HU_05" = p56_hu_05_final, "YL_P56_HU_09" = p56_hu_09_final, "YL_P56_HU_11" = p56_hu_11_final)

combinedData <- Reduce(function(x, y) { merge(x, y, all = TRUE, by = "Genes") } , dataCombined)
row.names(combinedData) <- combinedData$Genes
combinedData$Genes <- NULL

yl.data <- combinedData #[row.names(combinedData) %in% ref$GeneSymbol,]

yl.data[is.na(yl.data)] <- 0
print(dim(yl.data)) ## 21968 36100

p56.hu.data <- yl.data #[,colnames(yl.data) %in% cells2keep]
print(dim(p56.hu.data))
# 21968 36100

##------------------------------------
## Save count table
save(p56.hu.data, file = "YL_P56_HU_CB_DF_COUNTS_UPDATED.RData")

##------------------------------------
## Initialize the Seurat object with the raw (non-normalized data).
seuObj <- CreateSeuratObject(counts = p56.hu.data, project = "YL_P56_HU_CB_DF")
print(dim(seuObj)) ## 21967 36100

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

print(table(seuObj@active.ident))
# HU_P56_02 HU_P56_05 HU_P56_09 HU_P56_11 
#      9150      6189     13469      7292

##------------------------------------
## Visualize Data QC
p9 <- VlnPlot(seuObj, features = c("nFeature_RNA", "nCount_RNA", "pMito_RNA"), ncol = 3, pt.size = 0)
p2 <- FeatureScatter(seuObj, feature1 = "nCount_RNA", feature2 = "pMito_RNA")
p3 <- FeatureScatter(seuObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

ggsave(filename = "SEURAT_YL_QC_1.pdf", plot = p9, width = 8, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_YL_QC_2.pdf", plot = p2, width = 5, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_YL_QC_3.pdf", plot = p3, width = 5, height = 4, units = "in", dpi = 150)

##------------------------------------
## Data Filtering
seuObj <- subset(seuObj, subset = nCount_RNA < 10000 & pMito_RNA < 0.5)
print(dim(seuObj)) ## 21968 27628

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
all.genes <- rownames(seuObj)
seuObj <- ScaleData(seuObj, features = all.genes, vars.to.regress = c("nCount_RNA", "pMito_RNA"))

##------------------------------------
## PCA & Jackstraw
## Compute PCA
seuObj <- RunPCA(seuObj, features = VariableFeatures(object = seuObj), npcs = 100)

## Compute and Score Jackstraw
seuObj <- JackStraw(seuObj, num.replicate = 100, dims = 100)
seuObj <- ScoreJackStraw(seuObj, dims = 1:100)

## Plot PCA loadings and Jackstraw scores
p6 <- ElbowPlot(seuObj, ndims = 100)
p56 <- JackStrawPlot(seuObj, dims = 1:100)

ggsave(filename = "SEURAT_YL_PCA_1.pdf", plot = p6, width = 6, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_YL_PCA_2.pdf", plot = p56, width = 12, height = 6, units = "in", dpi = 150)

##------------------------------------
## Data Clustering
## IDENTIFY PCS
pcScores <- seuObj@reductions$pca@jackstraw$overall.p.values
pcScoresNoSign <- pcScores[pcScores[,2] > 0.05,]
selpcs <- min(pcScoresNoSign[,1]) - 1
print(selpcs)


seuObj <- FindNeighbors(seuObj, dims = 1:selpcs)
seuObj <- FindClusters(seuObj, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 0.8, 1.4, 1.6, 1.8, 2.0))
seuObj <- RunUMAP(seuObj, dims = 1:selpcs)

##------------------------------------
## Save RData
save(seuObj, file = "YL_P56_HU_CB_DF_CLUST.RData")

##------------------------------------
## Subset
seuObj_replace_plus <- subset(x = seuObj, subset = Clock > 0 & `Clock-human` > 0 & Arntl > 0) # 4372
seuObj_replace_minus <- subset(x = seuObj, subset = Clock > 0 & `Clock-human` > 0 & Arntl <= 0) # 7177

seuObj_mis_plus <- subset(x = seuObj, subset = Clock <= 0 & `Clock-human` > 0 & Arntl > 0) # 3082
seuObj_mis_minus <- subset(x = seuObj, subset = Clock <= 0 & `Clock-human` > 0 & Arntl <= 0) # 8191

seuObj_ko_plus <- subset(x = seuObj, subset = Clock > 0 & `Clock-human` <= 0 & Arntl > 0) # 187
seuObj_ko_minus <- subset(x = seuObj, subset = Clock > 0 & `Clock-human` <= 0 & Arntl <= 0) # 985

seuObj_no_plus <- subset(x = seuObj, subset = Clock <= 0 & `Clock-human` <= 0 & Arntl > 0) # 2605
seuObj_no_minus <- subset(x = seuObj, subset = Clock <= 0 & `Clock-human` <= 0 & Arntl <= 0) # 8183

save(seuObj_replace_plus, file = "YL_P56_HU_CB_DF_replace_plus.RData")
save(seuObj_replace_minus, file = "YL_P56_HU_CB_DF_replace_minus.RData")
save(seuObj_mis_plus, file = "YL_P56_HU_CB_DF_mis_plus.RData")
save(seuObj_mis_minus, file = "YL_P56_HU_CB_DF_mis_minus.RData")
save(seuObj_ko_plus, file = "YL_P56_HU_CB_DF_ko_plus.RData")
save(seuObj_ko_minus, file = "YL_P56_HU_CB_DF_ko_minus.RData")
save(seuObj_no_plus, file = "YL_P56_HU_CB_DF_no_plus.RData")
save(seuObj_no_minus, file = "YL_P56_HU_CB_DF_no_minus.RData")

```

-----
