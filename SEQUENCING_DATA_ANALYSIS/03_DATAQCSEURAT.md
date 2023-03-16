# DATA QC BEFORE FILTERING (AMBIENT RNA & DOUBLETS REMOVAL)
- python/3.7.x-anaconda
- gcc/8.3.0
- hdf5_18/1.8.17
- R/4.0.2-gccmkl


```{r}
##-------------------------------------------------------
## SEURAT QC ANALYSIS
##-------------------------------------------------------
rm(list = ls())
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(ggrepel)
library(reticulate)
library(RColorBrewer)

##-------------------------------------------------------
## MAIN
##-------------------------------------------------------
print("-----> READING INPUT DATA")
# COLLAPSED AND COMBINED DATA COPIED FROM FOLLOWING LOCATION
load("/work/Neuroinformatics_Core/akulk1/MANUAL/YL_DATA_COMBINED/01_COMBINE_COUNTS/YL_P07_P56_MANUAL_COUNTS_COMBINED.RData")
# allData

##-------------------------------------------------------
## CREATE SEURAT OBJECT
##-------------------------------------------------------
print("-----> CREATING SEURAT OBJECT")
seuObjAll <- CreateSeuratObject(counts = allData, project = "YL")

##-------------------------------------------------------
## CREATE AND UPDATE META DATA
##-------------------------------------------------------
print("-----> UPDATING META DATA")
## CREATE META DATA
metaTemp1 <- as.data.frame(seuObjAll@meta.data)
colnames(metaTemp1) <- c("Dataset", "nUMI", "nGenes")

metaTemp2 <- as.data.frame(matrix(unlist(strsplit(row.names(metaTemp1), "_")), byrow = TRUE, ncol = 5))
row.names(metaTemp2) <- row.names(metaTemp1)
colnames(metaTemp2) <- c("Dataset", "Age", "Genotype", "Library", "Barcode")
metaTemp2$Index <- paste(metaTemp2$Dataset, metaTemp2$Age, metaTemp2$Genotype, metaTemp2$Library, sep = "_")
metaTemp2$CellBarcode <- row.names(metaTemp2)

# metaTemp3 <- metaTemp2
# metaTemp3$Barcode <- NULL
# metaTemp3$Index <- paste(metaTemp3$Dataset, metaTemp3$Age, metaTemp3$Genotype, metaTemp3$Library, sep = "_")
# metaTemp4 <- as.data.frame(matrix(unlist(strsplit(unique(metaTemp3$Index), "_")), byrow = TRUE, ncol = 4))
# row.names(metaTemp4) <- unique(metaTemp3$Index)
# colnames(metaTemp4) <- c("Dataset", "Age", "Genotype", "Library")
# write.table(metaTemp4, "YL_P07_P56_META_TEMP.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

metaTemp3 <- read.table("YL_P07_P56_META.txt", header = TRUE, sep = "\t")

metaTemp <- merge(metaTemp2[,c("Index", "CellBarcode")], metaTemp3, by = "Index", all.x = TRUE)
row.names(metaTemp) <- metaTemp$CellBarcode


# "Index", "CellBarcode", "Dataset", "Age", "Genotype", "Library", "Sex", "RNAPrep", "LibraryPrep", "Sequencing" 

## UPDATE SEURAT OBJECT
seuObjAll[["pMito_RNA"]] <- PercentageFeatureSet(seuObjAll, pattern = "^mt-")

metaId <- metaTemp$Index
names(metaId) <- row.names(metaTemp)
seuObjAll <- AddMetaData(object = seuObjAll, metadata = metaId, col.name = "SampleIndex")

metaCB <- metaTemp$CellBarcode
names(metaCB) <- row.names(metaTemp)
seuObjAll <- AddMetaData(object = seuObjAll, metadata = metaCB, col.name = "CellBarcode")

metaAge <- metaTemp$Age
names(metaAge) <- row.names(metaTemp)
seuObjAll <- AddMetaData(object = seuObjAll, metadata = metaAge, col.name = "Age")

metaAge <- metaTemp$Age
names(metaAge) <- row.names(metaTemp)
seuObjAll <- AddMetaData(object = seuObjAll, metadata = metaAge, col.name = "Age")

metaGeno <- metaTemp$Genotype
names(metaGeno) <- row.names(metaTemp)
seuObjAll <- AddMetaData(object = seuObjAll, metadata = metaGeno, col.name = "Genotype")

metaSex <- metaTemp$Sex
names(metaSex) <- row.names(metaTemp)
seuObjAll <- AddMetaData(object = seuObjAll, metadata = metaSex, col.name = "Sex")

metaRNA <- metaTemp$RNAPrep
names(metaRNA) <- row.names(metaTemp)
seuObjAll <- AddMetaData(object = seuObjAll, metadata = metaRNA, col.name = "RNAPrep")

metaLib <- metaTemp$LibraryPrep
names(metaLib) <- row.names(metaTemp)
seuObjAll <- AddMetaData(object = seuObjAll, metadata = metaLib, col.name = "LibPrep")

metaSeq <- metaTemp$Sequencing
names(metaSeq) <- row.names(metaTemp)
seuObjAll <- AddMetaData(object = seuObjAll, metadata = metaSeq, col.name = "SeqRun")


##-------------------------------------------------------
## SAVE QC DATA
##-------------------------------------------------------
print("-----> SAVING SEURAT OBJECT WITH UPDATED META DATA")
save(seuObjAll, file = "YL_P07_P56_MANUAL_DATA_QC.RData")

metaData <- as.data.frame(seuObjAll@meta.data)
save(metaData, file = "YL_P07_P56_MANUAL_META_QC.RData")


##-------------------------------------------------------
## QC PLOTS
##-------------------------------------------------------
print("-----> GENERATING QC PLOTS BEFORE FILTERING")

for(myfactor in c("SampleIndex", "Age", "Genotype", "Sex", "RNAPrep", "LibPrep", "SeqRun"))
	{
	print(myfactor)
	## QC Plots before removing mito genes
	## violin plots
	plot_nU <- VlnPlot(object = seuObjAll, features = "nCount_RNA", pt.size = 0, group.by = myfactor)
	plot_nG <- VlnPlot(object = seuObjAll, features = "nFeature_RNA", pt.size = 0, group.by = myfactor)
	plot_pM <- VlnPlot(object = seuObjAll, features = "pMito_RNA", pt.size = 0, group.by = myfactor)
	plotQC1 <- grid.arrange(plot_nU, plot_nG, plot_pM, ncol = 3)
	ggsave(paste(seuObjAll@project.name, "_QC_", myfactor, "_1.pdf", sep = ""), plot = plotQC1, width = 12, height = 4, units = "in", dpi = 300)
	## scatter plots
	plot_nUpM <- FeatureScatter(object = seuObjAll, feature1 = "nCount_RNA", feature2 = "pMito_RNA", pt.size = 0.1, group.by = myfactor)
	plot_nUnG <- FeatureScatter(object = seuObjAll, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.1, group.by = myfactor)
	plotQC2 <- grid.arrange(plot_nUpM, plot_nUnG, ncol = 2)
	ggsave(paste(seuObjAll@project.name, "_QC_", myfactor, "_2.pdf", sep = ""), plot = plotQC2, width = 10, height = 4, units = "in", dpi = 300)
	}


print("-----> FINISHED")
##-------------------------------------------------------
## END
##-------------------------------------------------------
```
