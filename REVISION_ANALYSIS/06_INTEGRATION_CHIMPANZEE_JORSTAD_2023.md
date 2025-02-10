# SEURAT INTEGRATION, CHIMPANZEE, JORSTAD 2023

## Filter CellBender+DoubletFinder Datasets
```{r}
## START
## R Libraries
rm(list = ls())
.libPaths("/project/RESOURCES/FROM_HOME/x86_64-pc-linux-gnu-library/4.1")
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(reshape2)
library(clustree)
library(dittoSeq)
library(scCustomize)
library(ShinyCell)
set.seed(10)

## FUNCTION TO FILTER CELLBARCODES (CB-RETAINED) AND GENES (REF-PROTEIN-CODING)
filterCounts <- function(cbPath, dfPath, samplePrefix) {
print(samplePrefix)

## load cb data
cbData <- Read_CellBender_h5_Mat(cbPath, use.names = TRUE, unique.features = TRUE, h5_group_name = NULL, feature_slot_name = "features")
print(dim(cbData))
colnames(cbData) <- paste(samplePrefix, colnames(cbData), sep = "_")
cbData <- as.data.frame(cbData)

## load df data
load(dfPath)
dfRetained <- row.names(seuObj_doublets_removed@meta.data)

## filter cb data
cbdfData <- cbData[, colnames(cbData) %in% dfRetained]
print(dim(cbdfData))

## return
return(cbdfData)
}

## LOAD CELLBENDER H5 FILE
cbdf_nwtx_0017_12 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0017_12_Out_filtered.h5",
"/project/D_DOUBLETFINDER_ALLEN/NW_TX0017_12_Seurat_DF.RData",
"NW_TX0017_12")

cbdf_nwtx_0017_13 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0017_13_Out_filtered.h5",
"/project/D_DOUBLETFINDER_ALLEN/NW_TX0017_13_Seurat_DF.RData",
"NW_TX0017_13")

cbdf_nwtx_0020_6 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0020_6_Out_filtered.h5",
"/project/D_DOUBLETFINDER_ALLEN/NW_TX0020_6_Seurat_DF.RData",
"NW_TX0020_6")

cbdf_nwtx_0020_7 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0020_7_Out_filtered.h5",
"/project/D_DOUBLETFINDER_ALLEN/NW_TX0020_7_Seurat_DF.RData",
"NW_TX0020_7")

cbdf_nwtx_0020_10 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0020_10_Out_filtered.h5",
"/project/D_DOUBLETFINDER_ALLEN/NW_TX0020_10_Seurat_DF.RData",
"NW_TX0020_10")

cbdf_nwtx_0020_13 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0020_13_Out_filtered.h5",
"/project/D_DOUBLETFINDER_ALLEN/NW_TX0020_13_Seurat_DF.RData",
"NW_TX0020_13")

cbdf_nwtx_0025_10 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0025_10_Out_filtered.h5",
"/project/D_DOUBLETFINDER_ALLEN/NW_TX0025_10_Seurat_DF.RData",
"NW_TX0025_10")

cbdf_nwtx_0025_11 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0025_11_Out_filtered.h5",
"/project/D_DOUBLETFINDER_ALLEN/NW_TX0025_11_Seurat_DF.RData",
"NW_TX0025_11")

cbdf_nwtx_0025_12 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0025_12_Out_filtered.h5",
"/project/D_DOUBLETFINDER_ALLEN/NW_TX0025_12_Seurat_DF.RData",
"NW_TX0025_12")

cbdf_nwtx_0027_6 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0027_6_Out_filtered.h5",
"/project/D_DOUBLETFINDER_ALLEN/NW_TX0027_6_Seurat_DF.RData",
"NW_TX0027_6")

cbdf_nwtx_0027_7 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0027_7_Out_filtered.h5",
"/project/D_DOUBLETFINDER_ALLEN/NW_TX0027_7_Seurat_DF.RData",
"NW_TX0027_7")

cbdf_nwtx_0027_8 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0027_8_Out_filtered.h5",
"/project/D_DOUBLETFINDER_ALLEN/NW_TX0027_8_Seurat_DF.RData",
"NW_TX0027_8")

cbdf_nwtx_0027_9 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0027_9_Out_filtered.h5",
"/project/D_DOUBLETFINDER_ALLEN/NW_TX0027_9_Seurat_DF.RData",
"NW_TX0027_9")

cbdf_nwtx_0028_1 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0028_1_Out_filtered.h5",
"/project/D_DOUBLETFINDER_ALLEN/NW_TX0028_1_Seurat_DF.RData",
"NW_TX0028_1")

cbdf_nwtx_0053_8 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0053_8_Out_filtered.h5",
"/project/D_DOUBLETFINDER_ALLEN/NW_TX0053_8_Seurat_DF.RData",
"NW_TX0053_8")

cbdf_nwtx_0053_9 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0053_9_Out_filtered.h5",
"/project/D_DOUBLETFINDER_ALLEN/NW_TX0053_9_Seurat_DF.RData",
"NW_TX0053_9")

cbdf_nwtx_0055_7 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0055_7_Out_filtered.h5",
"/project/D_DOUBLETFINDER_ALLEN/NW_TX0055_7_Seurat_DF.RData",
"NW_TX0055_7")

cbdf_nwtx_0055_8 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0055_8_Out_filtered.h5",
"/project/D_DOUBLETFINDER_ALLEN/NW_TX0055_8_Seurat_DF.RData",
"NW_TX0055_8")

cbdf_nwtx_0055_9 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0055_9_Out_filtered.h5",
"/project/D_DOUBLETFINDER_ALLEN/NW_TX0055_9_Seurat_DF.RData",
"NW_TX0055_9")

cbdf_nwtx_0067_8 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0067_8_Out_filtered.h5",
"/project/D_DOUBLETFINDER_ALLEN/NW_TX0067_8_Seurat_DF.RData",
"NW_TX0067_8")

cbdf_nwtx_0067_9 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0067_9_Out_filtered.h5",
"/project/D_DOUBLETFINDER_ALLEN/NW_TX0067_9_Seurat_DF.RData",
"NW_TX0067_9")

cbdf_nwtx_0067_10 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0067_10_Out_filtered.h5",
"/project/D_DOUBLETFINDER_ALLEN/NW_TX0067_10_Seurat_DF.RData",
"NW_TX0067_10")

cbdf_nwtx_0068_8 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0068_8_Out_filtered.h5",
"/project/D_DOUBLETFINDER_ALLEN/NW_TX0068_8_Seurat_DF.RData",
"NW_TX0068_8")

cbdf_nwtx_0068_9 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0068_9_Out_filtered.h5",
"/project/D_DOUBLETFINDER_ALLEN/NW_TX0068_9_Seurat_DF.RData",
"NW_TX0068_9")

cbdf_nwtx_0070_1 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0070_1_Out_filtered.h5",
"/project/D_DOUBLETFINDER_ALLEN/NW_TX0070_1_Seurat_DF.RData",
"NW_TX0070_1")

cbdf_nwtx_0070_6 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0070_6_Out_filtered.h5",
"/project/D_DOUBLETFINDER_ALLEN/NW_TX0070_6_Seurat_DF.RData",
"NW_TX0070_6")

save(cbdf_nwtx_0017_12, cbdf_nwtx_0017_13, cbdf_nwtx_0020_6, cbdf_nwtx_0020_7, cbdf_nwtx_0020_10, cbdf_nwtx_0020_13, cbdf_nwtx_0025_10, cbdf_nwtx_0025_11, cbdf_nwtx_0025_12, cbdf_nwtx_0027_6, cbdf_nwtx_0027_7, cbdf_nwtx_0027_8, cbdf_nwtx_0027_9, cbdf_nwtx_0028_1, cbdf_nwtx_0053_8, cbdf_nwtx_0053_9, cbdf_nwtx_0055_7, cbdf_nwtx_0055_8, cbdf_nwtx_0055_9, cbdf_nwtx_0067_8, cbdf_nwtx_0067_9, cbdf_nwtx_0067_10, cbdf_nwtx_0068_8, cbdf_nwtx_0068_9, cbdf_nwtx_0070_1, cbdf_nwtx_0070_6,
file = "ALLEN_CHIMPANZEE_CBDF_COUNTS.RData")
```

## Seurat Preprocessing
```{r}
## LOAD CB COUNT TABLES
rm(list = ls())
load("ALLEN_CHIMPANZEE_CBDF_COUNTS.RData")

## FUNCTION TO RUN SEURAT PREPROCESSING PER SAMPLE
preprocess_seurat <- function(countsTable, samplePrefix, sampleSpecies) {
print(samplePrefix)
seuobj <- CreateSeuratObject(counts = countsTable)
seuobj <- PercentageFeatureSet(seuobj, pattern = "^MT-", col.name = "pMito_RNA")
seuobj$Sample <- samplePrefix
seuobj$Species <- sampleSpecies
return(seuobj)
}

seu_nwtx_0017_12 <- preprocess_seurat(cbdf_nwtx_0017_12, "NW_TX0017_12", "CHIMPANZEE")
seu_nwtx_0017_13 <- preprocess_seurat(cbdf_nwtx_0017_13, "NW_TX0017_13", "CHIMPANZEE")
seu_nwtx_0020_6 <- preprocess_seurat(cbdf_nwtx_0020_6, "NW_TX0020_6", "CHIMPANZEE")
seu_nwtx_0020_7 <- preprocess_seurat(cbdf_nwtx_0020_7, "NW_TX0020_7", "CHIMPANZEE")
seu_nwtx_0020_10 <- preprocess_seurat(cbdf_nwtx_0020_10, "NW_TX0020_10", "CHIMPANZEE")
seu_nwtx_0020_13 <- preprocess_seurat(cbdf_nwtx_0020_13, "NW_TX0020_13", "CHIMPANZEE")
seu_nwtx_0025_10 <- preprocess_seurat(cbdf_nwtx_0025_10, "NW_TX0025_10", "CHIMPANZEE")
seu_nwtx_0025_11 <- preprocess_seurat(cbdf_nwtx_0025_11, "NW_TX0025_11", "CHIMPANZEE")
seu_nwtx_0025_12 <- preprocess_seurat(cbdf_nwtx_0025_12, "NW_TX0025_12", "CHIMPANZEE")
seu_nwtx_0027_6 <- preprocess_seurat(cbdf_nwtx_0027_6, "NW_TX0027_6", "CHIMPANZEE")
seu_nwtx_0027_7 <- preprocess_seurat(cbdf_nwtx_0027_7, "NW_TX0027_7", "CHIMPANZEE")
seu_nwtx_0027_8 <- preprocess_seurat(cbdf_nwtx_0027_8, "NW_TX0027_8", "CHIMPANZEE")
seu_nwtx_0027_9 <- preprocess_seurat(cbdf_nwtx_0027_9, "NW_TX0027_9", "CHIMPANZEE")
seu_nwtx_0028_1 <- preprocess_seurat(cbdf_nwtx_0028_1, "NW_TX0028_1", "CHIMPANZEE")
seu_nwtx_0053_8 <- preprocess_seurat(cbdf_nwtx_0053_8, "NW_TX0053_8", "CHIMPANZEE")
seu_nwtx_0053_9 <- preprocess_seurat(cbdf_nwtx_0053_9, "NW_TX0053_9", "CHIMPANZEE")
seu_nwtx_0055_7 <- preprocess_seurat(cbdf_nwtx_0055_7, "NW_TX0055_7", "CHIMPANZEE")
seu_nwtx_0055_8 <- preprocess_seurat(cbdf_nwtx_0055_8, "NW_TX0055_8", "CHIMPANZEE")
seu_nwtx_0055_9 <- preprocess_seurat(cbdf_nwtx_0055_9, "NW_TX0055_9", "CHIMPANZEE")
seu_nwtx_0067_8 <- preprocess_seurat(cbdf_nwtx_0067_8, "NW_TX0067_8", "CHIMPANZEE")
seu_nwtx_0067_9 <- preprocess_seurat(cbdf_nwtx_0067_9, "NW_TX0067_9", "CHIMPANZEE")
seu_nwtx_0067_10 <- preprocess_seurat(cbdf_nwtx_0067_10, "NW_TX0067_10", "CHIMPANZEE")
seu_nwtx_0068_8 <- preprocess_seurat(cbdf_nwtx_0068_8, "NW_TX0068_8", "CHIMPANZEE")
seu_nwtx_0068_9 <- preprocess_seurat(cbdf_nwtx_0068_9, "NW_TX0068_9", "CHIMPANZEE")
seu_nwtx_0070_1 <- preprocess_seurat(cbdf_nwtx_0070_1, "NW_TX0070_1", "CHIMPANZEE")
seu_nwtx_0070_6 <- preprocess_seurat(cbdf_nwtx_0070_6, "NW_TX0070_6", "CHIMPANZEE")

seu_list <- list("NW_TX0017_12" = seu_nwtx_0017_12, "NW_TX0017_13" = seu_nwtx_0017_13, "NW_TX0020_6" = seu_nwtx_0020_6, "NW_TX0020_7" = seu_nwtx_0020_7,
"NW_TX0020_10" = seu_nwtx_0020_10, "NW_TX0020_13" = seu_nwtx_0020_13, "NW_TX0025_10" = seu_nwtx_0025_10, "NW_TX0025_11" = seu_nwtx_0025_11,
"NW_TX0025_12" = seu_nwtx_0025_12, "NW_TX0027_6" = seu_nwtx_0027_6, "NW_TX0027_7" = seu_nwtx_0027_7, "NW_TX0027_8" = seu_nwtx_0027_8,
"NW_TX0027_9" = seu_nwtx_0027_9, "NW_TX0028_1" = seu_nwtx_0028_1, "NW_TX0053_8" = seu_nwtx_0053_8, "NW_TX0053_9" = seu_nwtx_0053_9,
"NW_TX0055_7" = seu_nwtx_0055_7, "NW_TX0055_8" = seu_nwtx_0055_8, "NW_TX0055_9" = seu_nwtx_0055_9, "NW_TX0067_8" = seu_nwtx_0067_8,
"NW_TX0067_9" = seu_nwtx_0067_9, "NW_TX0067_10" = seu_nwtx_0067_10, "NW_TX0068_8" = seu_nwtx_0068_8, "NW_TX0068_9" = seu_nwtx_0068_9,
"NW_TX0070_1" = seu_nwtx_0070_1, "NW_TX0070_6" = seu_nwtx_0070_6)

save(seu_list, file = "ALLEN_CHIMPANZEE_CBDF_SEURAT.RData")
```

## Merge & Integration
```{r}
## Integration
## load list of seurat objects
rm(list = ls())
load("ALLEN_CHIMPANZEE_CBDF_SEURAT.RData")

## merge list of seurat objects into a single seurat object
seuObjMerged <- merge(x = seu_list$NW_TX0017_12,
y = c(seu_list$NW_TX0017_13,
seu_list$NW_TX0020_6,
seu_list$NW_TX0020_7,
seu_list$NW_TX0020_10,
seu_list$NW_TX0020_13,
seu_list$NW_TX0025_10,
seu_list$NW_TX0025_11,
seu_list$NW_TX0025_12,
seu_list$NW_TX0027_6,
seu_list$NW_TX0027_7,
seu_list$NW_TX0027_8,
seu_list$NW_TX0027_9,
seu_list$NW_TX0028_1,
seu_list$NW_TX0053_8,
seu_list$NW_TX0053_9,
seu_list$NW_TX0055_7,
seu_list$NW_TX0055_8,
seu_list$NW_TX0055_9,
seu_list$NW_TX0067_8,
seu_list$NW_TX0067_9,
seu_list$NW_TX0067_10,
seu_list$NW_TX0068_8,
seu_list$NW_TX0068_9,
seu_list$NW_TX0070_1,
seu_list$NW_TX0070_6
),
project = "ALLEN_CHIMPANZEE")

## run standard anlaysis workflow
seuObjMerged <- NormalizeData(seuObjMerged)
seuObjMerged <- FindVariableFeatures(seuObjMerged)
seuObjMerged <- ScaleData(seuObjMerged)
seuObjMerged <- RunPCA(seuObjMerged)
seuObjMerged <- FindNeighbors(seuObjMerged, dims = 1:30, reduction = "pca")
seuObjMerged <- FindClusters(seuObjMerged, resolution = 0.8, cluster.name = "merged_clusters")

## save merged object
save(seuObjMerged, file = "ALLEN_CHIMPANZEE_CBDF_SEURAT_MERGED.RData")

## perform integration
seuObjInt <- IntegrateLayers(object = seuObjMerged, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", verbose = FALSE)

# re-join layers after integration
seuObjInt[["RNA"]] <- JoinLayers(seuObjInt[["RNA"]])
seuObjInt <- FindNeighbors(seuObjInt, reduction = "integrated.cca", dims = 1:30)
seuObjInt <- FindClusters(seuObjInt, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0), n.iter = 1000)
seuObjInt <- RunUMAP(seuObjInt, dims = 1:30, reduction = "integrated.cca", min.dist = 0.1, spread = 1, n.epochs = 1000)

## Save integrated object
save(seuObjInt, file = "ALLEN_CHIMPANZEE_CBDF_SEURAT_INTEGRATED.RData")
```

## Generate Shiny App
```{r}
## Generate Shiny App
rm(list = ls())
load("ALLEN_CHIMPANZEE_CBDF_SEURAT_INTEGRATED.RData")
# seuObjInt

seuObjInt[["RNA3"]] <- as(object = seuObjInt[["RNA"]], Class = "Assay")
shinyConf = createConfig(seuObjInt)
makeShinyApp(seuObjInt, shinyConf, gene.mapping = TRUE, gex.assay = "RNA3", gex.slot = "data", shiny.dir = "ALLEN_CHIMPANZEE/", shiny.title = "ALLEN_CHIMPANZEE")

## END
```

## Filter Data
```{R}
## Load Integrated Seurat Object
rm(list = ls())
load("ALLEN_CHIMPANZEE_CBDF_SEURAT_INTEGRATED.RData")
# seuObjInt

DefaultAssay(seuObjInt) <- "RNA"
Idents(seuObjInt) <- "RNA_snn_res.0.8"

## filtering criterion
nuhi <- 10000
pmhi <- 5

print(dim(seuObjInt))
# 19605 206726

filt_temp <- subset(seuObjInt, subset = nCount_RNA < nuhi & pMito_RNA < pmhi)
print(dim(filt_temp))
# 19605 196423

seuObjIntFilt <- subset(filt_temp, idents = c(30, 31, 32, 33, 36, 37), invert = TRUE)
print(dim(seuObjIntFilt))
# 19605 194469

## Re-cluster
seuObjIntFilt$RNA_snn_res.0.2 <- NULL
seuObjIntFilt$RNA_snn_res.0.4 <- NULL
seuObjIntFilt$RNA_snn_res.0.6 <- NULL
seuObjIntFilt$RNA_snn_res.0.8 <- NULL
seuObjIntFilt$RNA_snn_res.1 <- NULL
seuObjIntFilt$RNA_snn_res.1.2 <- NULL
seuObjIntFilt$RNA_snn_res.1.4 <- NULL
seuObjIntFilt$RNA_snn_res.1.6 <- NULL
seuObjIntFilt$RNA_snn_res.1.8 <- NULL
seuObjIntFilt$RNA_snn_res.2 <- NULL
seuObjIntFilt$seurat_clusters <- NULL

seuObjIntFilt <- FindNeighbors(seuObjIntFilt, reduction = "integrated.cca", dims = 1:30)
seuObjIntFilt <- FindClusters(seuObjIntFilt, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0), n.iter = 1000)
seuObjIntFilt <- RunUMAP(seuObjIntFilt, reduction = "integrated.cca", dims = 1:30, min.dist = 0.1, spread = 1, n.epochs = 1000)

seutree <- clustree(seuObjIntFilt, prefix = "RNA_snn_res.", node_colour = "sc3_stability") # + scale_color_brewer(palette = "Set1") + scale_edge_color_continuous(low = "blue", high = "red")
ggsave(filename = "ALLEN_CHIMPANZEE_CBDF_SEURAT_INTEGRATED_FILT_CLUSTREE.pdf", plot = seutree, width = 24, height = 24, units = "in", dpi = 150)

## Save integrated object
save(seuObjIntFilt, file = "ALLEN_CHIMPANZEE_CBDF_SEURAT_INTEGRATED_FILT.RData")
```


## Generate Shiny App
```{R}
rm(list = ls())
load("ALLEN_CHIMPANZEE_CBDF_SEURAT_INTEGRATED_FILT.RData")
# seuObjIntFilt

seuObjIntFilt[["RNA3"]] <- as(object = seuObjIntFilt[["RNA"]], Class = "Assay")
shinyConf = createConfig(seuObjIntFilt)
makeShinyApp(seuObjIntFilt, shinyConf, gene.mapping = TRUE, gex.assay = "RNA3", gex.slot = "data", shiny.dir = "ALLEN_CHIMPANZEE_FILT/", shiny.title = "ALLEN_CHIMPANZEE_FILT")

## END
```

