# SEURAT INTEGRATION, HUMAN, JORSTAD 2023

## Filter CellBender+DoubletFinder Datasets
```{R}
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
cbdf_nwtx_0024_8 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0024_8_Out_filtered.h5",
"/project/D_DOUBLETFINDER_ALLEN/NW_TX0024_8_Seurat_DF.RData",
"NW_TX0024_8")

cbdf_nwtx_0024_9 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0024_9_Out_filtered.h5",
"/project/D_DOUBLETFINDER_ALLEN/NW_TX0024_9_Seurat_DF.RData",
"NW_TX0024_9")

cbdf_nwtx_0024_12 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0024_12_Out_filtered.h5",
"/project/D_DOUBLETFINDER_ALLEN/NW_TX0024_12_Seurat_DF.RData",
"NW_TX0024_12")

cbdf_nwtx_0024_13 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_NW_TX0024_13_Out_filtered.h5",
"/project/D_DOUBLETFINDER_ALLEN/NW_TX0024_13_Seurat_DF.RData",
"NW_TX0024_13")

cbdf_10x_216_7 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_10X216_7_Out_filtered.h5",
"/project/D_DOUBLETFINDER_ALLEN/10X216_7_Seurat_DF.RData",
"10X216_7")

cbdf_10x_216_8 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_10X216_8_Out_filtered.h5",
"/project/D_DOUBLETFINDER_ALLEN/10X216_8_Seurat_DF.RData",
"10X216_8")

cbdf_10x_240_1 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_10X240_1_Out_filtered.h5",
"/project/D_DOUBLETFINDER_ALLEN/10X240_1_Seurat_DF.RData",
"10X240_1")

cbdf_10x_240_2 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_10X240_2_Out_filtered.h5",
"/project/D_DOUBLETFINDER_ALLEN/10X240_2_Seurat_DF.RData",
"10X240_2")

cbdf_10x_240_3 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_10X240_3_Out_filtered.h5",
"/project/D_DOUBLETFINDER_ALLEN/10X240_3_Seurat_DF.RData",
"10X240_3")

cbdf_10x_240_4 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_10X240_4_Out_filtered.h5",
"/project/D_DOUBLETFINDER_ALLEN/10X240_4_Seurat_DF.RData",
"10X240_4")

cbdf_10x_240_5 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_10X240_5_Out_filtered.h5",
"/project/D_DOUBLETFINDER_ALLEN/10X240_5_Seurat_DF.RData",
"10X240_5")

cbdf_10x_240_6 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_10X240_6_Out_filtered.h5",
"/project/D_DOUBLETFINDER_ALLEN/10X240_6_Seurat_DF.RData",
"10X240_6")

cbdf_10x_240_7 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_10X240_7_Out_filtered.h5",
"/project/D_DOUBLETFINDER_ALLEN/10X240_7_Seurat_DF.RData",
"10X240_7")

cbdf_10x_240_8 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_10X240_8_Out_filtered.h5",
"/project/D_DOUBLETFINDER_ALLEN/10X240_8_Seurat_DF.RData",
"10X240_8")

cbdf_10x_241_1 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_10X241_1_Out_filtered.h5",
"/project/D_DOUBLETFINDER_ALLEN/10X241_1_Seurat_DF.RData",
"10X241_1")

cbdf_10x_241_2 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_10X241_2_Out_filtered.h5",
"/project/D_DOUBLETFINDER_ALLEN/10X241_2_Seurat_DF.RData",
"10X241_2")

cbdf_10x_241_3 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_10X241_3_Out_filtered.h5",
"/project/D_DOUBLETFINDER_ALLEN/10X241_3_Seurat_DF.RData",
"10X241_3")

cbdf_10x_241_4 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_10X241_4_Out_filtered.h5",
"/project/D_DOUBLETFINDER_ALLEN/10X241_4_Seurat_DF.RData",
"10X241_4")

cbdf_10x_241_5 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_10X241_5_Out_filtered.h5",
"/project/D_DOUBLETFINDER_ALLEN/10X241_5_Seurat_DF.RData",
"10X241_5")

cbdf_10x_241_6 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_10X241_6_Out_filtered.h5",
"/project/D_DOUBLETFINDER_ALLEN/10X241_6_Seurat_DF.RData",
"10X241_6")

cbdf_10x_241_7 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_10X241_7_Out_filtered.h5",
"/project/D_DOUBLETFINDER_ALLEN/10X241_7_Seurat_DF.RData",
"10X241_7")

cbdf_10x_241_8 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_10X241_8_Out_filtered.h5",
"/project/D_DOUBLETFINDER_ALLEN/10X241_8_Seurat_DF.RData",
"10X241_8")

cbdf_10x_359_4 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_10X359_4_Out_filtered.h5",
"/project/D_DOUBLETFINDER_ALLEN/10X359_4_Seurat_DF.RData",
"10X359_4")

cbdf_10x_359_5 <- filterCounts("/project/C_CELLBENDER_ALLEN/CELLBENDER_10X359_5_Out_filtered.h5",
"/project/D_DOUBLETFINDER_ALLEN/10X359_5_Seurat_DF.RData",
"10X359_5")

save(cbdf_nwtx_0024_8, cbdf_nwtx_0024_9, cbdf_nwtx_0024_12, cbdf_nwtx_0024_13,
cbdf_10x_216_8,
cbdf_10x_240_1, cbdf_10x_240_2, cbdf_10x_240_3, cbdf_10x_240_4, cbdf_10x_240_5, cbdf_10x_240_6, cbdf_10x_240_7, cbdf_10x_240_8,
cbdf_10x_241_1, cbdf_10x_241_2, cbdf_10x_241_3, cbdf_10x_241_4, cbdf_10x_241_5, cbdf_10x_241_6, cbdf_10x_241_7, cbdf_10x_241_8,
cbdf_10x_359_4, cbdf_10x_359_5,
file = "ALLEN_HUMAN_CBDF_COUNTS.RData")
```

## Seurat Preprocessing
```{R}

## LOAD CB COUNT TABLES
rm(list = ls())
load("ALLEN_HUMAN_CBDF_COUNTS.RData")
# cbdf_nwtx_0024_8, cbdf_nwtx_0024_9, cbdf_nwtx_0024_12, cbdf_nwtx_0024_13,
# cbdf_10x_216_8,
# cbdf_10x_240_1, cbdf_10x_240_2, cbdf_10x_240_3, cbdf_10x_240_4, cbdf_10x_240_5, cbdf_10x_240_6, cbdf_10x_240_7, cbdf_10x_240_8,
# cbdf_10x_241_1, cbdf_10x_241_2, cbdf_10x_241_3, cbdf_10x_241_4, cbdf_10x_241_5, cbdf_10x_241_6, cbdf_10x_241_7, cbdf_10x_241_8,
# cbdf_10x_359_4, cbdf_10x_359_5,

## FUNCTION TO RUN SEURAT PREPROCESSING PER SAMPLE
preprocess_seurat <- function(countsTable, samplePrefix, sampleSpecies) {
print(samplePrefix)
seuobj <- CreateSeuratObject(counts = countsTable)
seuobj <- PercentageFeatureSet(seuobj, pattern = "^MT-", col.name = "pMito_RNA")
seuobj$Sample <- samplePrefix
seuobj$Species <- sampleSpecies
return(seuobj)
}

seu_nwtx_0024_8 <- preprocess_seurat(cbdf_nwtx_0024_8, "NW_TX0024_8", "HUMAN")
seu_nwtx_0024_9 <- preprocess_seurat(cbdf_nwtx_0024_9, "NW_TX0024_9", "HUMAN")
seu_nwtx_0024_12 <- preprocess_seurat(cbdf_nwtx_0024_12, "NW_TX0024_12", "HUMAN")
seu_nwtx_0024_13 <- preprocess_seurat(cbdf_nwtx_0024_13, "NW_TX0024_13", "HUMAN")
seu_10x_216_8 <- preprocess_seurat(cbdf_10x_216_8, "10X216_8", "HUMAN")
seu_10x_240_1 <- preprocess_seurat(cbdf_10x_240_1, "10X240_1", "HUMAN")
seu_10x_240_2 <- preprocess_seurat(cbdf_10x_240_2, "10X240_2", "HUMAN")
seu_10x_240_3 <- preprocess_seurat(cbdf_10x_240_3, "10X240_3", "HUMAN")
seu_10x_240_4 <- preprocess_seurat(cbdf_10x_240_4, "10X240_4", "HUMAN")
seu_10x_240_5 <- preprocess_seurat(cbdf_10x_240_5, "10X240_5", "HUMAN")
seu_10x_240_6 <- preprocess_seurat(cbdf_10x_240_6, "10X240_6", "HUMAN")
seu_10x_240_7 <- preprocess_seurat(cbdf_10x_240_7, "10X240_7", "HUMAN")
seu_10x_240_8 <- preprocess_seurat(cbdf_10x_240_8, "10X240_8", "HUMAN")
seu_10x_241_1 <- preprocess_seurat(cbdf_10x_241_1, "10X241_1", "HUMAN")
seu_10x_241_2 <- preprocess_seurat(cbdf_10x_241_2, "10X241_2", "HUMAN")
seu_10x_241_3 <- preprocess_seurat(cbdf_10x_241_3, "10X241_3", "HUMAN")
seu_10x_241_4 <- preprocess_seurat(cbdf_10x_241_4, "10X241_4", "HUMAN")
seu_10x_241_5 <- preprocess_seurat(cbdf_10x_241_5, "10X241_5", "HUMAN")
seu_10x_241_6 <- preprocess_seurat(cbdf_10x_241_6, "10X241_6", "HUMAN")
seu_10x_241_7 <- preprocess_seurat(cbdf_10x_241_7, "10X241_7", "HUMAN")
seu_10x_241_8 <- preprocess_seurat(cbdf_10x_241_8, "10X241_8", "HUMAN")
seu_10x_359_4 <- preprocess_seurat(cbdf_10x_359_4, "10X359_4", "HUMAN")
seu_10x_359_5 <- preprocess_seurat(cbdf_10x_359_5, "10X359_5", "HUMAN")

seu_list <- list("NW_TX0024_8" = seu_nwtx_0024_8,
"NW_TX0024_9" = seu_nwtx_0024_9,
"NW_TX0024_12" = seu_nwtx_0024_12,
"NW_TX0024_13" = seu_nwtx_0024_13,
"NM_10X216_8" = seu_10x_216_8,
"NM_10X240_1" = seu_10x_240_1,
"NM_10X240_2" = seu_10x_240_2,
"NM_10X240_3" = seu_10x_240_3,
"NM_10X240_4" = seu_10x_240_4,
"NM_10X240_5" = seu_10x_240_5,
"NM_10X240_6" = seu_10x_240_6,
"NM_10X240_7" = seu_10x_240_7,
"NM_10X240_8" = seu_10x_240_8,
"NM_10X241_1" = seu_10x_241_1,
"NM_10X241_2" = seu_10x_241_2,
"NM_10X241_3" = seu_10x_241_3,
"NM_10X241_4" = seu_10x_241_4,
"NM_10X241_5" = seu_10x_241_5,
"NM_10X241_6" = seu_10x_241_6,
"NM_10X241_7" = seu_10x_241_7,
"NM_10X241_8" = seu_10x_241_8,
"NM_10X359_4" = seu_10x_359_4,
"NM_10X359_5" = seu_10x_359_5
)

save(seu_list, file = "ALLEN_HUMAN_CBDF_SEURAT.RData")
```

## Merge & Integration
```{R}
## load list of seurat objects
rm(list = ls())
load("ALLEN_HUMAN_CBDF_SEURAT.RData")
names(seu_list)

## merge list of seurat objects into a single seurat object
seuObjMerged <- merge(x = seu_list$NW_TX0024_8,
y = c(seu_list$NW_TX0024_9,
seu_list$NW_TX0024_12,
seu_list$NW_TX0024_13,
seu_list$NM_10X216_8,
seu_list$NM_10X240_1,
seu_list$NM_10X240_2,
seu_list$NM_10X240_3,
seu_list$NM_10X240_4,
seu_list$NM_10X240_5,
seu_list$NM_10X240_6,
seu_list$NM_10X240_7,
seu_list$NM_10X240_8,
seu_list$NM_10X241_1,
seu_list$NM_10X241_2,
seu_list$NM_10X241_3,
seu_list$NM_10X241_4,
seu_list$NM_10X241_5,
seu_list$NM_10X241_6,
seu_list$NM_10X241_7,
seu_list$NM_10X241_8,
seu_list$NM_10X359_4,
seu_list$NM_10X359_5
),
project = "ALLEN_HUMAN")

## run standard anlaysis workflow
seuObjMerged <- NormalizeData(seuObjMerged)
seuObjMerged <- FindVariableFeatures(seuObjMerged)
seuObjMerged <- ScaleData(seuObjMerged)
seuObjMerged <- RunPCA(seuObjMerged)
seuObjMerged <- FindNeighbors(seuObjMerged, dims = 1:30, reduction = "pca")
seuObjMerged <- FindClusters(seuObjMerged, resolution = 0.8, cluster.name = "merged_clusters")

## save merged object
save(seuObjMerged, file = "ALLEN_HUMAN_CBDF_SEURAT_MERGED.RData")

## perform integration
seuObjInt <- IntegrateLayers(object = seuObjMerged, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", verbose = FALSE)

# re-join layers after integration
seuObjInt[["RNA"]] <- JoinLayers(seuObjInt[["RNA"]])
seuObjInt <- FindNeighbors(seuObjInt, reduction = "integrated.cca", dims = 1:30)
seuObjInt <- FindClusters(seuObjInt, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0), n.iter = 1000)
seuObjInt <- RunUMAP(seuObjInt, dims = 1:30, reduction = "integrated.cca", min.dist = 0.1, spread = 1, n.epochs = 1000)

## Save integrated object
save(seuObjInt, file = "ALLEN_HUMAN_CBDF_SEURAT_INTEGRATED.RData")
```

## Generate Shiny App
```{R}
## Generate Shiny App
rm(list = ls())
load("ALLEN_HUMAN_CBDF_SEURAT_INTEGRATED.RData")
# seuObjInt
seuObjInt[["RNA3"]] <- as(object = seuObjInt[["RNA"]], Class = "Assay")
shinyConf = createConfig(seuObjInt)
makeShinyApp(seuObjInt, shinyConf, gene.mapping = TRUE, gex.assay = "RNA3", gex.slot = "data", shiny.dir = "ALLEN_HUMAN/", shiny.title = "ALLEN_HUMAN")

## END
```

## Data Filtering
```{R}
## Load Integrated Seurat Object
rm(list = ls())
load("ALLEN_HUMAN_CBDF_SEURAT_INTEGRATED.RData")

DefaultAssay(seuObjInt) <- "RNA"
Idents(seuObjInt) <- "RNA_snn_res.0.8"

## filtering criterion
nuhi <- 10000
pmhi <- 5

print(dim(seuObjInt))
# 20046 174340

filt_temp <- subset(seuObjInt, subset = nCount_RNA < nuhi & pMito_RNA < pmhi)
print(dim(filt_temp))
# 20046 124594

seuObjIntFilt <- subset(filt_temp, idents = c(5, 6, 7, 8, 9, 13, 14, 15, 16, 17, 18, 19, 20, 21, 25, 27), invert = TRUE) ## post filtering (earlier step, few cluster have no cells left, hence removed from this filtering)
print(dim(seuObjIntFilt))
# 20046 109813

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
ggsave(filename = "ALLEN_HUMAN_CBDF_SEURAT_INTEGRATED_FILT_CLUSTREE.pdf", plot = seutree, width = 24, height = 24, units = "in", dpi = 150)

## Save integrated object
save(seuObjIntFilt, file = "ALLEN_HUMAN_CBDF_SEURAT_INTEGRATED_FILT.RData")
```


## Generate Shiny App
```{R}
rm(list = ls())
load("ALLEN_HUMAN_CBDF_SEURAT_INTEGRATED_FILT.RData")
# seuObjIntFilt

seuObjIntFilt[["RNA3"]] <- as(object = seuObjIntFilt[["RNA"]], Class = "Assay")
shinyConf = createConfig(seuObjIntFilt)
makeShinyApp(seuObjIntFilt, shinyConf, gene.mapping = TRUE, gex.assay = "RNA3", gex.slot = "data", shiny.dir = "ALLEN_HUMAN_FILT/", shiny.title = "ALLEN_HUMAN_FILT")
```


## Additional Data Filtering
```{R}
rm(list = ls())
load("ALLEN_HUMAN_CBDF_SEURAT_INTEGRATED_FILT.RData")

DefaultAssay(seuObjIntFilt) <- "RNA"
Idents(seuObjIntFilt) <- "RNA_snn_res.0.8"

## filtering criterion
nuhi <- 10000
pmhi <- 5

print(dim(seuObjIntFilt))
# 20046 174340

filt_temp <- seuObjIntFilt
rm(seuObjIntFilt)

seuObjIntFilt <- subset(filt_temp, idents = c(14, 28, 30, 33, 34, 35, 36, 37, 39, 40, 41), invert = TRUE) ## post filtering (earlier step, few cluster have no cells left, hence removed from this filtering)
print(dim(seuObjIntFilt))
# 20046 109813

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
ggsave(filename = "ALLEN_HUMAN_CBDF_SEURAT_INTEGRATED_FILT_CLUSTREE.pdf", plot = seutree, width = 24, height = 24, units = "in", dpi = 150)

## Save integrated object
save(seuObjIntFilt, file = "ALLEN_HUMAN_CBDF_SEURAT_INTEGRATED_FILT2.RData")
```


## Generate Shiny App
```{R}
rm(list = ls())
load("ALLEN_HUMAN_CBDF_SEURAT_INTEGRATED_FILT2.RData")
# seuObjIntFilt

seuObjIntFilt[["RNA3"]] <- as(object = seuObjIntFilt[["RNA"]], Class = "Assay")
shinyConf = createConfig(seuObjIntFilt)
makeShinyApp(seuObjIntFilt, shinyConf, gene.mapping = TRUE, gex.assay = "RNA3", gex.slot = "data", shiny.dir = "ALLEN_HUMAN_FILT2/", shiny.title = "ALLEN_HUMAN_FILT2")

## END
```
