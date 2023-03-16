# REMOVE AMBIENT RNA USING CELLBENDER
- python/3.7.x-anaconda
- cellbender

## CELLBENDER FOR P07 | WT
```{sh}
ls | grep -P "^YL_P07_WT" | xargs -I % -n 1 -P 48 sh -c 'echo %; cellbender remove-background --input % --output %/CellBender_Out.h5 --expected-cells 10000'
```

## CELLBENDER FOR P07 | KO
```{sh}
ls | grep -P "^YL_P07_KO" | xargs -I % -n 1 -P 48 sh -c 'echo %; cellbender remove-background --input % --output %/CellBender_Out.h5 --expected-cells 10000'
```

## CELLBENDER FOR P07 | HU
```{sh}
ls | grep -P "^YL_P07_HU" | xargs -I % -n 1 -P 48 sh -c 'echo %; cellbender remove-background --input % --output %/CellBender_Out.h5 --expected-cells 10000'
```
## CELLBENDER FOR P56 | WT
```{sh}
ls | grep -P "^YL_P56_WT" | xargs -I % -n 1 -P 48 sh -c 'echo %; cellbender remove-background --input % --output %/CellBender_Out.h5 --expected-cells 10000'
```

## CELLBENDER FOR P56 | KO
```{sh}
ls | grep -P "^YL_P56_KO" | xargs -I % -n 1 -P 48 sh -c 'echo %; cellbender remove-background --input % --output %/CellBender_Out.h5 --expected-cells 10000'
```

## CELLBENDER FOR P56 | HU
```{sh}
ls | grep -P "^YL_P56_HU" | xargs -I % -n 1 -P 48 sh -c 'echo %; cellbender remove-background --input % --output %/CellBender_Out.h5 --expected-cells 10000'
```

## FILTER AMBIENT RNA
- R/4.0.2-gccmkl
- gcc/6.1.0 
- geos/gcc/3.7.2
- gdal/gcc/2.4.1
- proj/gcc/6.0.0
- hdf5_18/1.8.11
- python/3.6.4-anaconda

```{R}
##-------------------------------------------------------
## SEURAT ANALYSIS | FILTER AMBIENT RNA
##-------------------------------------------------------
## LIBRARIES
rm(list = ls())
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(ggrepel)
library(reticulate)
library(SeuratData)
library(SeuratDisk)
library(scCustomize)


## LOAD SEURAT OBJECT
load("CLOCKSEQ_P07_P56_WT_HU_KO_QC.RData")
# dt.all

dt.updated <- UpdateSeuratObject(object = dt.all)

## CELLBENDER MANUAL
#### P07 WT
p07wt02 <- read.table("/work/Neuroinformatics_Core/akulk1/MANUAL/YL_DATA/CELLBENDER_MANUAL/P07_WT/YL_P07_WT_02_CB/CellBender_Out_cell_barcodes.csv", sep = ",") #        6990
colnames(p07wt02) <- "CellBarcode"
p07wt02$Index <- paste("YL_P07_WT_02", p07wt02$CellBarcode, sep = "_")

p07wt04 <- read.table("/work/Neuroinformatics_Core/akulk1/MANUAL/YL_DATA/CELLBENDER_MANUAL/CELLBENDER_AGAIN/YL_P07_WT_04_CB/CellBender_Out_cell_barcodes.csv", sep = ",") # 8
colnames(p07wt04) <- "CellBarcode"
p07wt04$Index <- paste("YL_P07_WT_04", p07wt04$CellBarcode, sep = "_")

p07wt11 <- read.table("/work/Neuroinformatics_Core/akulk1/MANUAL/YL_DATA/CELLBENDER_MANUAL/P07_WT/YL_P07_WT_11_CB/CellBender_Out_cell_barcodes.csv", sep = ",") #       11406
colnames(p07wt11) <- "CellBarcode"
p07wt11$Index <- paste("YL_P07_WT_11", p07wt11$CellBarcode, sep = "_")

p07wt16 <- read.table("/work/Neuroinformatics_Core/akulk1/MANUAL/YL_DATA/CELLBENDER_MANUAL/P07_WT/YL_P07_WT_16_CB/CellBender_Out_cell_barcodes.csv", sep = ",") #        6405
colnames(p07wt16) <- "CellBarcode"
p07wt16$Index <- paste("YL_P07_WT_16", p07wt16$CellBarcode, sep = "_")


#### P07 HU
p07hu03 <- read.table("/work/Neuroinformatics_Core/akulk1/MANUAL/YL_DATA/CELLBENDER_MANUAL/P07_HU/YL_P07_HU_03_CB/CellBender_Out_cell_barcodes.csv", sep = ",") # 10538
colnames(p07hu03) <- "CellBarcode"
p07hu03$Index <- paste("YL_P07_HU_03", p07hu03$CellBarcode, sep = "_")

p07hu05 <- read.table("/work/Neuroinformatics_Core/akulk1/MANUAL/YL_DATA/CELLBENDER_MANUAL/P07_HU/YL_P07_HU_05_CB/CellBender_Out_cell_barcodes.csv", sep = ",") #  8668
colnames(p07hu05) <- "CellBarcode"
p07hu05$Index <- paste("YL_P07_HU_05", p07hu05$CellBarcode, sep = "_")

p07hu12 <- read.table("/work/Neuroinformatics_Core/akulk1/MANUAL/YL_DATA/CELLBENDER_MANUAL/P07_HU/YL_P07_HU_12_CB/CellBender_Out_cell_barcodes.csv", sep = ",") # 11589
colnames(p07hu12) <- "CellBarcode"
p07hu12$Index <- paste("YL_P07_HU_12", p07hu12$CellBarcode, sep = "_")

p07hu17 <- read.table("/work/Neuroinformatics_Core/akulk1/MANUAL/YL_DATA/CELLBENDER_MANUAL/P07_HU/YL_P07_HU_17_CB/CellBender_Out_cell_barcodes.csv", sep = ",") #    13
colnames(p07hu17) <- "CellBarcode"
p07hu17$Index <- paste("YL_P07_HU_17", p07hu17$CellBarcode, sep = "_")


#### P07 KO
p07ko01 <- read.table("/work/Neuroinformatics_Core/akulk1/MANUAL/YL_DATA/CELLBENDER_MANUAL/P07_KO/YL_P07_KO_01_CB/CellBender_Out_cell_barcodes.csv", sep = ",") #  9512
colnames(p07ko01) <- "CellBarcode"
p07ko01$Index <- paste("YL_P07_KO_01", p07ko01$CellBarcode, sep = "_")

p07ko06 <- read.table("/work/Neuroinformatics_Core/akulk1/MANUAL/YL_DATA/CELLBENDER_MANUAL/P07_KO/YL_P07_KO_06_CB/CellBender_Out_cell_barcodes.csv", sep = ",") #    16
colnames(p07ko06) <- "CellBarcode"
p07ko06$Index <- paste("YL_P07_KO_06", p07ko06$CellBarcode, sep = "_")

p07ko09 <- read.table("/work/Neuroinformatics_Core/akulk1/MANUAL/YL_DATA/CELLBENDER_MANUAL/P07_KO/YL_P07_KO_09_CB/CellBender_Out_cell_barcodes.csv", sep = ",") # 10599
colnames(p07ko09) <- "CellBarcode"
p07ko09$Index <- paste("YL_P07_KO_09", p07ko09$CellBarcode, sep = "_")

p07ko18 <- read.table("/work/Neuroinformatics_Core/akulk1/MANUAL/YL_DATA/CELLBENDER_MANUAL/P07_KO/YL_P07_KO_18_CB/CellBender_Out_cell_barcodes.csv", sep = ",") #  7934
colnames(p07ko18) <- "CellBarcode"
p07ko18$Index <- paste("YL_P07_KO_18", p07ko18$CellBarcode, sep = "_")


#### P56 WT
p56wt01 <- read.table("/work/Neuroinformatics_Core/akulk1/MANUAL/YL_DATA/CELLBENDER_MANUAL/P56_WT/YL_P56_WT_01_CB/CellBender_Out_cell_barcodes.csv", sep = ",") # 10029
colnames(p56wt01) <- "CellBarcode"
p56wt01$Index <- paste("YL_P56_WT_01", p56wt01$CellBarcode, sep = "_")

p56wt08 <- read.table("/work/Neuroinformatics_Core/akulk1/MANUAL/YL_DATA/CELLBENDER_MANUAL/P56_WT/YL_P56_WT_08_CB/CellBender_Out_cell_barcodes.csv", sep = ",") # 12924
colnames(p56wt08) <- "CellBarcode"
p56wt08$Index <- paste("YL_P56_WT_08", p56wt08$CellBarcode, sep = "_")

p56wt10 <- read.table("/work/Neuroinformatics_Core/akulk1/MANUAL/YL_DATA/CELLBENDER_MANUAL/P56_WT/YL_P56_WT_10_CB/CellBender_Out_cell_barcodes.csv", sep = ",") # 12592
colnames(p56wt10) <- "CellBarcode"
p56wt10$Index <- paste("YL_P56_WT_10", p56wt10$CellBarcode, sep = "_")

p56wt12 <- read.table("/work/Neuroinformatics_Core/akulk1/MANUAL/YL_DATA/CELLBENDER_MANUAL/P56_WT/YL_P56_WT_12_CB/CellBender_Out_cell_barcodes.csv", sep = ",") #  8882
colnames(p56wt12) <- "CellBarcode"
p56wt12$Index <- paste("YL_P56_WT_12", p56wt12$CellBarcode, sep = "_")


#### P56 HU
p56hu02 <- read.table("/work/Neuroinformatics_Core/akulk1/MANUAL/YL_DATA/CELLBENDER_MANUAL/P56_HU/YL_P56_HU_02_CB/CellBender_Out_cell_barcodes.csv", sep = ",") #  9904
colnames(p56hu02) <- "CellBarcode"
p56hu02$Index <- paste("YL_P56_HU_02", p56hu02$CellBarcode, sep = "_")

p56hu05 <- read.table("/work/Neuroinformatics_Core/akulk1/MANUAL/YL_DATA/CELLBENDER_MANUAL/P56_HU/YL_P56_HU_05_CB/CellBender_Out_cell_barcodes.csv", sep = ",") #  7608
colnames(p56hu05) <- "CellBarcode"
p56hu05$Index <- paste("YL_P56_HU_05", p56hu05$CellBarcode, sep = "_")

p56hu09 <- read.table("/work/Neuroinformatics_Core/akulk1/MANUAL/YL_DATA/CELLBENDER_MANUAL/P56_HU/YL_P56_HU_09_CB/CellBender_Out_cell_barcodes.csv", sep = ",") # 14343
colnames(p56hu09) <- "CellBarcode"
p56hu09$Index <- paste("YL_P56_HU_09", p56hu09$CellBarcode, sep = "_")

p56hu11 <- read.table("/work/Neuroinformatics_Core/akulk1/MANUAL/YL_DATA/CELLBENDER_MANUAL/P56_HU/YL_P56_HU_11_CB/CellBender_Out_cell_barcodes.csv", sep = ",") #  8256
colnames(p56hu11) <- "CellBarcode"
p56hu11$Index <- paste("YL_P56_HU_11", p56hu11$CellBarcode, sep = "_")


#### P56 KO
p56ko03 <- read.table("/work/Neuroinformatics_Core/akulk1/MANUAL/YL_DATA/CELLBENDER_MANUAL/P56_KO/YL_P56_KO_03_CB/CellBender_Out_cell_barcodes.csv", sep = ",") #  9304
colnames(p56ko03) <- "CellBarcode"
p56ko03$Index <- paste("YL_P56_KO_03", p56ko03$CellBarcode, sep = "_")

p56ko04 <- read.table("/work/Neuroinformatics_Core/akulk1/MANUAL/YL_DATA/CELLBENDER_MANUAL/P56_KO/YL_P56_KO_04_CB/CellBender_Out_cell_barcodes.csv", sep = ",") # 10521
colnames(p56ko04) <- "CellBarcode"
p56ko04$Index <- paste("YL_P56_KO_04", p56ko04$CellBarcode, sep = "_")

p56ko07 <- read.table("/work/Neuroinformatics_Core/akulk1/MANUAL/YL_DATA/CELLBENDER_MANUAL/P56_KO/YL_P56_KO_07_CB/CellBender_Out_cell_barcodes.csv", sep = ",") # 11679
colnames(p56ko07) <- "CellBarcode"
p56ko07$Index <- paste("YL_P56_KO_07", p56ko07$CellBarcode, sep = "_")

p56ko13 <- read.table("/work/Neuroinformatics_Core/akulk1/MANUAL/YL_DATA/CELLBENDER_MANUAL/P56_KO/YL_P56_KO_13_CB/CellBender_Out_cell_barcodes.csv", sep = ",") #  8596
colnames(p56ko13) <- "CellBarcode"
p56ko13$Index <- paste("YL_P56_KO_13", p56ko13$CellBarcode, sep = "_")


cbData <- rbind(p07wt02, p07wt04, p07wt11, p07wt16, 
                p07hu03, p07hu05, p07hu12, p07hu17,
                p07ko01, p07ko06, p07ko09, p07ko18, 
                p56wt01, p56wt08, p56wt10, p56wt12,
                p56hu02, p56hu05, p56hu09, p56hu11,
                p56ko03, p56ko04, p56ko07, p56ko13)

common <- intersect(row.names(dt.updated@meta.data), cbData$Index)  # 140638

cbDataCommon <- cbData[cbData$Index %in% common,]
cbDataCommon$CellBender <- rep("Retained", nrow(cbDataCommon))
row.names(cbDataCommon) <- cbDataCommon$Index

newMeta <- merge(dt.updated@meta.data, cbDataCommon, by = "row.names", all.x = TRUE)
newMeta$CellBender[is.na(newMeta$CellBender)] <- "Discarded"

cbMetaNew <- newMeta$CellBender
names(cbMetaNew) <- newMeta$Row.names

dt.updated$CellBender <- cbMetaNew

save(dt.updated, file = "CLOCKSEQ_P07_P56_WT_HU_KO_CELLBENDER.RData")


## FILTER AMBIENT RNA
dt.cb <- subset(dt.updated, subset = CellBender == "Retained")

save(dt.cb, file = "CLOCKSEQ_P07_P56_WT_HU_KO_CELLBENDER_RETAINED.RData")

```



-----
