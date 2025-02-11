load(".RData")
tab <- seuObj

metaTemp <- as.data.frame(tab@meta.data)
metaTemp$CellBarcode <- rownames(metaTemp) 
metaTemp$clusterid <- metaTemp$RNA_snn_res.0.8
annotations <- read.table("annotation.txt", header = TRUE, sep = "\t")
annotations2 <- merge(metaTemp, annotations, by = "clusterid")
annotations3 <- annotations2$celltype
names(annotations3) <- annotations2$CellBarcode
tab[["CellType"]] <- annotations3

DefaultAssay(tab) <- "RNA"
cpm <- NormalizeData(tab,scale.factor = 1e6)
cpm.sub <- subset(cpm, subset = CellType %in% "#") #    ENDO  Oligos   EX  ASTRO   IN  MICRO Oligos
log_data <- GetAssayData(cpm.sub)
log2_data <- log(expm1(log_data) + 1, 2)
t = t(log2_data)
exp_m <- t[,c("CLOCK")]
result[n,i] <- mean(exp_m)
