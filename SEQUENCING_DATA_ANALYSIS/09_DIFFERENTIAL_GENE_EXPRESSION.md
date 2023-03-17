# DIFFERENTIAL GENE EXPRESSION ANALYSIS

```{R}
library(plyr)
library(dplyr)
library(tidyverse)
library(tidyr)
library(Seurat)
library(patchwork)
library(Matrix.utils)
library(ggpubr)
library(reshape2)
library(data.table)
library(rio)
library(scater)

#load clustered dataset and subset the cell types and genotypes for comparison
load("CLUSTERED.RData")
mycelltype <- c("cell type")
sel <- subset(clock.combined, subset = CellType %in% mycelltype)
myGenoType <- c("HU", "WT")
sel <- subset(sel, subset = Genotype %in% myGenoType)
sub = sel
Idents(sub) = sub$Genotype

## Prepare data for MAST ##
mat <- sub@assays$RNA@data
metause = sub@meta.data

# Keep genes with 10% expression in at least one species
cells1 = WhichCells(sub, idents = 'HU')
cells2 = WhichCells(sub, idents = 'WT')
pass1 = rowSums(mat[,cells1])/length(cells1) > 0.1
pass2 = rowSums(mat[,cells2])/length(cells2) > 0.1
mat = mat[pass1 | pass2,]

# Create MAST object
sca <- MAST::FromMatrix(exprsArray = as.matrix(x = mat),
                        cData = metause,
                        fData = data.frame(rownames(mat)))
# Calculate fold change
avg_logfc <- log(rowMeans(expm1(mat[,cells1])) + 1) - log(rowMeans(expm1(mat[,cells2])) + 1)

# MAST package with covariates and Hurdle model
mastfix = MAST::zlm(~Genotype + RNAPrep + Sex, sca)
summaryCond <- summary(object = mastfix, doLRT = 'GenotypeWT')
summaryDt <- summaryCond$datatable
p_val <- summaryDt[summaryDt$component == "H", 4]
genes.return <- summaryDt[summaryDt$component == "H", 1]
to.return <- data.frame(p_val, row.names = genes.return$primerid)
colnames(to.return)[1] = 'p_value'
fix_res = to.return
fix_res$adj_p_value = p.adjust(fix_res$p_value, method = 'BH')
fix_res$avg_logfc = avg_logfc

# MAST package with covariates and glmer and mixed model
mastfix_glmer = MAST::zlm(~Genotype + (1|Sample) + RNAPrep + Sex, sca, ebayes = F, method = 'glmer')
summaryCond <- summary(object = mastfix_glmer, doLRT = 'GenotypeWT')
summaryDt <- summaryCond$datatable
p_val <- summaryDt[summaryDt$component == "H", 4]
genes.return <- summaryDt[summaryDt$component == "H", 1]
to.return <- data.frame(p_val, row.names = genes.return$primerid)
colnames(to.return)[1] = 'p_value'
fix_res_glmer = to.return
fix_res_glmer$adj_p_value = p.adjust(fix_res$p_value, method = 'BH')
fix_res_glmer$avg_logfc = avg_logfc

######use Seurat to calculate folder change and expressed percentage of each gene
DefaultAssay(sel) <- "RNA"
latentVars <- c("nCount_RNA", "pMito_RNA", "Sex", "RNAPrep") 
sel@active.ident <- factor(sel@meta.data$Genotype)
names(sel@active.ident) <- row.names(sel@meta.data)
avg.exp <- log1p(AverageExpression(sel, verbose = FALSE)$RNA)
avg.exp$Gene <- row.names(avg.exp)
print(dim(avg.exp))
clock.deg.MAST <- FindMarkers(sel, ident.1 = "HU", ident.2 = "WT", test.use = "MAST", verbose = TRUE, logfc.threshold = 0, only.pos = FALSE, latent.vars = latentVars)
clock.deg.MAST <- subset(clock.deg.MAST, avg_log2FC > 0.1 | avg_log2FC < -0.1)

#merge the results of different methods by gene name and save the results for further analyses
hurdle <- fix_res
seurat <- clock.deg.MAST
glmer <- mastfix_glmer
colnames(hurdle) <- c("pval_hurdle", "adj_pval_hurdle", "logFC_hurdle")
colnames(seurat) <- c("pval_seurat", "logFC_seurat", "pct1", "pct2", "adj_pval_seurat")
colnames(glmer) <- c("pval_glmer", "adj_pval_glmer", "logFC_glmer")
hurdle$Gene <- rownames(hurdle)
seurat$Gene <- rownames(seurat)
glmer$Gene <- rownames(glmer)
combined <- merge(glmer, seurat, by = "Gene")
combined <- merge(combined, hurdle, by = "Gene")
write.table(combined, "DEGs_celltype.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

```

-----