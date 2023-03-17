# GENE ONTOLOGY ENRICHMENT

```{R}
library(reshape2)
library(ggplot2)
library(WGCNA)

#combine the output of TOPPgene results by using half lists of genes separately 
tab1 <- read.table(file = "result_from_using_of_half_expressed_genes", sep = "\t", header = T, fill = T, quote="")
tab1 <- tab1[, c(1,2,3,9,10)]
tab2 <- read.table(file = "result_from_using_of_the_other_half_expressed_genes", sep = "\t", header = T, fill = T, quote="")
tab2 <- tab2[, c(1,2,3,9,10)]
tab1 <- merge(aggregate(Hit.Count.in.Genome~ID, data = tab1, FUN = max), tab1)
tab2 <- merge(aggregate(Hit.Count.in.Genome~ID, data = tab2, FUN = max), tab2)
tab <- rbind(tab1, tab2)
tab_order <- tab[order(tab$ID),]

#Combine the hit number of genes in each GO term from the two batches of enrichment
matrix <- data.frame(matrix(ncol = 5, nrow = 0))
m = 1
n = 1
l = length(tab$ID)+1
while (n < l) {
  if (n == length(tab$ID)) {
    matrix[m,1] <- tab_order[n,1]
    matrix[m,2] <- tab_order[n,2]
    matrix[m,3] <- tab_order[n,3]
    matrix[m,4] <- tab_order[n,4]
    matrix[m,5] <- tab_order[n,5]
  } else if (tab_order[n,1] == tab_order[n+1,1]) {
    count_q = tab_order[n,5] + tab_order[n+1,5]
    count_t = tab_order[n,2] + tab_order[n+1,2]
    matrix[m,1] <- tab_order[n,1]
    matrix[m,2] <- count_t
    matrix[m,3] <- tab_order[n,3]
    matrix[m,4] <- tab_order[n,4]
    matrix[m,5] <- count_q
    n = n + 1
  } else {
    matrix[m,1] <- tab_order[n,1]
    matrix[m,2] <- tab_order[n,2]
    matrix[m,3] <- tab_order[n,3]
    matrix[m,4] <- tab_order[n,4]
    matrix[m,5] <- tab_order[n,5]
  }
  
  m = m + 1
  n = n + 1
}

#conduct hypergeometric test and calculate adjusted p value
for (i in 1:length(matrix$X1)) {
  q = matrix[i,5] - 1
  m = matrix[i,2]
  n = number_of_expressed_genes - m
  k = 49
  stats = 1-phyper(q, m, n, k, lower.tail = TRUE, log.p = FALSE)
  matrix[i,6] <- stats
}
colnames(matrix) <- c("ID", "#_GO", "Category","Name", "Hit", "pval")
matrix$holm <- p.adjust(matrix$pval, "holm")
matrix$hochberg <- p.adjust(matrix$pval, "hochberg")
matrix$hommel <- p.adjust(matrix$pval, "hommel")
matrix$bonferroni <- p.adjust(matrix$pval, "bonferroni")
matrix$BH <- p.adjust(matrix$pval, "BH")
matrix$BY <- p.adjust(matrix$pval, "BY")
matrix$fdr <- p.adjust(matrix$pval, "fdr")
matrix$dataset <- "P56_KW_DOWN"

write.table(matrix, "GO_enrichment.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
```

-----
