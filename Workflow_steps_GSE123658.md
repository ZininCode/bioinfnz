# DESeq2 Analysis Workflow for GSE123658 (T1D vs Healthy)

This workflow documents the complete analysis of the GSE123658 dataset using DESeq2. It includes raw count processing, metadata extraction, normalization, differential expression analysis, and optional result visualization and export.  
The dataset compares **Type 1 Diabetes (T1D)** vs **Healthy** human blood samples.

---

## 1. Load Required Libraries

```r
library(DESeq2)
library(GEOquery)
library(biomaRt)
library(pheatmap)
```

---

## 2. Read the Raw Count Matrix

```r
count_matrix <- read.delim("GSE123658_read_counts.gene_level.txt.gz", row.names = 1)
```

*Loads raw gene-level counts with Ensembl gene IDs as row names.*

---

## 3. Explore Count Matrix Structure

```r
dim(count_matrix)           # Dimensions (genes × samples)
head(count_matrix[, 1:5])   # Preview values
```

*Inspect dimensions and data format (integer counts per gene per sample).*

---

## 4. Load GEO Metadata

```r
gse <- getGEO("GSE123658", GSEMatrix = TRUE)
pheno_data <- pData(gse[[1]])
```

*Downloads sample metadata.*

---

## 5. Explore Metadata Fields

```r
View(pheno_data)
head(pheno_data$title)
head(pheno_data$supplementary_file_1)
```

*Determine which fields contain condition labels (e.g. “T1D”) and sample identifiers.*

---

## 6. Build Sample Info Table (file_map)

```r
file_map <- data.frame(
  gsm = pheno_data$geo_accession,
  supp_file = basename(pheno_data$supplementary_file_1),
  stringsAsFactors = FALSE
)

file_map$short_id <- sub(".*_(.*)\.geneLevel.*", "\1", file_map$supp_file)
file_map$condition <- ifelse(grepl("T1D", pheno_data$title, ignore.case = TRUE), "T1D", "Healthy")
```

*Extract internal short IDs (used as column names in count_matrix) and condition labels.*

---

## 7. Filter to Matching Samples Only

```r
file_map_filtered <- file_map[file_map$short_id %in% colnames(count_matrix), ]
```

*Retain only metadata rows that match the column names of the count matrix.*

---

## 8. Create DESeq2 col_data (Metadata Table)

```r
col_data <- data.frame(
  condition = factor(file_map_filtered$condition),
  row.names = file_map_filtered$short_id,
  stringsAsFactors = FALSE
)
```

---

## 9. Subset Count Matrix to Matching Columns

```r
count_matrix <- count_matrix[, rownames(col_data)]
```

*Ensure the count matrix only contains samples with matching metadata.*

---

## 10. Check Compatibility

```r
all(colnames(count_matrix) == rownames(col_data))  # Should return TRUE
```

---

## 11. Create DESeq2 Dataset Object

```r
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = col_data,
  design = ~ condition
)
```

*This object now holds both raw counts and condition metadata.*

---

## 12. Run DESeq2 Analysis

```r
dds <- DESeq(dds)
```

*Performs size factor estimation, dispersion modeling, and negative binomial testing.*

---

## 13. Extract and View Results

```r
res <- results(dds)
res_ordered <- res[order(res$padj), ]
summary(res_ordered)
```

---

## 14. Filter for Significant Genes

```r
sig_genes <- res_ordered[which(res_ordered$padj < 0.05), ]
dim(sig_genes)  # Number of significant genes
```

---

## 15. Plot Normalized Counts of Top Gene

```r
plotCounts(dds, gene = rownames(sig_genes)[1], intgroup = "condition")
```

---

## 16. Optional: Annotate Gene IDs with Gene Names

```r
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_info <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = rownames(sig_genes),
  mart = ensembl
)

sig_annotated <- merge(
  as.data.frame(sig_genes),
  gene_info,
  by.x = "row.names",
  by.y = "ensembl_gene_id",
  all.x = TRUE
)
```

---

## 17. Optional: Volcano Plot

```r
with(res, plot(log2FoldChange, -log10(padj), pch = 20, main = "Volcano Plot"))
```

---

## 18. Optional: Export Significant Genes

```r
write.csv(as.data.frame(sig_genes), file = "significant_genes_T1D_vs_Healthy.csv")
```

---

## Notes

- Count data are **raw gene-level counts** from RNA-seq.
- Sample conditions ("T1D", "Healthy") are extracted from the sample titles.
- Only 47 of 82 samples had clear condition + matching file ID, and were retained.
- Further analysis can include heatmaps, clustering, or pathway enrichment.
