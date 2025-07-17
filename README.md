# 🔬 Differential Gene Expression in T1D vs Healthy Whole Blood Samples (GSE123658)

## 📘 Background
This project analyzes RNA-seq data from the GEO dataset [**GSE123658**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123658), which contains **whole blood samples** from:

- **33 patients with Type 1 Diabetes (T1D)**
- **14 healthy controls**

The dataset includes **raw gene-level read counts**, allowing for downstream statistical analysis of gene expression differences between these two groups.

---

## ⚙️ Methods

We used the **DESeq2** package in R for the following steps:

1. **Data Preparation**
   - Parsed raw count data from GEO supplementary files
   - Matched sample metadata (condition: T1D or Healthy)

2. **DESeq2 Analysis**
   - Constructed a DESeq2 object using raw counts and sample condition
   - Performed normalization, dispersion estimation, and differential expression testing
   - Filtered genes with adjusted p-value (padj) < 0.05

3. **Visualization**
   - Plotted normalized expression of the most significant gene
   - (Optional) Volcano plot and heatmap available in scripts

---

## 📁 Output Files

| File                                 | Description                                             |
|--------------------------------------|---------------------------------------------------------|
| `DESeq2_full_results.csv`            | All genes tested, including log2 fold change and p-values |
| `DESeq2_significant_results.csv`     | Subset of genes with padj < 0.05                        |
| `DESeq2_significant_with_symbols.csv`| Significant genes with gene symbols annotated           |
| `volcano_plot.png` (optional)        | Volcano plot of all genes (if generated)               |
| `heatmap_top50.png` (optional)       | Heatmap of top 50 DE genes (if generated)              |

---

## 📊 Summary of Results

- **Genes analyzed:** 16,781
- **Significant (padj < 0.05):** 2,508 genes  
  - Upregulated in T1D: ~12%  
  - Downregulated in T1D: ~9.6%  

The top differentially expressed genes include several with known immune functions.
