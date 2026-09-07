# Differential Gene Expression Analysis in Lung Cancer

## Project Overview

This project focuses on differential gene expression analysis in lung cancer using publicly available RNA-seq data from The Cancer Genome Atlas Lung Adenocarcinoma project (TCGA-LUAD). The study compares primary lung tumor samples with solid tissue normal samples to identify significantly upregulated and downregulated genes.

The analysis was performed using DESeq2 in R. Functional enrichment analysis was carried out using Gene Ontology Biological Process and KEGG pathway analysis to identify important biological processes and pathways involved in lung cancer progression.

## Dataset

- Database: TCGA
- Project: TCGA-LUAD
- Cancer type: Lung adenocarcinoma
- Data type: RNA-seq gene expression quantification
- Workflow type: STAR - Counts
- Tumor samples: 60
- Normal samples: 59
- Total samples: 119

## Tools and Packages Used

- RStudio
- TCGAbiolinks
- SummarizedExperiment
- DESeq2
- ggplot2
- pheatmap
- dplyr
- clusterProfiler
- org.Hs.eg.db
- enrichplot

## Methodology

1. RNA-seq data collection from TCGA-LUAD
2. Preparation of count matrix and sample metadata
3. Differential gene expression analysis using DESeq2
4. Identification of significant DEGs using adjusted p-value and log2 fold change thresholds
5. Visualization using PCA plot, sample distance heatmap, volcano plot, and heatmap
6. Gene annotation using gene symbols
7. GO Biological Process enrichment analysis
8. KEGG pathway enrichment analysis
9. Biological interpretation of significant genes and pathways

## DEG Analysis Criteria

Genes were considered significant based on:

- Adjusted p-value < 0.05
- Absolute log2 fold change > 1

## Results

A total of 12,324 significant differentially expressed genes were identified between lung tumor and normal samples.

| Category | Number of Genes |
|---|---:|
| Significant DEGs | 12,324 |
| Upregulated genes | 8,616 |
| Downregulated genes | 3,708 |

After gene symbol annotation:

| Category | Number of Gene Symbols |
|---|---:|
| Upregulated gene symbols | 8,564 |
| Downregulated gene symbols | 3,691 |
| All DEG gene symbols | 12,253 |

## Visualizations

The following plots were generated:

- PCA plot
- Sample distance heatmap
- Volcano plot
- Heatmap of top 50 differentially expressed genes
- GO enrichment dotplot
- KEGG enrichment dotplot

## Functional Enrichment Analysis

GO Biological Process enrichment showed involvement of:

- Extracellular matrix organization
- Extracellular structure organization
- Cell-cell adhesion
- Humoral immune response
- Immunoglobulin-mediated immune response
- B cell-mediated immunity

KEGG pathway enrichment showed involvement of:

- Neuroactive ligand-receptor interaction
- Cadherin signaling
- Cytoskeleton-related pathways
- Hormone signaling
- Neutrophil extracellular trap formation

## Conclusion

This project identified significant transcriptomic alterations between lung cancer and normal lung tissue samples. The results suggest that lung cancer progression is associated with extracellular matrix remodeling, altered cell adhesion, immune response changes, cytoskeletal reorganization, and signaling pathway regulation. The identified genes and pathways may be useful for further biomarker discovery and therapeutic target research.

## Repository Structure

```text
scripts/     R scripts used for analysis
results/     DEG and enrichment result files
plots/       PCA, volcano, heatmap, GO and KEGG plots
report/      Final project report
