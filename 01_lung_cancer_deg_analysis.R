# Differential Gene Expression Analysis in Lung Cancer
# Project: TCGA-LUAD RNA-seq analysis using DESeq2, GO, and KEGG
# Author: Poreddy Gopika Reddy

# -----------------------------
# 0. Project setup
# -----------------------------

# Set this path to your local project folder
setwd("C:/Users/DELL/Desktop/Lung_Cancer_DEG_Project")

dir.create("data", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("plots", showWarnings = FALSE)
dir.create("report", showWarnings = FALSE)
dir.create("scripts", showWarnings = FALSE)

# -----------------------------
# 1. Load required packages
# -----------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

required_packages <- c(
  "TCGAbiolinks",
  "SummarizedExperiment",
  "DESeq2",
  "ggplot2",
  "pheatmap",
  "dplyr",
  "clusterProfiler",
  "org.Hs.eg.db",
  "enrichplot"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# -----------------------------
# 2. Download TCGA-LUAD RNA-seq data
# -----------------------------

query_all <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)

results <- getResults(query_all)
table(results$sample_type)

# Select a computationally manageable subset:
# 60 Primary Tumor samples + all Solid Tissue Normal samples
normal_samples <- results$cases[results$sample_type == "Solid Tissue Normal"]
tumor_samples <- results$cases[results$sample_type == "Primary Tumor"]

set.seed(123)
tumor_subset <- sample(tumor_samples, 60)
selected_samples <- c(tumor_subset, normal_samples)

query_small <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode = selected_samples
)

options(timeout = 100000)
GDCdownload(query_small, method = "api", files.per.chunk = 5)

luad_data <- GDCprepare(query_small)
saveRDS(luad_data, file = "data/TCGA_LUAD_subset_STAR_counts.rds")

# -----------------------------
# 3. Metadata preparation
# -----------------------------

metadata_simple <- data.frame(
  sample_id = colnames(luad_data),
  patient_id = substr(colnames(luad_data), 1, 12),
  sample_type = as.character(colData(luad_data)$sample_type),
  stringsAsFactors = FALSE
)

metadata_simple$condition <- ifelse(
  metadata_simple$sample_type == "Primary Tumor",
  "Tumor",
  "Normal"
)

table(metadata_simple$sample_type)
table(metadata_simple$condition)

write.csv(
  metadata_simple,
  file = "results/TCGA_LUAD_subset_metadata.csv",
  row.names = FALSE
)

# -----------------------------
# 4. Count matrix preparation
# -----------------------------

assayNames(luad_data)

count_matrix <- assay(luad_data, "unstranded")
dim(count_matrix)

metadata_deseq <- metadata_simple
rownames(metadata_deseq) <- metadata_deseq$sample_id
metadata_deseq <- metadata_deseq[, "condition", drop = FALSE]
metadata_deseq$condition <- factor(metadata_deseq$condition, levels = c("Normal", "Tumor"))

all(colnames(count_matrix) == rownames(metadata_deseq))

write.csv(count_matrix, file = "results/TCGA_LUAD_count_matrix.csv")
write.csv(metadata_deseq, file = "results/TCGA_LUAD_metadata_DESeq2.csv")

saveRDS(count_matrix, file = "data/TCGA_LUAD_count_matrix.rds")
saveRDS(metadata_deseq, file = "data/TCGA_LUAD_metadata_DESeq2.rds")

# -----------------------------
# 5. Differential expression analysis using DESeq2
# -----------------------------

dds <- DESeqDataSetFromMatrix(
  countData = round(count_matrix),
  colData = metadata_deseq,
  design = ~ condition
)

# Remove low-count genes
dds <- dds[rowSums(counts(dds)) >= 10, ]

dds <- DESeq(dds)

res <- results(
  dds,
  contrast = c("condition", "Tumor", "Normal")
)

res_ordered <- res[order(res$padj), ]
res_df <- as.data.frame(res_ordered)
res_df$gene_id <- rownames(res_df)

summary(res)

deg_results <- res_df %>%
  filter(!is.na(padj)) %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1)

upregulated_genes <- deg_results %>%
  filter(log2FoldChange > 1)

downregulated_genes <- deg_results %>%
  filter(log2FoldChange < -1)

nrow(deg_results)
nrow(upregulated_genes)
nrow(downregulated_genes)

saveRDS(dds, file = "data/TCGA_LUAD_DESeq2_dds.rds")
saveRDS(res, file = "data/TCGA_LUAD_DESeq2_results.rds")

write.csv(res_df, file = "results/TCGA_LUAD_all_DESeq2_results.csv", row.names = FALSE)
write.csv(deg_results, file = "results/TCGA_LUAD_significant_DEGs.csv", row.names = FALSE)
write.csv(upregulated_genes, file = "results/TCGA_LUAD_upregulated_genes.csv", row.names = FALSE)
write.csv(downregulated_genes, file = "results/TCGA_LUAD_downregulated_genes.csv", row.names = FALSE)

# -----------------------------
# 6. PCA plot
# -----------------------------

vsd <- vst(dds, blind = FALSE)

pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA Plot of TCGA-LUAD Tumor and Normal Samples") +
  theme_minimal()

ggsave(
  filename = "plots/TCGA_LUAD_PCA_plot.png",
  plot = pca_plot,
  width = 7,
  height = 5,
  dpi = 300
)

# -----------------------------
# 7. Sample distance heatmap
# -----------------------------

sample_dist <- dist(t(assay(vsd)))
sample_dist_matrix <- as.matrix(sample_dist)

annotation_col <- data.frame(condition = colData(vsd)$condition)
rownames(annotation_col) <- colnames(vsd)

png(
  filename = "plots/TCGA_LUAD_sample_distance_heatmap.png",
  width = 1200,
  height = 1000,
  res = 150
)

pheatmap(
  sample_dist_matrix,
  annotation_col = annotation_col,
  annotation_row = annotation_col,
  main = "Sample-to-Sample Distance Heatmap"
)

dev.off()

# -----------------------------
# 8. Volcano plot
# -----------------------------

volcano_data <- res_df %>%
  filter(!is.na(padj)) %>%
  mutate(
    regulation = case_when(
      padj < 0.05 & log2FoldChange > 1 ~ "Upregulated",
      padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
      TRUE ~ "Not Significant"
    )
  )

table(volcano_data$regulation)

volcano_plot <- ggplot(
  volcano_data,
  aes(x = log2FoldChange, y = -log10(padj), color = regulation)
) +
  geom_point(alpha = 0.7, size = 1.2) +
  scale_color_manual(
    values = c(
      "Upregulated" = "red",
      "Downregulated" = "blue",
      "Not Significant" = "grey"
    )
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(
    title = "Volcano Plot of Differentially Expressed Genes in TCGA-LUAD",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Gene Regulation"
  ) +
  theme_minimal()

ggsave(
  filename = "plots/TCGA_LUAD_volcano_plot.png",
  plot = volcano_plot,
  width = 8,
  height = 6,
  dpi = 300
)

# -----------------------------
# 9. Heatmap of top 50 DEGs
# -----------------------------

top50_genes <- deg_results %>%
  arrange(padj) %>%
  head(50)

top50_matrix <- assay(vsd)[top50_genes$gene_id, ]
top50_matrix_scaled <- t(scale(t(top50_matrix)))

annotation_col <- data.frame(condition = metadata_deseq$condition)
rownames(annotation_col) <- rownames(metadata_deseq)

png(
  filename = "plots/TCGA_LUAD_top50_DEGs_heatmap.png",
  width = 1200,
  height = 1000,
  res = 150
)

pheatmap(
  top50_matrix_scaled,
  annotation_col = annotation_col,
  show_colnames = FALSE,
  show_rownames = TRUE,
  main = "Heatmap of Top 50 Differentially Expressed Genes"
)

dev.off()

# -----------------------------
# 10. Gene annotation
# -----------------------------

gene_annotation <- as.data.frame(rowData(luad_data))
write.csv(
  gene_annotation,
  file = "results/TCGA_LUAD_gene_annotation.csv",
  row.names = TRUE
)

gene_annotation_simple <- data.frame(
  gene_id = rownames(gene_annotation),
  gene_name = gene_annotation$gene_name,
  gene_type = gene_annotation$gene_type,
  stringsAsFactors = FALSE
)

res_annotated <- res_df %>%
  left_join(gene_annotation_simple, by = "gene_id")

deg_annotated <- deg_results %>%
  left_join(gene_annotation_simple, by = "gene_id")

upregulated_annotated <- upregulated_genes %>%
  left_join(gene_annotation_simple, by = "gene_id")

downregulated_annotated <- downregulated_genes %>%
  left_join(gene_annotation_simple, by = "gene_id")

write.csv(res_annotated, file = "results/TCGA_LUAD_all_DESeq2_results_annotated.csv", row.names = FALSE)
write.csv(deg_annotated, file = "results/TCGA_LUAD_significant_DEGs_annotated.csv", row.names = FALSE)
write.csv(upregulated_annotated, file = "results/TCGA_LUAD_upregulated_genes_annotated.csv", row.names = FALSE)
write.csv(downregulated_annotated, file = "results/TCGA_LUAD_downregulated_genes_annotated.csv", row.names = FALSE)

up_gene_symbols <- upregulated_annotated %>%
  filter(!is.na(gene_name), gene_name != "") %>%
  pull(gene_name) %>%
  unique()

down_gene_symbols <- downregulated_annotated %>%
  filter(!is.na(gene_name), gene_name != "") %>%
  pull(gene_name) %>%
  unique()

all_deg_symbols <- deg_annotated %>%
  filter(!is.na(gene_name), gene_name != "") %>%
  pull(gene_name) %>%
  unique()

length(up_gene_symbols)
length(down_gene_symbols)
length(all_deg_symbols)

write.table(up_gene_symbols, file = "results/upregulated_gene_symbols.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(down_gene_symbols, file = "results/downregulated_gene_symbols.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(all_deg_symbols, file = "results/all_DEG_gene_symbols.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

# -----------------------------
# 11. GO and KEGG enrichment analysis
# -----------------------------

deg_entrez <- bitr(
  all_deg_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

deg_entrez_unique <- unique(deg_entrez$ENTREZID)

go_enrich <- enrichGO(
  gene = deg_entrez_unique,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

go_results <- as.data.frame(go_enrich)

write.csv(
  go_results,
  file = "results/TCGA_LUAD_GO_enrichment_results.csv",
  row.names = FALSE
)

go_summary <- go_results[, c(
  "ID",
  "Description",
  "GeneRatio",
  "BgRatio",
  "pvalue",
  "p.adjust",
  "qvalue",
  "Count"
)]

write.csv(
  go_summary,
  file = "results/TCGA_LUAD_GO_enrichment_summary.csv",
  row.names = FALSE
)

go_dotplot <- dotplot(go_enrich, showCategory = 15) +
  ggtitle("GO Biological Process Enrichment of DEGs")

ggsave(
  filename = "plots/TCGA_LUAD_GO_dotplot.png",
  plot = go_dotplot,
  width = 9,
  height = 7,
  dpi = 300
)

kegg_enrich <- enrichKEGG(
  gene = deg_entrez_unique,
  organism = "hsa",
  pvalueCutoff = 0.05
)

kegg_enrich <- setReadable(
  kegg_enrich,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)

kegg_results <- as.data.frame(kegg_enrich)

write.csv(
  kegg_results,
  file = "results/TCGA_LUAD_KEGG_enrichment_results.csv",
  row.names = FALSE
)

kegg_summary <- kegg_results[, c(
  "ID",
  "Description",
  "GeneRatio",
  "BgRatio",
  "pvalue",
  "p.adjust",
  "qvalue",
  "Count"
)]

write.csv(
  kegg_summary,
  file = "results/TCGA_LUAD_KEGG_enrichment_summary.csv",
  row.names = FALSE
)

kegg_dotplot <- dotplot(kegg_enrich, showCategory = 15) +
  ggtitle("KEGG Pathway Enrichment of DEGs")

ggsave(
  filename = "plots/TCGA_LUAD_KEGG_dotplot.png",
  plot = kegg_dotplot,
  width = 9,
  height = 7,
  dpi = 300
)

# -----------------------------
# 12. Final file check
# -----------------------------

list.files("plots")
list.files("results")

cat("\nAnalysis completed successfully.\n")
cat("Significant DEGs:", nrow(deg_results), "\n")
cat("Upregulated genes:", nrow(upregulated_genes), "\n")
cat("Downregulated genes:", nrow(downregulated_genes), "\n")
