##**************************
## BINF*6110 - Assignment 4
##
## Student Name: Cyrus Akbarally
##
## Student Number: 1099054
##
## 2026-04-06
##
##**************************

rm(list = ls())

## _ Install packages----
## install.packages("tidyverse")
## install.packages("viridis")
## install.packages("vegan")
## install.packages("dplyr")
## install.packages("ggplot2")
## install.packages("patchwork")
## install.packages("pheatmap")
## install.packages("BiocManager")
## install.packages("Seurat")
## install.packages("ggridges")
## devtools::install_github("immunogenomics/presto")
## BiocManager::install("SingleR")
## BiocManager::install("celldex")
## BiocManager::install("DESeq2")
## BiocManager::install("ashr")
## BiocManager::install("clusterProfiler")
## BiocManager::install("org.Mm.eg.db")
## BiocManager::install("enrichplot")

#--


## _ Packages used----
library(tidyverse)
library(readr)
library(viridis)
library(vegan)
library(dplyr)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(Seurat)
library(presto)
library(SingleR)
library(celldex)
library(DESeq2)
library(ashr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(ggridges)
## setwd("~/BINF6110 Genomic Methods for Bioinformatics/Assignment_4")

# _ Import data----
# Seurat object provided by Dr. Elias Taylor, University of Guelph (2026)
mouse_seurat <- readRDS("seurat_ass4.rds")

# _ Inspect metadata----
View(mouse_seurat@meta.data)


#--


### 1. UMAP OF MOUSE DATA----

# Normalization and Scaling
mouse_seurat_UMAP <- NormalizeData(mouse_seurat)
mouse_seurat_UMAP <- FindVariableFeatures(mouse_seurat_UMAP)
mouse_seurat_UMAP <- ScaleData(mouse_seurat_UMAP)

mouse_seurat_UMAP <- RunPCA(mouse_seurat_UMAP)

ElbowPlot(mouse_seurat_UMAP) + ggtitle("Elbow Plot of Mouse Influenza Mucosa Data")

mouse_seurat_UMAP <- FindNeighbors(mouse_seurat_UMAP, dims = 1:15, reduction = "pca")
mouse_seurat_UMAP <- FindClusters(mouse_seurat_UMAP, resolution = 0.5, cluster.name = "mouse_clusters")

mouse_seurat_UMAP <- RunUMAP(mouse_seurat_UMAP, dims = 1:15, reduction = "pca")
DimPlot(mouse_seurat_UMAP, reduction = "umap", group.by = "mouse_clusters") + ggtitle("UMAP of Mouse Influenza Mucosa Clusters")

# Cell counts for each cluster
table(mouse_seurat_UMAP$mouse_clusters)

# Obtain all unique values for categorical variables (for reference)
unique(mouse_seurat_UMAP$orig.ident)
unique(mouse_seurat_UMAP$biosample_id)
unique(mouse_seurat_UMAP$organ_custom)
unique(mouse_seurat_UMAP$disease__ontology_label)
unique(mouse_seurat_UMAP$time)
unique(mouse_seurat_UMAP$mouse_id)


#--


### 2. ANNOTATION OF CELLS BY GENE MARKERS----


## Using SingleR to automatically annotate all cells and clusters


## Annotate By Cell
mouse_annotate_cell <- GetAssayData(mouse_seurat_UMAP, assay = "RNA", layer = "data")
mouse_ref <- celldex::MouseRNAseqData()

mouse_singleR_cell <- SingleR(test = mouse_annotate_cell, ref = mouse_ref, labels = mouse_ref$label.main)

## Visualize annotated cells with annotated UMAP
mouse_seurat_UMAP$Label_SingleR_Cell <- mouse_singleR_cell$labels
DimPlot(mouse_seurat_UMAP, group.by = "Label_SingleR_Cell") + ggtitle("Mouse Cell Annotation - SingleR")



## Annotate By Cluster
mouse_annotate_cluster <- AggregateExpression(mouse_seurat_UMAP, group.by = "mouse_clusters", assays = "RNA", slot = "data")$RNA

mouse_singleR_cluster <- SingleR(test = mouse_annotate_cluster, ref = mouse_ref, labels = mouse_ref$label.main)

rownames(mouse_singleR_cluster) <- sub("^g", "", rownames(mouse_singleR_cluster))

cluster_map <- setNames(mouse_singleR_cluster$labels, rownames(mouse_singleR_cluster))

cell_clusters <- as.character(mouse_seurat_UMAP$mouse_clusters)

cluster_labels <- cluster_map[cell_clusters]
names(cluster_labels) <- colnames(mouse_seurat_UMAP)

## Visualize annotated clusters with annotated UMAP
mouse_seurat_UMAP$Label_SingleR_Cluster <- cluster_labels
DimPlot(mouse_seurat_UMAP, group.by = "Label_SingleR_Cluster") + ggtitle("Mouse Cluster Annotation - SingleR")


## Summary table of annotations to analyze labels for verification by manual annotation
mouse_singleR_comparison <- table(mouse_seurat_UMAP$mouse_clusters, mouse_seurat_UMAP$Label_SingleR_Cluster)

mouse_singleR_comparison <- as.data.frame(mouse_singleR_comparison, stringsAsFactors = FALSE)
names(mouse_singleR_comparison) <- c("Cluster", "Label", "Cell Frequency")

mouse_singleR_comparison <- mouse_singleR_comparison %>%
  group_by(Cluster) %>%
  slice_max(`Cell Frequency`) %>%
  ungroup() %>%
  arrange(as.numeric(Cluster))

## Summary table of number of clusters and total cell frequencies per label
mouse_singleR_label_summary <- as.data.frame(print(mouse_singleR_comparison %>%
                                                     group_by(Label) %>%
                                                     summarise(
                                                       n_clusters = n_distinct(Cluster),
                                                       total_cells = sum(`Cell Frequency`)
                                                       ) %>%
                                                     mutate(percent = (total_cells / sum(total_cells)) * 100)))


#--


### 3. FEATURE PLOTS OF TOP GENE MARKERS----

# Determine all gene markers for each cluster
mouse_markers <- FindAllMarkers(mouse_seurat_UMAP)
head(mouse_markers)

# Filter the top genes per cluster by adjusted p-value, average log2fold change, and percent expression

# Top 3 genes per cluster (can also be used for manual annotation of cell types)
top_genes_3 <- mouse_markers %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 0.5, p_val_adj < 0.05, pct.1 > 0.250, pct.2 < 0.100) %>%
  arrange(desc(avg_log2FC), p_val_adj) %>%
  slice_head(n = 3)

top_gene_list_3 <- top_genes_3$gene

# Topmost gene per cluster
top_genes <- mouse_markers %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 0.5, p_val_adj < 0.05, pct.1 > 0.250, pct.2 < 0.100) %>%
  arrange(desc(avg_log2FC), p_val_adj) %>%
  slice_head(n = 1)

top_gene_list <- top_genes$gene


#Topmost gene by cell type
top_genes$Cell_Type <- mouse_singleR_comparison$Label
top_genes_cell_type <- top_genes %>%
  group_by(Cell_Type) %>%
  slice_head(n = 1)

top_genes_feature <- top_genes_cell_type$gene

# Feature plots of top genes for each cell type
FeaturePlot(mouse_seurat_UMAP, features = top_genes_feature[c(1, 5:7, 8, 10)], reduction = "umap", raster = TRUE, ncol = 3) + plot_annotation(title = "Feature Plots of Genes of Interest - Immune Cells", theme = theme(plot.title = element_text(face = "bold", hjust = 0.5)))

FeaturePlot(mouse_seurat_UMAP, features = top_genes_feature[c(2:4, 9)], reduction = "umap", raster = TRUE, ncol = 2) + plot_annotation(title = "Feature Plots of Genes of Interest - Structural/Specialized Cells", theme = theme(plot.title = element_text(face = "bold", hjust = 0.5)))


#--


### 4. DIFFERENTIAL EXPRESSION OF A MOUSE CLUSTER----

## Subset cluster 4 (B cells) and pseudobulk for differential expression analysis
GSEA_DE_cluster <- subset(mouse_seurat_UMAP, idents = 4)
table(GSEA_DE_cluster$mouse_id)

mouse_pb <- AggregateExpression(GSEA_DE_cluster, assays = "RNA", group.by = c("mouse_id", "disease__ontology_label", "organ_custom", "time"), slot = "counts", return.seurat = TRUE)

View(mouse_pb@meta.data)


## Obtaining matrix and metadata for DESeq2
cluster_DE_counts <- GetAssayData(mouse_pb, layer = "counts")
mouse_meta <- mouse_pb@meta.data

## Remove 11 rows with NA from both matrix and metadata
na_rows <- which(!complete.cases(mouse_meta))
mouse_meta <- mouse_meta[-na_rows, ]
cluster_DE_counts <- cluster_DE_counts[, -na_rows]

## Sanity check to see if TRUE
all(colnames(cluster_DE_counts) == rownames(mouse_meta))

## Check for alignment
mouse_meta_clean <- mouse_meta[complete.cases(mouse_meta), , drop = FALSE]
cluster_DE_counts_clean <- cluster_DE_counts[, rownames(mouse_meta_clean), drop = FALSE]

## Sanity check to see if TRUE
all(colnames(cluster_DE_counts_clean) == rownames(mouse_meta_clean))

mouse_meta_clean$mouse_sample <- sub("-.*", "", mouse_meta_clean$mouse_id)


## Differential Expression Analysis with DESeq2
mouse_meta_clean$group <- interaction(mouse_meta_clean$disease__ontology_label, mouse_meta_clean$organ_custom, mouse_meta_clean$time)

dds <- DESeqDataSetFromMatrix(countData = cluster_DE_counts_clean, colData = mouse_meta_clean, design = ~ disease__ontology_label + organ_custom + time)
dds <- DESeq(dds)

resultsNames(dds)


## Sample pairwise comparisons to compare gene expression differences between each variable

## Differential expression based on time (time)
res_dds_time <- results(dds, name = "group_influenza.LNG.D14_vs_influenza.LNG.D02")
res_dds_time

## Differential expression based on tissue type (organ_custom)
res_dds_tissue1 <- results(dds, name = "group_influenza.OM.D02_vs_influenza.LNG.D02")
res_dds_tissue1

res_dds_tissue2 <- results(dds, name = "group_influenza.RM.D02_vs_influenza.LNG.D02")
res_dds_tissue2

## Differential expression before and post-infection (disease__ontology_label)
res_dds_disease <- results(dds, name = "group_normal.LNG.Naive_vs_influenza.LNG.D02")
res_dds_disease


## Visualize genes with heatmaps

## Apply ashr shrinkage to log2FoldChange for visualization
resLFCtime <- lfcShrink(dds, coef = "group_influenza.LNG.D14_vs_influenza.LNG.D02", type = "ashr")
resLFCtissue1 <- lfcShrink(dds, coef = "group_influenza.OM.D02_vs_influenza.LNG.D02", type = "ashr")
resLFCtissue2 <- lfcShrink(dds, coef = "group_influenza.RM.D02_vs_influenza.LNG.D02", type = "ashr")
resLFCdisease <- lfcShrink(dds, coef = "group_normal.LNG.Naive_vs_influenza.LNG.D02", type = "ashr")

## Remove NA rows
resLFCtime <- na.omit(resLFCtime)
resLFCtissue1 <- na.omit(resLFCtissue1)
resLFCtissue2 <- na.omit(resLFCtissue2)
resLFCdisease <- na.omit(resLFCdisease)

## Order genes based on descending log2FoldChange
top_genes_time <- head(order(resLFCtime$log2FoldChange), decreasing = TRUE)
gene_names_time <- rownames(resLFCtime)[top_genes_time]

top_genes_tissue1 <- head(order(resLFCtissue1$log2FoldChange), decreasing = TRUE)
gene_names_tissue1 <- rownames(resLFCtissue1)[top_genes_tissue1]

top_genes_tissue2 <- head(order(resLFCtissue2$log2FoldChange), decreasing = TRUE)
gene_names_tissue2 <- rownames(resLFCtissue2)[top_genes_tissue2]

top_genes_disease <- head(order(resLFCdisease$log2FoldChange), decreasing = TRUE)
gene_names_disease <- rownames(resLFCdisease)[top_genes_disease]

vsd <- vst(dds)

mat_time <- assay(vsd)[gene_names_time, ]
mat_tissue1 <- assay(vsd)[gene_names_tissue1, ]
mat_tissue2 <- assay(vsd)[gene_names_tissue2, ]
mat_disease <- assay(vsd)[gene_names_disease, ]

annotation_df <- mouse_meta_clean[, "group", drop = FALSE]
colnames(annotation_df) <- c("group")

## Time heatmap
pheatmap(mat_time,
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = annotation_df,
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "Differential Expression Heatmap for Time After Infection (D014 vs D02 (Reference))"
)


## Tissue (1) heatmap
pheatmap(mat_tissue1,
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = annotation_df,
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "Differential Expression Heatmap for Tissue Type (OM vs LNG (Reference))"
)

## Tissue (2) heatmap
pheatmap(mat_tissue2,
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = annotation_df,
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "Differential Expression Heatmap for Tissue Type (RM vs LNG (Reference))"
)

## Disease Ontology heatmap
pheatmap(mat_disease,
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = annotation_df,
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "Differential Expression Heatmap for Disease Ontology Before and Post-Infection (Normal (Naive) vs Influenza (D02)(Reference))"
)


#--


### 5. GO-GSEA FUNCTIONAL ANNOTATION OF A MOUSE CLUSTER----


## View the original results with no shrinkage as data frames----
res_dftime_no <- as.data.frame(res_dds_time)

res_dftissue1_no <- as.data.frame(res_dds_tissue1)
res_dftissue2_no <- as.data.frame(res_dds_tissue2)

res_dfdisease_no <- as.data.frame(res_dds_disease)


## Remove NA rows
res_dftime_no <- na.omit(res_dftime_no)

res_dftissue1_no <- na.omit(res_dftissue1_no)
res_dftissue2_no <- na.omit(res_dftissue2_no)

res_dfdisease_no <- na.omit(res_dfdisease_no)


res_dftime_no$SYMBOL <- rownames(res_dftime_no)

res_dftissue1_no$SYMBOL <- rownames(res_dftissue1_no)
res_dftissue2_no$SYMBOL <- rownames(res_dftissue2_no)

res_dfdisease_no$SYMBOL <- rownames(res_dfdisease_no)


## Convert gene names to Entrez IDs
gene_map_time <- bitr(res_dftime_no$SYMBOL,
                  fromType = "SYMBOL",
                  toType = c("ENTREZID"),
                  OrgDb = org.Mm.eg.db
)

gene_map_tissue1 <- bitr(res_dftissue1_no$SYMBOL,
                      fromType = "SYMBOL",
                      toType = c("ENTREZID"),
                      OrgDb = org.Mm.eg.db
)

gene_map_tissue2 <- bitr(res_dftissue2_no$SYMBOL,
                         fromType = "SYMBOL",
                         toType = c("ENTREZID"),
                         OrgDb = org.Mm.eg.db
)

gene_map_disease <- bitr(res_dfdisease_no$SYMBOL,
                  fromType = "SYMBOL",
                  toType = c("ENTREZID"),
                  OrgDb = org.Mm.eg.db
)

res_dftime_no <- inner_join(res_dftime_no, gene_map_time, by = "SYMBOL")
res_dftissue1_no <- inner_join(res_dftissue1_no, gene_map_time, by = "SYMBOL")
res_dftissue2_no <- inner_join(res_dftissue2_no, gene_map_tissue2, by = "SYMBOL")
res_dfdisease_no <- inner_join(res_dfdisease_no, gene_map_disease, by = "SYMBOL")


## Ranked lists by Wald statistic from DESeq2----
wald_time <- res_dftime_no$stat
genes_time <- res_dftime_no$ENTREZID
ranked_time <- setNames(wald_time, genes_time)
ranked_time <- na.omit(ranked_time)
ranked_time <- sort(ranked_time, decreasing = TRUE)

wald_tissue1 <- res_dftissue1_no$stat
genes_tissue1 <- res_dftissue1_no$ENTREZID
ranked_tissue1 <- setNames(wald_tissue1, genes_tissue1)
ranked_tissue1 <- na.omit(ranked_tissue1)
ranked_tissue1 <- sort(ranked_tissue1, decreasing = TRUE)

wald_tissue2 <- res_dftissue2_no$stat
genes_tissue2 <- res_dftissue2_no$ENTREZID
ranked_tissue2 <- setNames(wald_tissue2, genes_tissue2)
ranked_tissue2 <- na.omit(ranked_tissue2)
ranked_tissue2 <- sort(ranked_tissue2, decreasing = TRUE)


wald_disease <- res_dfdisease_no$stat
genes_disease <- res_dfdisease_no$ENTREZID
ranked_disease <- setNames(wald_disease, genes_disease)
ranked_disease <- na.omit(ranked_disease)
ranked_disease <- sort(ranked_disease, decreasing = TRUE)


## GO-GSEA comparison----

gsea_time <- gseGO(ranked_time,
                   OrgDb = org.Mm.eg.db,
                   ont = "BP",
                   keyType = "ENTREZID",
)

gsea_tissue1 <- gseGO(ranked_tissue1,
                      OrgDb = org.Mm.eg.db,
                      ont = "BP",
                      keyType = "ENTREZID",
)

gsea_tissue2 <- gseGO(ranked_tissue2,
                      OrgDb = org.Mm.eg.db,
                      ont = "BP",
                      keyType = "ENTREZID",
)

gsea_disease <- gseGO(ranked_disease,
                      OrgDb = org.Mm.eg.db,
                      ont = "BP",
                      keyType = "ENTREZID",
)

## Dot plots for GO-GSEA visualization----

dotplot(gsea_time, showCategory = 15, split = ".sign") +
  facet_grid(. ~ .sign) +
  ggtitle("GO Biological Process - Time")

dotplot(gsea_tissue1, showCategory = 15, split = ".sign") +
  facet_grid(. ~ .sign) +
  ggtitle("GO Biological Process - OM vs LNG Tissues")

dotplot(gsea_tissue2, showCategory = 15, split = ".sign") +
  facet_grid(. ~ .sign) +
  ggtitle("GO Biological Process - RM vs LNG Tissues")

dotplot(gsea_disease, showCategory = 15, split = ".sign") +
  facet_grid(. ~ .sign) +
  ggtitle("GO Biological Process - Disease Ontology (Before and Post-Infection)")


## Ridge plots for GO-GSEA visualization----

ridgeplot(gsea_time) +
  labs(title = "GO Biological Process - Time")

ridgeplot(gsea_tissue1) +
  labs(title = "GO Biological Process - OM vs LNG Tissues")

ridgeplot(gsea_tissue2) +
  labs(title = "GO Biological Process - RM vs LNG Tissues")

ridgeplot(gsea_disease) +
  labs(title = "GO Biological Process - Disease Ontology (Before and Post-Infection)")
