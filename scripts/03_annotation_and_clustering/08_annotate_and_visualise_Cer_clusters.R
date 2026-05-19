############################################################
# 08_annotate_and_visualise_Cer_clusters.R
# Purpose:
# - Assign biological labels to Seurat clusters using
#   published marker module scores
# - Visualise clusters on UMAP and t-SNE
############################################################

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

setwd("~/colemdjl-exphaem/MSc_DL_Murali")

############################################################
# 1. Load clustered + scored object
############################################################

cer <- readRDS("data/processed/Cer_with_named_marker_scores.rds")
DefaultAssay(cer) <- "RNA"

############################################################
# 2. Identify score columns (robust)
############################################################

score_cols <- grep("^Score_", colnames(cer@meta.data), value = TRUE)
stopifnot(length(score_cols) > 0)

############################################################
# 3. Mean marker score per Seurat cluster
############################################################

cluster_scores <- cer@meta.data %>%
  select(seurat_clusters, all_of(score_cols)) %>%
  group_by(seurat_clusters) %>%
  summarise(
    across(all_of(score_cols), mean, na.rm = TRUE),
    .groups = "drop"
  )

############################################################
# 4. Pick best-scoring cell type per cluster
############################################################

best_labels <- cluster_scores %>%
  pivot_longer(
    cols = all_of(score_cols),
    names_to = "cell_type",
    values_to = "mean_score"
  ) %>%
  group_by(seurat_clusters) %>%
  slice_max(mean_score, n = 1, with_ties = FALSE) %>%
  ungroup()

# Clean names: Score_GC → GC
best_labels$cell_type <- gsub("^Score_", "", best_labels$cell_type)

############################################################
# 5. Attach biological labels to Seurat object
############################################################

cluster_to_celltype <- setNames(
  best_labels$cell_type,
  as.character(best_labels$seurat_clusters)
)

cer$cell_type <- unname(
  cluster_to_celltype[as.character(cer$seurat_clusters)]
)

# Check distribution
table(cer$cell_type)

# Check for NAs (important!)
sum(is.na(cer$cell_type))

# Inspect a few cells
head(cer@meta.data[, c("seurat_clusters", "cell_type")])

saveRDS(
  cer,
  "data/processed/Cer_with_celltype_annotations.rds"
)

############################################################
# 6. Save annotation tables
############################################################

write.csv(
  cluster_scores,
  "data/processed/Cer_cluster_mean_marker_scores.csv",
  row.names = FALSE
)

write.csv(
  best_labels,
  "data/processed/Cer_cluster_best_celltype.csv",
  row.names = FALSE
)

############################################################
# 7. Visualisation
############################################################

p_umap_clusters <- DimPlot(
  cer,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) + ggtitle("Cerebellum – Seurat clusters (UMAP)")

p_umap_celltype <- DimPlot(
  cer,
  reduction = "umap",
  group.by = "cell_type",
  label = TRUE,
  repel = TRUE
) + ggtitle("Cerebellum – annotated cell types (UMAP)")

p_tsne_clusters <- DimPlot(
  cer,
  reduction = "tsne",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) + ggtitle("Cerebellum – Seurat clusters (t-SNE)")

p_tsne_celltype <- DimPlot(
  cer,
  reduction = "tsne",
  group.by = "cell_type",
  label = TRUE,
  repel = TRUE
) + ggtitle("Cerebellum – annotated cell types (t-SNE)")

############################################################
# 8. Save figures
############################################################

ggsave("figures/Cer_UMAP_clusters_vs_celltype.png",
       p_umap_clusters + p_umap_celltype,
       width = 12, height = 6, dpi = 300)

ggsave("figures/Cer_tSNE_clusters_vs_celltype.png",
       p_tsne_clusters + p_tsne_celltype,
       width = 12, height = 6, dpi = 300)

############################################################
# 9. Save final annotated object
############################################################

