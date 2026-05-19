############################################################
# 08_annotate_and_visualise_Thal_clusters.R
# Purpose:
# - Assign biological labels to Seurat clusters
#   using mean Thalamus marker module scores
# - Attach BOTH cluster-level and cell-level labels
# - Visualise clusters on UMAP and t-SNE
############################################################

source("scripts/00_setup.R")

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

############################################################
# 1. Load clustered Thalamus object
############################################################

thal <- readRDS(
  file.path(DATA_PROC, "Thal_integrated_clustered.rds")
)

DefaultAssay(thal) <- "RNA"

############################################################
# 2. Identify Thalamus marker score columns
############################################################

score_cols <- grep("^Thal_marker", colnames(thal@meta.data), value = TRUE)
stopifnot(length(score_cols) > 0)
stopifnot("seurat_clusters" %in% colnames(thal@meta.data))

############################################################
# 3. Mean marker score per Seurat cluster
############################################################

cluster_scores <- thal@meta.data %>%
  select(seurat_clusters, all_of(score_cols)) %>%
  group_by(seurat_clusters) %>%
  summarise(
    across(all_of(score_cols), mean, na.rm = TRUE),
    .groups = "drop"
  )

############################################################
# 4. Pick best-scoring marker per cluster
############################################################

best_labels <- cluster_scores %>%
  pivot_longer(
    cols = all_of(score_cols),
    names_to = "marker",
    values_to = "mean_score"
  ) %>%
  group_by(seurat_clusters) %>%
  slice_max(mean_score, n = 1, with_ties = FALSE) %>%
  ungroup()

############################################################
# 5. Explicit marker → biological label mapping (Thalamus)
############################################################

marker_to_celltype <- c(
  Thal_marker1  = "Thal-CycPro",
  Thal_marker2  = "Thal-Endo",
  Thal_marker3  = "Thal-ExN-1",
  Thal_marker4  = "Thal-ExN-2",
  Thal_marker5  = "Thal-ExN-3",
  Thal_marker6  = "Thal-IP",
  Thal_marker7  = "Thal-InN-1",
  Thal_marker8  = "Thal-InN-2",
  Thal_marker9  = "Thal-InN-3",
  Thal_marker10 = "Thal-InN-4",
  Thal_marker11 = "Thal-InN-5",
  Thal_marker12 = "Thal-InN-6",
  Thal_marker13 = "Thal-InN-7",
  Thal_marker14 = "Thal-InN-8",
  Thal_marker15 = "Thal-MG",
  Thal_marker16 = "Thal-OPC",
  Thal_marker17 = "Thal-RG-1",
  Thal_marker18 = "Thal-RG-2",
  Thal_marker19 = "Thal-oRG",
  Thal_marker20 = "Thal-RG-4",
  Thal_marker21 = "Thal-RG-5",
  Thal_marker22 = "Thal-RG-6",
  Thal_marker23 = "Thal-RG-7"
)

best_labels$cluster_celltype <- marker_to_celltype[best_labels$marker]
stopifnot(!any(is.na(best_labels$cluster_celltype)))

############################################################
# 6. Attach CLUSTER-LEVEL biological labels
############################################################

cluster_to_celltype <- setNames(
  best_labels$cluster_celltype,
  as.character(best_labels$seurat_clusters)
)

thal$cluster_celltype <- unname(
  cluster_to_celltype[as.character(thal$seurat_clusters)]
)

############################################################
# 7. Attach CELL-LEVEL biological labels (STANDARDISED)
############################################################

# This makes Thal identical to FC / Cer for downstream scripts
thal$cell_type <- thal$cluster_celltype

############################################################
# 8. Sanity checks
############################################################

message(">>> Cluster-level cell-type distribution:")
print(table(thal$cluster_celltype))

message(">>> Cell-level cell-type distribution:")
print(table(thal$cell_type))

message(">>> NA cluster labels:")
print(sum(is.na(thal$cluster_celltype)))

message(">>> NA cell-level labels:")
print(sum(is.na(thal$cell_type)))

############################################################
# 9. Save object + tables
############################################################

saveRDS(
  thal,
  file.path(DATA_PROC, "Thal_with_cluster_celltype_annotations.rds")
)

saveRDS(
  thal,
  file.path(DATA_PROC, "Thal_with_celltype_annotations.rds")
)

write.csv(
  cluster_scores,
  file.path(DATA_PROC, "Thal_cluster_mean_marker_scores.csv"),
  row.names = FALSE
)

write.csv(
  best_labels,
  file.path(DATA_PROC, "Thal_cluster_best_celltype.csv"),
  row.names = FALSE
)

############################################################
# 10. Visualisation
############################################################

p_umap_clusters <- DimPlot(
  thal,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) + ggtitle("Thalamus – Seurat clusters (UMAP)")

p_umap_celltype <- DimPlot(
  thal,
  reduction = "umap",
  group.by = "cell_type",
  label = TRUE,
  repel = TRUE
) + ggtitle("Thalamus – cell types (UMAP)")

p_tsne_clusters <- DimPlot(
  thal,
  reduction = "tsne",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) + ggtitle("Thalamus – Seurat clusters (t-SNE)")

p_tsne_celltype <- DimPlot(
  thal,
  reduction = "tsne",
  group.by = "cell_type",
  label = TRUE,
  repel = TRUE
) + ggtitle("Thalamus – cell types (t-SNE)")

ggsave(
  file.path(FIG_DIR, "Thal_UMAP_clusters_vs_celltype.png"),
  p_umap_clusters + p_umap_celltype,
  width = 12, height = 6, dpi = 300
)

ggsave(
  file.path(FIG_DIR, "Thal_tSNE_clusters_vs_celltype.png"),
  p_tsne_clusters + p_tsne_celltype,
  width = 12, height = 6, dpi = 300
)

message(">>> Thalamus cluster annotation + visualisation complete")
############################################################
# END SCRIPT
############################################################
