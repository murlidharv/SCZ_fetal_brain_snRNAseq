############################################################
# 08_annotate_and_visualise_FC_clusters.R
# Purpose:
# - Assign biological labels to Seurat clusters
#   using mean FC marker module scores
# - Visualise clusters on UMAP and t-SNE
############################################################

source("scripts/00_setup.R")

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

############################################################
# 1. Load clustered FC object
############################################################

fc <- readRDS(
  file.path(DATA_PROC, "FC_integrated_clustered.rds")
)

DefaultAssay(fc) <- "RNA"

############################################################
# 2. Identify FC marker score columns
############################################################

score_cols <- grep("^FC_marker", colnames(fc@meta.data), value = TRUE)
stopifnot(length(score_cols) > 0)
stopifnot("seurat_clusters" %in% colnames(fc@meta.data))

############################################################
# 3. Mean marker score per Seurat cluster
############################################################

cluster_scores <- fc@meta.data %>%
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
# 5. Explicit marker → biological label mapping
############################################################

marker_to_celltype <- c(
  FC_marker1  = "FC-CycPro",
  FC_marker2  = "FC-RG-1",
  FC_marker3  = "FC-ExN-1",
  FC_marker4  = "FC-ExN-2",
  FC_marker5  = "FC-ExN-3",
  FC_marker6  = "FC-ExN-4",
  FC_marker7  = "FC-ExN-5",
  FC_marker8  = "FC-InN-1",
  FC_marker9  = "FC-InN-2",
  FC_marker10 = "FC-InN-3",
  FC_marker11 = "FC-InN-4",
  FC_marker12 = "FC-IP",
  FC_marker13 = "FC-RG-2",
  FC_marker14 = "FC-OPC",
  FC_marker15 = "FC-Endo",
  FC_marker16 = "FC-MG"
)

best_labels$cluster_celltype <- marker_to_celltype[best_labels$marker]

stopifnot(!any(is.na(best_labels$cluster_celltype)))

############################################################
# 6. Attach CLUSTER-LEVEL biological labels (SAFE)
############################################################

cluster_to_celltype <- setNames(
  best_labels$cluster_celltype,
  as.character(best_labels$seurat_clusters)
)

fc$cluster_celltype <- unname(
  cluster_to_celltype[as.character(fc$seurat_clusters)]
)

############################################################
# 7. Sanity checks
############################################################

message(">>> Cluster-level cell-type distribution:")
print(table(fc$cluster_celltype))

message(">>> NA cluster labels:")
print(sum(is.na(fc$cluster_celltype)))

head(fc@meta.data[, c("seurat_clusters", "FC_celltype", "cluster_celltype")])

############################################################
# 8. Save object + tables
############################################################

saveRDS(
  fc,
  file.path(DATA_PROC, "FC_with_cluster_celltype_annotations.rds")
)

write.csv(
  cluster_scores,
  file.path(DATA_PROC, "FC_cluster_mean_marker_scores.csv"),
  row.names = FALSE
)

write.csv(
  best_labels,
  file.path(DATA_PROC, "FC_cluster_best_celltype.csv"),
  row.names = FALSE
)

############################################################
# 9. Visualisation
############################################################

p_umap_clusters <- DimPlot(
  fc,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) + ggtitle("Frontal cortex – Seurat clusters (UMAP)")

p_umap_celltype <- DimPlot(
  fc,
  reduction = "umap",
  group.by = "cluster_celltype",
  label = TRUE,
  repel = TRUE
) + ggtitle("Frontal cortex – cluster-level cell types (UMAP)")

p_tsne_clusters <- DimPlot(
  fc,
  reduction = "tsne",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) + ggtitle("Frontal cortex – Seurat clusters (t-SNE)")

p_tsne_celltype <- DimPlot(
  fc,
  reduction = "tsne",
  group.by = "cluster_celltype",
  label = TRUE,
  repel = TRUE
) + ggtitle("Frontal cortex – cluster-level cell types (t-SNE)")

ggsave(
  file.path(FIG_DIR, "FC_UMAP_clusters_vs_celltype.png"),
  p_umap_clusters + p_umap_celltype,
  width = 12, height = 6, dpi = 300
)

ggsave(
  file.path(FIG_DIR, "FC_tSNE_clusters_vs_celltype.png"),
  p_tsne_clusters + p_tsne_celltype,
  width = 12, height = 6, dpi = 300
)

message(">>> FC cluster annotation + visualisation complete")
