############################################################
# 08_annotate_and_visualise_Hipp_clusters.R
# Purpose:
# - Assign biological labels to Seurat clusters
#   using mean Hippocampus marker module scores
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
# 1. Load clustered Hippocampus object
############################################################

hipp <- readRDS(
  file.path(DATA_PROC, "Hipp_integrated_annotated.rds")
)

DefaultAssay(hipp) <- "RNA"

############################################################
# 2. Identify Hipp marker score columns
############################################################

score_cols <- grep("^Hipp_marker", colnames(hipp@meta.data), value = TRUE)
stopifnot(length(score_cols) > 0)
stopifnot("seurat_clusters" %in% colnames(hipp@meta.data))

############################################################
# 3. Mean marker score per Seurat cluster
############################################################

cluster_scores <- hipp@meta.data %>%
  dplyr::select(seurat_clusters, all_of(score_cols)) %>%
  dplyr::group_by(seurat_clusters) %>%
  dplyr::summarise(
    dplyr::across(all_of(score_cols), mean, na.rm = TRUE),
    .groups = "drop"
  )

############################################################
# 4. Pick best-scoring marker per cluster
############################################################

best_labels <- cluster_scores %>%
  tidyr::pivot_longer(
    cols = all_of(score_cols),
    names_to = "marker",
    values_to = "mean_score"
  ) %>%
  dplyr::group_by(seurat_clusters) %>%
  dplyr::slice_max(mean_score, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()

############################################################
# 5. Explicit marker → biological label mapping (Hipp)
############################################################

marker_to_celltype <- c(
  Hipp_marker1  = "Hipp-CR-1",
  Hipp_marker2  = "Hipp-CR-2",
  Hipp_marker3  = "Hipp-CycPro",
  Hipp_marker4  = "Hipp-Endo",
  Hipp_marker5  = "Hipp-ExN-1",
  Hipp_marker6  = "Hipp-ExN-2",
  Hipp_marker7  = "Hipp-ExN-3",
  Hipp_marker8  = "Hipp-ExN-4",
  Hipp_marker9  = "Hipp-ExN-5",
  Hipp_marker10 = "Hipp-ExN-6",
  Hipp_marker11 = "Hipp-InN-1",
  Hipp_marker12 = "Hipp-InN-2",
  Hipp_marker13 = "Hipp-InN-3",
  Hipp_marker14 = "Hipp-InN-4",
  Hipp_marker15 = "Hipp-MG",
  Hipp_marker16 = "Hipp-OPC",
  Hipp_marker17 = "Hipp-RG-1",
  Hipp_marker18 = "Hipp-RG-2",
  Hipp_marker19 = "Hipp-oRG"
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

hipp$cluster_celltype <- unname(
  cluster_to_celltype[as.character(hipp$seurat_clusters)]
)

############################################################
# 7. Attach CELL-LEVEL biological labels (STANDARDISED)
############################################################

hipp$cell_type <- hipp$cluster_celltype
hipp$Hipp_celltype <- hipp$cluster_celltype

############################################################
# 8. Sanity checks
############################################################

message(">>> Cluster-level cell-type distribution:")
print(table(hipp$cluster_celltype))

message(">>> Cell-level cell-type distribution:")
print(table(hipp$cell_type))

message(">>> NA cluster labels:")
print(sum(is.na(hipp$cluster_celltype)))

############################################################
# 9. Save object + tables
############################################################

saveRDS(
  hipp,
  file.path(DATA_PROC, "Hipp_with_cluster_celltype_annotations.rds")
)

saveRDS(
  hipp,
  file.path(DATA_PROC, "Hipp_integrated_annotated.rds")
)

write.csv(
  cluster_scores,
  file.path(DATA_PROC, "Hipp_cluster_mean_marker_scores.csv"),
  row.names = FALSE
)

write.csv(
  best_labels,
  file.path(DATA_PROC, "Hipp_cluster_best_celltype.csv"),
  row.names = FALSE
)

############################################################
# 10. Visualisation
############################################################

p_umap_clusters <- DimPlot(
  hipp,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) + ggtitle("Hippocampus – Seurat clusters (UMAP)")

p_umap_celltype <- DimPlot(
  hipp,
  reduction = "umap",
  group.by = "cell_type",
  label = TRUE,
  repel = TRUE
) + ggtitle("Hippocampus – cell types (UMAP)")

p_tsne_clusters <- DimPlot(
  hipp,
  reduction = "tsne",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) + ggtitle("Hippocampus – Seurat clusters (t-SNE)")

p_tsne_celltype <- DimPlot(
  hipp,
  reduction = "tsne",
  group.by = "cell_type",
  label = TRUE,
  repel = TRUE
) + ggtitle("Hippocampus – cell types (t-SNE)")

ggsave(
  file.path(FIG_DIR, "Hipp_UMAP_clusters_vs_celltype.png"),
  p_umap_clusters + p_umap_celltype,
  width = 12, height = 6, dpi = 300
)

ggsave(
  file.path(FIG_DIR, "Hipp_tSNE_clusters_vs_celltype.png"),
  p_tsne_clusters + p_tsne_celltype,
  width = 12, height = 6, dpi = 300
)

message(">>> Hippocampus cluster annotation + visualisation complete")
############################################################
# END SCRIPT
############################################################
