############################################################
# 08_annotate_and_visualise_GE_clusters.R
# Purpose:
# - Assign biological labels to Seurat clusters
#   using mean Ganglionic Eminence marker module scores
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
# 1. Load clustered Ganglionic Eminence object
############################################################

ge <- readRDS(
  file.path(DATA_PROC, "GE_integrated_clustered.rds")
)

DefaultAssay(ge) <- "RNA"

############################################################
# 2. Identify GE marker score columns
############################################################

score_cols <- grep("^GE_marker", colnames(ge@meta.data), value = TRUE)
stopifnot(length(score_cols) > 0)
stopifnot("seurat_clusters" %in% colnames(ge@meta.data))

############################################################
# 3. Mean marker score per Seurat cluster
############################################################

cluster_scores <- ge@meta.data %>%
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
# 5. Explicit marker → biological label mapping (GE)
############################################################

marker_to_celltype <- c(
  GE_marker1  = "GE-CycPro",
  GE_marker2  = "GE-InN-1",
  GE_marker3  = "GE-InN-2",
  GE_marker4  = "GE-InN-3",
  GE_marker5  = "GE-InN-4",
  GE_marker6  = "GE-InN-5",
  GE_marker7  = "GE-InN-6",
  GE_marker8  = "GE-InN-7",
  GE_marker9  = "GE-RG-1",
  GE_marker10 = "GE-RG-2",
  GE_marker11 = "GE-oRG"
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

ge$cluster_celltype <- unname(
  cluster_to_celltype[as.character(ge$seurat_clusters)]
)

############################################################
# 7. Attach CELL-LEVEL biological labels (STANDARDISED)
############################################################

# This makes GE identical to FC / Cer / Thal for downstream scripts
ge$cell_type <- ge$cluster_celltype

############################################################
# 8. Sanity checks
############################################################

message(">>> Cluster-level cell-type distribution:")
print(table(ge$cluster_celltype))

message(">>> Cell-level cell-type distribution:")
print(table(ge$cell_type))

message(">>> NA cluster labels:")
print(sum(is.na(ge$cluster_celltype)))

message(">>> NA cell-level labels:")
print(sum(is.na(ge$cell_type)))

############################################################
# 9. Save object + tables
############################################################

saveRDS(
  ge,
  file.path(DATA_PROC, "GE_with_cluster_celltype_annotations.rds")
)

saveRDS(
  ge,
  file.path(DATA_PROC, "GE_with_celltype_annotations.rds")
)

write.csv(
  cluster_scores,
  file.path(DATA_PROC, "GE_cluster_mean_marker_scores.csv"),
  row.names = FALSE
)

write.csv(
  best_labels,
  file.path(DATA_PROC, "GE_cluster_best_celltype.csv"),
  row.names = FALSE
)

############################################################
# 10. Visualisation
############################################################

p_umap_clusters <- DimPlot(
  ge,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) + ggtitle("Ganglionic Eminence – Seurat clusters (UMAP)")

p_umap_celltype <- DimPlot(
  ge,
  reduction = "umap",
  group.by = "cell_type",
  label = TRUE,
  repel = TRUE
) + ggtitle("Ganglionic Eminence – cell types (UMAP)")

p_tsne_clusters <- DimPlot(
  ge,
  reduction = "tsne",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) + ggtitle("Ganglionic Eminence – Seurat clusters (t-SNE)")

p_tsne_celltype <- DimPlot(
  ge,
  reduction = "tsne",
  group.by = "cell_type",
  label = TRUE,
  repel = TRUE
) + ggtitle("Ganglionic Eminence – cell types (t-SNE)")

ggsave(
  file.path(FIG_DIR, "GE_UMAP_clusters_vs_celltype.png"),
  p_umap_clusters + p_umap_celltype,
  width = 12, height = 6, dpi = 300
)

ggsave(
  file.path(FIG_DIR, "GE_tSNE_clusters_vs_celltype.png"),
  p_tsne_clusters + p_tsne_celltype,
  width = 12, height = 6, dpi = 300
)

message(">>> Ganglionic Eminence cluster annotation + visualisation complete")
############################################################
# END SCRIPT
############################################################
