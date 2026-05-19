# ============================================================
# 07_cluster_FC_integrated.R
# Clustering integrated frontal cortex object
# ============================================================

source("scripts/00_setup.R")
library(Seurat)

message(">>> Loading integrated FC object")

fc <- readRDS(
  file.path(DATA_PROC, "FC_integrated_annotated.rds")
)

# ------------------------------------------------------------
# Clustering on integrated assay
# ------------------------------------------------------------

DefaultAssay(fc) <- "integrated"

fc <- ScaleData(fc, verbose = FALSE)
fc <- RunPCA(fc, npcs = 30, verbose = FALSE)

# Inspect elbow if you want:
# ElbowPlot(fc)

fc <- FindNeighbors(fc, dims = 1:20)
fc <- FindClusters(fc, resolution = 0.5)

# ------------------------------------------------------------
# UMAP + t-SNE
# ------------------------------------------------------------

fc <- RunUMAP(fc, dims = 1:20)
fc <- RunTSNE(fc, dims = 1:20)

# ------------------------------------------------------------
# Save clustered object
# ------------------------------------------------------------

out <- file.path(DATA_PROC, "FC_integrated_clustered.rds")
saveRDS(fc, out)

message(">>> DONE")
message(out)
