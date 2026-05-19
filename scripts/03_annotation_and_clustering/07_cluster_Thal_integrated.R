# ============================================================
# 07_cluster_Thal_integrated.R
# Clustering integrated thalamus object
# ============================================================

source("scripts/00_setup.R")
library(Seurat)

message(">>> Loading integrated Thal object")

thal <- readRDS(
  file.path(DATA_PROC, "Thal_integrated_annotated.rds")
)

# ------------------------------------------------------------
# Clustering on integrated assay
# ------------------------------------------------------------

DefaultAssay(thal) <- "integrated"

thal <- ScaleData(thal, verbose = FALSE)
thal <- RunPCA(thal, npcs = 30, verbose = FALSE)

# Inspect elbow if you want:
# ElbowPlot(thal)

thal <- FindNeighbors(thal, dims = 1:20)
thal <- FindClusters(thal, resolution = 0.5)

# ------------------------------------------------------------
# UMAP + t-SNE
# ------------------------------------------------------------

thal <- RunUMAP(thal, dims = 1:20)
thal <- RunTSNE(thal, dims = 1:20)

# ------------------------------------------------------------
# Save clustered object
# ------------------------------------------------------------

out <- file.path(DATA_PROC, "Thal_integrated_clustered.rds")
saveRDS(thal, out)

message(">>> DONE")
message(out)

