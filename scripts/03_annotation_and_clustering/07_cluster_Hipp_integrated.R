# ============================================================
# 07_cluster_Hipp_integrated.R
# Clustering integrated hippocampus object
# ============================================================

source("scripts/00_setup.R")
library(Seurat)

message(">>> Loading integrated Hipp object")

hipp <- readRDS(
  file.path(DATA_PROC, "Hipp_allSamples_integrated.rds")
)

# ------------------------------------------------------------
# Clustering on integrated assay
# ------------------------------------------------------------

DefaultAssay(hipp) <- "integrated"

hipp <- ScaleData(hipp, verbose = FALSE)
hipp <- RunPCA(hipp, npcs = 30, verbose = FALSE)

# Optional:
# ElbowPlot(hipp)

hipp <- FindNeighbors(hipp, dims = 1:20)
hipp <- FindClusters(hipp, resolution = 0.5)

# ------------------------------------------------------------
# UMAP + t-SNE
# ------------------------------------------------------------

hipp <- RunUMAP(hipp, dims = 1:20)
hipp <- RunTSNE(hipp, dims = 1:20)

# ------------------------------------------------------------
# Save clustered object
# ------------------------------------------------------------

out <- file.path(DATA_PROC, "Hipp_integrated_clustered.rds")
saveRDS(hipp, out)

message(">>> DONE")
message(out)
