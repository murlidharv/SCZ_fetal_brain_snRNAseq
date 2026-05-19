# ============================================================
# 07_cluster_GE_integrated.R
# Clustering integrated ganglionic eminence object
# ============================================================

source("scripts/00_setup.R")
library(Seurat)

message(">>> Loading integrated GE object")

ge <- readRDS(
  file.path(DATA_PROC, "GE_integrated_annotated.rds")
)

# ------------------------------------------------------------
# Clustering on integrated assay
# ------------------------------------------------------------

DefaultAssay(ge) <- "integrated"

ge <- ScaleData(ge, verbose = FALSE)
ge <- RunPCA(ge, npcs = 30, verbose = FALSE)

# Inspect elbow if you want:
# ElbowPlot(ge)

ge <- FindNeighbors(ge, dims = 1:20)
ge <- FindClusters(ge, resolution = 0.5)

# ------------------------------------------------------------
# UMAP + t-SNE
# ------------------------------------------------------------

ge <- RunUMAP(ge, dims = 1:20)
ge <- RunTSNE(ge, dims = 1:20)

# ------------------------------------------------------------
# Save clustered object
# ------------------------------------------------------------

out <- file.path(DATA_PROC, "GE_integrated_clustered.rds")
saveRDS(ge, out)

message(">>> DONE")
message(out)
