# ============================================================
# 08_cluster_Cer_integrated_with_tSNE.R
# Clustering of merged cerebellum dataset
# UMAP + t-SNE
# Seurat v5 compatible
# ============================================================

source("scripts/00_setup.R")

library(Seurat)
library(ggplot2)

options(bitmapType = "cairo")

# ------------------------------------------------------------
# Load integrated cerebellum object
# ------------------------------------------------------------

cer <- readRDS(
  file.path(DATA_PROC, "Cer_allSamples_integrated.rds")
)

# ------------------------------------------------------------
# Set assay for clustering
# ------------------------------------------------------------

DefaultAssay(cer) <- "integrated"

# ------------------------------------------------------------
# Scaling
# ------------------------------------------------------------

cer <- ScaleData(cer, verbose = FALSE)

# ------------------------------------------------------------
# PCA
# ------------------------------------------------------------

cer <- RunPCA(
  cer,
  npcs = 50,
  verbose = FALSE
)

# ------------------------------------------------------------
# Choose dimensions
# ------------------------------------------------------------

dims_use <- 1:30

# ------------------------------------------------------------
# Nearest neighbors
# ------------------------------------------------------------

cer <- FindNeighbors(
  cer,
  dims = dims_use
)

# ------------------------------------------------------------
# Clustering
# ------------------------------------------------------------

cer <- FindClusters(
  cer,
  resolution = 0.4
)

# ------------------------------------------------------------
# UMAP
# ------------------------------------------------------------

cer <- RunUMAP(
  cer,
  dims = dims_use,
  reduction = "pca"
)

# ------------------------------------------------------------
# t-SNE
# ------------------------------------------------------------
# t-SNE is run on PCA space (standard practice)

cer <- RunTSNE(
  cer,
  dims = dims_use,
  reduction = "pca",
  perplexity = 30,
  check_duplicates = FALSE
)

# ------------------------------------------------------------
# Save clustered object
# ------------------------------------------------------------

saveRDS(
  cer,
  file = file.path(DATA_PROC, "Cer_clustered.rds")
)

# ------------------------------------------------------------
# Plot: UMAP clusters
# ------------------------------------------------------------

p_umap <- DimPlot(
  cer,
  reduction = "umap",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Cerebellum clusters (UMAP)")

ggsave(
  filename = file.path(FIG_DIR, "Cer_UMAP_clusters.png"),
  plot = p_umap,
  width = 6,
  height = 5,
  dpi = 300
)

# ------------------------------------------------------------
# Plot: t-SNE clusters
# ------------------------------------------------------------

p_tsne <- DimPlot(
  cer,
  reduction = "tsne",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Cerebellum clusters (t-SNE)")

ggsave(
  filename = file.path(FIG_DIR, "Cer_tSNE_clusters.png"),
  plot = p_tsne,
  width = 6,
  height = 5,
  dpi = 300
)

# ------------------------------------------------------------
# Optional: sample mixing QC
# ------------------------------------------------------------

if ("sample" %in% colnames(cer@meta.data)) {
  
  p_umap_s <- DimPlot(
    cer,
    reduction = "umap",
    group.by = "sample"
  ) +
    ggtitle("UMAP by sample")
  
  ggsave(
    filename = file.path(FIG_DIR, "Cer_UMAP_by_sample.png"),
    plot = p_umap_s,
    width = 6,
    height = 5,
    dpi = 300
  )
  
  p_tsne_s <- DimPlot(
    cer,
    reduction = "tsne",
    group.by = "sample"
  ) +
    ggtitle("t-SNE by sample")
  
  ggsave(
    filename = file.path(FIG_DIR, "Cer_tSNE_by_sample.png"),
    plot = p_tsne_s,
    width = 6,
    height = 5,
    dpi = 300
  )
}

message(">>> Clustering with UMAP and t-SNE complete")
