# ======================================
# 03_normalization.R
# ======================================

source("scripts/00_setup.R")

# -------------------------------
# Sample ID
# -------------------------------
sample_id <- "993_Hipp_B2"

rds_in  <- file.path(DATA_PROC, paste0(sample_id, "_filtered.rds"))
rds_out <- file.path(DATA_PROC, paste0(sample_id, "_normalized.rds"))

stopifnot(file.exists(rds_in))

# -------------------------------
# Load filtered object
# -------------------------------
seurat.obj <- readRDS(rds_in)

#DefaultAssay(seurat.obj) <- "RNA" 

# -------------------------------
# Normalization
# -------------------------------
seurat.obj <- NormalizeData(
  seurat.obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

# -------------------------------
# Identify highly variable genes
# -------------------------------
seurat.obj <- FindVariableFeatures(
  seurat.obj,
  selection.method = "vst",
  nfeatures = 2000
)

# -------------------------------
# (Optional) Inspect HVGs
# -------------------------------
VariableFeaturePlot(seurat.obj)

top10 <- head(VariableFeatures(seurat.obj), 10)
LabelPoints(
  plot = VariableFeaturePlot(seurat.obj),
  points = top10,
  repel = TRUE
)

# -------------------------------
# Save normalized object
# -------------------------------
saveRDS(seurat.obj, rds_out)
