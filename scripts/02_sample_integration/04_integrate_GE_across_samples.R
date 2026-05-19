# ============================================================
# GE Integration Script (Seurat v5)
# ============================================================

source("scripts/00_setup.R")
library(Seurat)

message(">>> Loading normalized GE objects")

ge_510 <- readRDS(file.path(DATA_PROC, "510_GE_B1_normalized.rds"))
ge_611 <- readRDS(file.path(DATA_PROC, "611_GE_B1_normalized.rds"))
ge_993 <- readRDS(file.path(DATA_PROC, "993_GE_B2_normalized.rds"))

# ---------------------------------------------------------
# Seurat v5: collapse RNA layers
# ---------------------------------------------------------

message(">>> Joining RNA layers")

ge_510 <- JoinLayers(ge_510, assay = "RNA")
ge_611 <- JoinLayers(ge_611, assay = "RNA")
ge_993 <- JoinLayers(ge_993, assay = "RNA")

# ---------------------------------------------------------
# Ensure correct assay
# ---------------------------------------------------------

DefaultAssay(ge_510) <- "RNA"
DefaultAssay(ge_611) <- "RNA"
DefaultAssay(ge_993) <- "RNA"

# ---------------------------------------------------------
# Add sample labels
# ---------------------------------------------------------

ge_510$sample <- "510"
ge_611$sample <- "611"
ge_993$sample <- "993"

ge.list <- list(ge_510, ge_611, ge_993)
names(ge.list) <- c("510", "611", "993")

# ---------------------------------------------------------
# Integration
# ---------------------------------------------------------

message(">>> Selecting features")

features <- SelectIntegrationFeatures(
  object.list = ge.list,
  nfeatures = 2000
)

message(">>> Finding anchors")

anchors <- FindIntegrationAnchors(
  object.list = ge.list,
  anchor.features = features,
  dims = 1:20
)

message(">>> Integrating data")

ge_integrated <- IntegrateData(
  anchorset = anchors,
  dims = 1:20
)

# ---------------------------------------------------------
# Save
# ---------------------------------------------------------

out <- file.path(DATA_PROC, "GE_allSamples_integrated.rds")
saveRDS(ge_integrated, out)

message(">>> DONE")
message(out)
