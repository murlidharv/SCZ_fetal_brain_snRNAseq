# ============================================================
# Hipp Integration Script (Seurat v5)
# ============================================================

source("scripts/00_setup.R")
library(Seurat)

message(">>> Loading normalized Hipp objects")

hipp_510 <- readRDS(file.path(DATA_PROC, "510_Hipp_B2_normalized.rds"))
hipp_611 <- readRDS(file.path(DATA_PROC, "611_Hipp_B2_normalized.rds"))
hipp_993 <- readRDS(file.path(DATA_PROC, "993_Hipp_B2_normalized.rds"))

# ---------------------------------------------------------
# Seurat v5: collapse RNA layers
# ---------------------------------------------------------

message(">>> Joining RNA layers")

hipp_510 <- JoinLayers(hipp_510, assay = "RNA")
hipp_611 <- JoinLayers(hipp_611, assay = "RNA")
hipp_993 <- JoinLayers(hipp_993, assay = "RNA")

# ---------------------------------------------------------
# Ensure correct assay
# ---------------------------------------------------------

DefaultAssay(hipp_510) <- "RNA"
DefaultAssay(hipp_611) <- "RNA"
DefaultAssay(hipp_993) <- "RNA"

# ---------------------------------------------------------
# Add sample labels
# ---------------------------------------------------------

hipp_510$sample <- "510"
hipp_611$sample <- "611"
hipp_993$sample <- "993"

hipp.list <- list(hipp_510, hipp_611, hipp_993)
names(hipp.list) <- c("510", "611", "993")

# ---------------------------------------------------------
# Integration
# ---------------------------------------------------------

message(">>> Selecting features")

features <- SelectIntegrationFeatures(
  object.list = hipp.list,
  nfeatures = 2000
)

message(">>> Finding anchors")

anchors <- FindIntegrationAnchors(
  object.list = hipp.list,
  anchor.features = features,
  dims = 1:20
)

message(">>> Integrating data")

hipp_integrated <- IntegrateData(
  anchorset = anchors,
  dims = 1:20
)

# ---------------------------------------------------------
# Save
# ---------------------------------------------------------

out <- file.path(DATA_PROC, "Hipp_allSamples_integrated.rds")
saveRDS(hipp_integrated, out)

message(">>> DONE")
message(out)
