############################################################
# 07_integrate_Thal_allSamples.R
# Purpose:
# - Integrate normalized Thalamus samples (510, 611, 993)
# - Seurat v5 compatible (JoinLayers)
############################################################

source("scripts/00_setup.R")
library(Seurat)

message(">>> Loading normalized Thal objects")

thal_510 <- readRDS(file.path(DATA_PROC, "510_Thal_B1_normalized.rds"))
thal_611 <- readRDS(file.path(DATA_PROC, "611_Thal_B2_normalized.rds"))
thal_993 <- readRDS(file.path(DATA_PROC, "993_Thal_B2_normalized.rds"))

# ---------------------------------------------------------
# Seurat v5: collapse RNA layers
# ---------------------------------------------------------

message(">>> Joining RNA layers")

thal_510 <- JoinLayers(thal_510, assay = "RNA")
thal_611 <- JoinLayers(thal_611, assay = "RNA")
thal_993 <- JoinLayers(thal_993, assay = "RNA")

# ---------------------------------------------------------
# Ensure correct assay
# ---------------------------------------------------------

DefaultAssay(thal_510) <- "RNA"
DefaultAssay(thal_611) <- "RNA"
DefaultAssay(thal_993) <- "RNA"

# ---------------------------------------------------------
# Add sample labels
# ---------------------------------------------------------

thal_510$sample <- "510"
thal_611$sample <- "611"
thal_993$sample <- "993"

thal.list <- list(thal_510, thal_611, thal_993)
names(thal.list) <- c("510", "611", "993")

# ---------------------------------------------------------
# Integration
# ---------------------------------------------------------

message(">>> Selecting integration features")
features <- SelectIntegrationFeatures(
  object.list = thal.list,
  nfeatures = 2000
)

message(">>> Finding integration anchors")
anchors <- FindIntegrationAnchors(
  object.list = thal.list,
  anchor.features = features,
  dims = 1:20
)

message(">>> Integrating data")
thal_integrated <- IntegrateData(
  anchorset = anchors,
  dims = 1:20
)

# ---------------------------------------------------------
# Save
# ---------------------------------------------------------

out <- file.path(DATA_PROC, "Thal_allSamples_integrated.rds")
saveRDS(thal_integrated, out)

message(">>> DONE")
message(out)
