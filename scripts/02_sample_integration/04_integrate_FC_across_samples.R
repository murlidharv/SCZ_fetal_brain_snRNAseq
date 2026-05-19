source("scripts/00_setup.R")
library(Seurat)

message(">>> Loading normalized FC objects")

fc_510 <- readRDS(file.path(DATA_PROC, "510_FC_B1_normalized.rds"))
fc_611 <- readRDS(file.path(DATA_PROC, "611_FC_B1_normalized.rds"))
fc_993 <- readRDS(file.path(DATA_PROC, "993_FC_B2_normalized.rds"))

# ---------------------------------------------------------
# Seurat v5: collapse RNA layers
# ---------------------------------------------------------

message(">>> Joining RNA layers")

fc_510 <- JoinLayers(fc_510, assay = "RNA")
fc_611 <- JoinLayers(fc_611, assay = "RNA")
fc_993 <- JoinLayers(fc_993, assay = "RNA")

# ---------------------------------------------------------
# Ensure correct assay
# ---------------------------------------------------------

DefaultAssay(fc_510) <- "RNA"
DefaultAssay(fc_611) <- "RNA"
DefaultAssay(fc_993) <- "RNA"

# ---------------------------------------------------------
# Add sample labels
# ---------------------------------------------------------

fc_510$sample <- "510"
fc_611$sample <- "611"
fc_993$sample <- "993"

fc.list <- list(fc_510, fc_611, fc_993)
names(fc.list) <- c("510", "611", "993")

# ---------------------------------------------------------
# Integration
# ---------------------------------------------------------

message(">>> Selecting features")
features <- SelectIntegrationFeatures(
  object.list = fc.list,
  nfeatures = 2000
)

message(">>> Finding anchors")
anchors <- FindIntegrationAnchors(
  object.list = fc.list,
  anchor.features = features,
  dims = 1:20
)

message(">>> Integrating data")
fc_integrated <- IntegrateData(
  anchorset = anchors,
  dims = 1:20
)

# ---------------------------------------------------------
# Save
# ---------------------------------------------------------

out <- file.path(DATA_PROC, "FC_allSamples_integrated.rds")
saveRDS(fc_integrated, out)

message(">>> DONE")
message(out)
