# =========================================================
# Thal annotation using published top-20 markers
# =========================================================

source("scripts/00_setup.R")
library(Seurat)
library(dplyr)

message(">>> Loading integrated Thal object")

thal <- readRDS(file.path(DATA_PROC, "Thal_allSamples_integrated.rds"))
DefaultAssay(thal) <- "RNA"

# ---------------------------------------------------------
# Load published markers
# ---------------------------------------------------------

marker_file <- file.path(DATA_PROC, "Thal_top20_cluster_markers.csv")
markers <- read.csv(marker_file, stringsAsFactors = FALSE)

# ---------------------------------------------------------
# Clean marker table
# ---------------------------------------------------------

# Drop the junk header row
markers <- markers %>%
  filter(cluster != "Cluster ID")

# Basic sanity check
stopifnot(all(c("cluster", "gene") %in% colnames(markers)))

message(">>> Marker clusters found:")
print(unique(markers$cluster))

# ---------------------------------------------------------
# Build marker list for AddModuleScore
# ---------------------------------------------------------

marker_list <- markers %>%
  group_by(cluster) %>%
  summarise(
    genes = list(unique(gene)),
    .groups = "drop"
  ) %>%
  deframe()

message(">>> Number of markers per cluster:")
print(lapply(marker_list, length))

# ---------------------------------------------------------
# Add module scores
# ---------------------------------------------------------

message(">>> Adding module scores")

thal <- AddModuleScore(
  object = thal,
  features = marker_list,
  assay = "RNA",
  name = "Thal_marker"
)

score_cols <- grep("^Thal_marker", colnames(thal@meta.data), value = TRUE)

# ---------------------------------------------------------
# Assign Thal cell-type labels
# ---------------------------------------------------------

message(">>> Assigning cell-type labels")

score_mat <- thal@meta.data[, score_cols, drop = FALSE]

thal$Thal_celltype <- names(marker_list)[
  apply(score_mat, 1, which.max)
]

# ---------------------------------------------------------
# Sanity checks
# ---------------------------------------------------------

message(">>> Thal cell-type distribution:")
print(table(thal$Thal_celltype))

# ---------------------------------------------------------
# Save annotated object
# ---------------------------------------------------------

out <- file.path(DATA_PROC, "Thal_integrated_annotated.rds")
saveRDS(thal, out)

message(">>> DONE")
message(out)
