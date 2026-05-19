# =========================================================
# GE annotation using published top-20 markers
# =========================================================

source("scripts/00_setup.R")
library(Seurat)
library(dplyr)

message(">>> Loading integrated GE object")

ge <- readRDS(
  file.path(DATA_PROC, "GE_integrated_clustered.rds")
)
DefaultAssay(ge) <- "RNA"

# ---------------------------------------------------------
# Load published markers
# ---------------------------------------------------------

marker_file <- file.path(
  DATA_PROC,
  "GE_top20_cluster_markers.csv"
)

markers <- read.csv(
  marker_file,
  stringsAsFactors = FALSE
)

# ---------------------------------------------------------
# Clean marker table
# ---------------------------------------------------------

# Drop junk header row if present
markers <- markers %>%
  filter(cluster != "Cluster ID")

# Basic sanity check
stopifnot(
  all(c("cluster", "gene") %in% colnames(markers))
)

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

ge <- AddModuleScore(
  object  = ge,
  features = marker_list,
  assay   = "RNA",
  name    = "GE_marker"
)

score_cols <- grep(
  "^GE_marker",
  colnames(ge@meta.data),
  value = TRUE
)

stopifnot(length(score_cols) > 0)

# ---------------------------------------------------------
# Assign GE cell-type labels
# ---------------------------------------------------------

message(">>> Assigning cell-type labels")

score_mat <- ge@meta.data[, score_cols, drop = FALSE]

ge$GE_celltype <- names(marker_list)[
  apply(score_mat, 1, which.max)
]

# ---------------------------------------------------------
# Sanity checks
# ---------------------------------------------------------

message(">>> GE cell-type distribution:")
print(table(ge$GE_celltype))

# ---------------------------------------------------------
# Save annotated object
# ---------------------------------------------------------

out <- file.path(
  DATA_PROC,
  "GE_integrated_annotated.rds"
)

saveRDS(ge, out)

message(">>> DONE")
message(out)
