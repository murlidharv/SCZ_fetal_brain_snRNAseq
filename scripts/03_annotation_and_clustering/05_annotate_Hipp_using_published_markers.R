# =========================================================
# Hipp annotation using published top-20 markers
# =========================================================

source("scripts/00_setup.R")
library(Seurat)
library(dplyr)

message(">>> Loading integrated Hipp object")

hipp <- readRDS(
  file.path(DATA_PROC, "Hipp_integrated_clustered.rds")
)

DefaultAssay(hipp) <- "RNA"

# ---------------------------------------------------------
# Load published markers
# ---------------------------------------------------------

marker_file <- file.path(
  DATA_PROC,
  "Hipp_top20_cluster_markers.csv"
)

markers <- read.csv(
  marker_file,
  stringsAsFactors = FALSE
)

# ---------------------------------------------------------
# Clean marker table
# ---------------------------------------------------------

markers <- markers %>%
  filter(cluster != "Cluster ID")

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

hipp <- AddModuleScore(
  object   = hipp,
  features = marker_list,
  assay    = "RNA",
  name     = "Hipp_marker"
)

score_cols <- grep(
  "^Hipp_marker",
  colnames(hipp@meta.data),
  value = TRUE
)

stopifnot(length(score_cols) > 0)

# ---------------------------------------------------------
# Assign Hipp cell-type labels
# ---------------------------------------------------------

message(">>> Assigning cell-type labels")

score_mat <- hipp@meta.data[, score_cols, drop = FALSE]

hipp$Hipp_celltype <- names(marker_list)[
  apply(score_mat, 1, which.max)
]

# ---------------------------------------------------------
# Sanity checks
# ---------------------------------------------------------

message(">>> Hipp cell-type distribution:")
print(table(hipp$Hipp_celltype))

# ---------------------------------------------------------
# Save annotated object
# ---------------------------------------------------------

out <- file.path(
  DATA_PROC,
  "Hipp_integrated_annotated.rds"
)

saveRDS(hipp, out)

message(">>> DONE")
message(out)
