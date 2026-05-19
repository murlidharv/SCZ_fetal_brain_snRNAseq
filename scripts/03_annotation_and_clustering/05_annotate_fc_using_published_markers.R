# =========================================================
# FC annotation using published top-20 markers
# =========================================================

source("scripts/00_setup.R")
library(Seurat)
library(dplyr)

message(">>> Loading integrated FC object")

fc <- readRDS(file.path(DATA_PROC, "FC_allSamples_integrated.rds"))
DefaultAssay(fc) <- "RNA"

# ---------------------------------------------------------
# Load published markers
# ---------------------------------------------------------

marker_file <- file.path(DATA_PROC, "FC_top20_cluster_markers.csv")
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

fc <- AddModuleScore(
  object = fc,
  features = marker_list,
  assay = "RNA",
  name = "FC_marker"
)

score_cols <- grep("^FC_marker", colnames(fc@meta.data), value = TRUE)

# ---------------------------------------------------------
# Assign FC cell-type labels
# ---------------------------------------------------------

message(">>> Assigning cell-type labels")

score_mat <- fc@meta.data[, score_cols, drop = FALSE]

fc$FC_celltype <- names(marker_list)[
  apply(score_mat, 1, which.max)
]

# ---------------------------------------------------------
# Sanity checks
# ---------------------------------------------------------

message(">>> FC cell-type distribution:")
print(table(fc$FC_celltype))



# ---------------------------------------------------------
# Save annotated object
# ---------------------------------------------------------

out <- file.path(DATA_PROC, "FC_integrated_annotated.rds")
saveRDS(fc, out)

message(">>> DONE")
message(out)




