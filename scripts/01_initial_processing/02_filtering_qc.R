source("scripts/00_setup.R")

sample_id <- "993_Hipp_B2"

rds_in  <- file.path(DATA_PROC, paste0(sample_id, "_qc.rds"))
rds_out <- file.path(DATA_PROC, paste0(sample_id, "_filtered.rds"))

stopifnot(file.exists(rds_in))

seurat.obj <- readRDS(rds_in)

table(seurat.obj$scDblFinder_class)

seurat.obj <- subset(
  seurat.obj,
  subset =
    nFeature_RNA >= 1000 &
    nFeature_RNA <= 5000 &
    percent.mt <= 5 &
    scDblFinder_class == "singlet"&
    percent.ribo <= 10
  )
table(seurat.obj$scDblFinder_class)
saveRDS(seurat.obj, rds_out)

# ===============================
# Post-filter QC plot
# ===============================
VlnPlot(
  seurat.obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)

