source("scripts/00_setup.R")
library(scDblFinder)
library(SingleCellExperiment)

sample_id <- "993_Hipp_B2"

h5_path  <- file.path(DATA_RAW, paste0(sample_id, ".h5"))
rds_path <- file.path(DATA_PROC, paste0(sample_id, "_qc.rds"))

stopifnot(file.exists(h5_path))

seurat.obj <- Read10X_h5(h5_path) |>
  CreateSeuratObject(
    project = sample_id,
    min.cells = 3,
    min.features = 200
  )

seurat.obj[["percent.mt"]] <- PercentageFeatureSet(
  seurat.obj,
  pattern = "^MT-"
)

seurat.obj[["percent.ribo"]] <- PercentageFeatureSet(
  seurat.obj,
  pattern = "^RPS|^RPL"
)
str(seurat.obj)

p1 <- VlnPlot(
  seurat.obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)

p2 <- FeatureScatter(
  seurat.obj,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA",
) + geom_smooth(method = "lm")


ggsave(file.path(FIG_DIR, paste0(sample_id, "_vln.png")), p1, width = 10, height = 4)
ggsave(file.path(FIG_DIR, paste0(sample_id, "_scatter.png")), p2, width = 5, height = 4)

#--------------------------------------------------
# Light pre-QC (required before DoubletFinder)
#--------------------------------------------------

seurat.obj <- subset(
  seurat.obj,
  subset =
    nFeature_RNA > 300 &
    nCount_RNA > 500 &
    percent.mt < 20
)


# seurat.obj already filtered for nFeature / nCount / percent.mt
DefaultAssay(seurat.obj) <- "RNA"

sce <- as.SingleCellExperiment(seurat.obj)

sce <- scDblFinder(
  sce,
  samples = NULL   # correct for single 10x library
)

table(colData(sce)$scDblFinder.class)
seurat.obj$scDblFinder_class <- colData(sce)$scDblFinder.class
seurat.obj$scDblFinder_score <- colData(sce)$scDblFinder.score

seurat.obj@meta.data

VlnPlot(
  seurat.obj,
  features = "nCount_RNA",
  group.by = "scDblFinder_class"
)




saveRDS(seurat.obj, rds_path)
