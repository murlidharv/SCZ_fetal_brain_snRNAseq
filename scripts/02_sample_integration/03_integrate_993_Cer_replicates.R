## -------------------------------------------------
## Anchor-based integration of biological replicates
## Sample: 993 | Region: Cerebellum
## Replicates: B1, B2
## -------------------------------------------------

source("scripts/00_setup.R")
library(Seurat)

## -------------------------------
## Sample IDs
## -------------------------------

sample_ids <- c("993_Cer_B1", "993_Cer_B2")

rds_in <- file.path(
  DATA_PROC,
  paste0(sample_ids, "_filtered.rds")
)

rds_out <- file.path(
  DATA_PROC,
  "993_Cer_integrated.rds"
)

stopifnot(all(file.exists(rds_in)))

## -------------------------------
## Load filtered objects
## -------------------------------

cer_b1 <- readRDS(rds_in[1])
cer_b2 <- readRDS(rds_in[2])

## -------------------------------
## Add metadata
## -------------------------------

cer_b1$replicate <- "B1"
cer_b2$replicate <- "B2"

cer_b1$sample <- "993"
cer_b2$sample <- "993"

cer_b1$region <- "Cer"
cer_b2$region <- "Cer"

## -------------------------------
## Create object list
## -------------------------------

cer.list <- list(
  Cer_B1 = cer_b1,
  Cer_B2 = cer_b2
)

## -------------------------------
## Normalize & find HVGs separately
## (required for anchor finding)
## -------------------------------

cer.list <- lapply(cer.list, function(obj) {
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(
    obj,
    selection.method = "vst",
    nfeatures = 2000,
    verbose = FALSE
  )
  obj
})

## -------------------------------
## Find integration anchors
## -------------------------------

cer.anchors <- FindIntegrationAnchors(
  object.list = cer.list,
  dims = 1:30
)

## -------------------------------
## Integrate data
## -------------------------------

cer_993 <- IntegrateData(
  anchorset = cer.anchors,
  dims = 1:30
)

## -------------------------------
## Save integrated object
## -------------------------------

saveRDS(cer_993, rds_out)

## -------------------------------------------------
## End of script
## -------------------------------------------------
#options(bitmapType = "cairo")

#table(cer_993$replicate)
#p <- VlnPlot(cer_993, features = "nCount_RNA", group.by = "replicate")
#p
