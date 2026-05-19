# =========================================================
# 14_UMAP_publication_figures.R
# Publication-quality UMAP figures
# =========================================================

# ---------------------------------------------------------
# LOAD SETUP
# ---------------------------------------------------------

source("00_setup.R")

# ---------------------------------------------------------
# LOAD LIBRARIES
# ---------------------------------------------------------

library(Seurat)
library(ggplot2)
library(patchwork)

# ---------------------------------------------------------
# LOAD OBJECTS
# ---------------------------------------------------------

FC <- readRDS(file.path(
  DATA_PROC,
  "FC_with_cluster_celltype_annotations.rds"
))

GE <- readRDS(file.path(
  DATA_PROC,
  "GE_with_cluster_celltype_annotations.rds"
))

Hipp <- readRDS(file.path(
  DATA_PROC,
  "Hipp_with_cluster_celltype_annotations.rds"
))

Thal <- readRDS(file.path(
  DATA_PROC,
  "Thal_with_celltype_annotations.rds"
))

Cer <- readRDS(file.path(
  DATA_PROC,
  "Cer_with_celltype_annotations.rds"
))

# ---------------------------------------------------------
# SET CELLTYPE IDENTITIES
# ---------------------------------------------------------

Idents(FC)   <- "cluster_celltype"
Idents(GE)   <- "cluster_celltype"
Idents(Hipp) <- "cluster_celltype"
Idents(Thal) <- "cluster_celltype"
Idents(Cer) <- "cell_type"
# ---------------------------------------------------------
# STANDARD UMAP FUNCTION
# ---------------------------------------------------------

plot_umap <- function(obj, title_text) {
  
  p <- DimPlot(
    obj,
    reduction = "umap",
    label = TRUE,
    repel = TRUE,
    label.size = 5
  ) +
    
    NoLegend() +
    
    ggtitle(title_text) +
    
    theme_classic() +
    
    theme(
      
      plot.title = element_text(
        face = "bold",
        size = 16,
        hjust = 0.5
      ),
      
      axis.title = element_text(
        size = 13
      ),
      
      axis.text = element_text(
        size = 11
      )
      
    )
  
  return(p)
}

# ---------------------------------------------------------
# GENERATE PANELS
# ---------------------------------------------------------

p1 <- plot_umap(FC,   "A. Frontal cortex")
p2 <- plot_umap(GE,   "B. Ganglionic eminence")
p3 <- plot_umap(Hipp, "C. Hippocampus")
p4 <- plot_umap(Thal, "D. Thalamus")
p5 <- plot_umap(Cer,  "E. Cerebellum")

# ---------------------------------------------------------
# COMBINE FIGURES
# ---------------------------------------------------------

# Figure 1
fig1 <- p1 | p2

# Figure 2
fig2 <- p3 | p4

# Figure 3
fig3 <- p5

# ---------------------------------------------------------
# SAVE FIGURE 1
# ---------------------------------------------------------

ggsave(
  filename = file.path(
    FIG_DIR,
    "Figure1_FC_GE_UMAP.pdf"
  ),
  
  plot = fig1,
  
  width = 16,
  height = 7,
  dpi = 600
)

ggsave(
  filename = file.path(
    FIG_DIR,
    "Figure1_FC_GE_UMAP.png"
  ),
  
  plot = fig1,
  
  width = 16,
  height = 7,
  dpi = 600
)

# ---------------------------------------------------------
# SAVE FIGURE 2
# ---------------------------------------------------------

ggsave(
  filename = file.path(
    FIG_DIR,
    "Figure2_Hipp_Thal_UMAP.pdf"
  ),
  
  plot = fig2,
  
  width = 16,
  height = 7,
  dpi = 600
)

ggsave(
  filename = file.path(
    FIG_DIR,
    "Figure2_Hipp_Thal_UMAP.png"
  ),
  
  plot = fig2,
  
  width = 16,
  height = 7,
  dpi = 600
)

# ---------------------------------------------------------
# SAVE FIGURE 3
# ---------------------------------------------------------

ggsave(
  filename = file.path(
    FIG_DIR,
    "Figure3_Cerebellum_UMAP.pdf"
  ),
  
  plot = fig3,
  
  width = 10,
  height = 8,
  dpi = 600
)

ggsave(
  filename = file.path(
    FIG_DIR,
    "Figure3_Cerebellum_UMAP.png"
  ),
  
  plot = fig3,
  
  width = 10,
  height = 8,
  dpi = 600
)

# ---------------------------------------------------------
# OPTIONAL SUPPLEMENTARY FIGURE
# ---------------------------------------------------------

supplementary_fig <- (p1 | p2) /
  (p3 | p4) /
  p5

ggsave(
  filename = file.path(
    FIG_DIR,
    "Supplementary_All_UMAPs.pdf"
  ),
  
  plot = supplementary_fig,
  
  width = 18,
  height = 18,
  dpi = 600
)

# ---------------------------------------------------------
# DONE
# ---------------------------------------------------------

message(">>> UMAP publication figures generated successfully")
message(">>> Saved to: ", FIG_DIR)

sessionInfo()