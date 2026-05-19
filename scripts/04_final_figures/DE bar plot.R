# ============================================================
# DE_bar_plot.R
# Schizophrenia-associated gene enrichment barplot
# TRUE hatched bars using ggpattern
# ============================================================

# ------------------------------------------------------------
# Load setup
# ------------------------------------------------------------

source("~/colemdjl-exphaem/MSc_DL_Murali/scripts/00_setup.R")

# ------------------------------------------------------------
# Install packages if needed
# ------------------------------------------------------------

# install.packages("ggpattern")
# install.packages("ggh4x")

# ------------------------------------------------------------
# Libraries
# ------------------------------------------------------------

library(ggplot2)
library(ggpattern)
library(ggh4x)
library(dplyr)
library(readr)

# ============================================================
# Load CSV
# ============================================================

input_file <- file.path(
  DATA_PROC,
  "SCZ_oRG_neuron_statistics_table.csv"
)

df_raw <- read_csv(input_file)

# ============================================================
# Create CLEAN plotting dataframe
# ============================================================

df <- data.frame(
  
  gene = c(
    "GRIN2A","GRIN2A","GRIN2A","GRIN2A",
    "KLC1","KLC1","KLC1","KLC1",
    "SP4","SP4","SP4","SP4",
    "STAG1","STAG1","STAG1","STAG1"
  ),
  
  region = c(
    "FC","GE","Hipp","Thal",
    "FC","GE","Hipp","Thal",
    "FC","GE","Hipp","Thal",
    "FC","GE","Hipp","Thal"
  ),
  
  log2FC = c(
    -2.6, 3.4, 1.5, 0.3,
    -0.8,-0.5,-0.6,-0.8,
    0.6,-0.5,0.5,-0.8,
    0.8,0.5,1.4,1.2
  ),
  
  FDR = c(
    0.001,1e-50,9.3e-28,1,
    0.015,1,3.89e-09,6.3e-05,
    1,1,5.89e-19,1.32e-13,
    0.0055,1.23e-06,5.60e-72,4.86e-125
  ),
  
  Comparison = c(
    "oRG","RGpooled","RGpooled","RGpooled",
    "RGpooled","RGpooled","RGpooled","RGpooled",
    "RGpooled","RGpooled","RGpooled","RGpooled",
    "RGpooled","RGpooled","RGpooled","RGpooled"
  )
)

# ============================================================
# Create significance stars
# ============================================================

df$Stars <- case_when(
  df$FDR < 0.001 ~ "***",
  df$FDR < 0.01  ~ "**",
  df$FDR < 0.05  ~ "*",
  TRUE ~ ""
)

# ============================================================
# Factor ordering
# ============================================================

df$gene <- factor(
  df$gene,
  levels = c("GRIN2A", "KLC1", "SP4", "STAG1")
)

df$region <- factor(
  df$region,
  levels = c("FC", "GE", "Hipp", "Thal")
)

# ============================================================
# Create plot
# ============================================================

p <- ggplot(
  df,
  aes(
    x = region,
    y = log2FC,
    fill = FDR < 0.05,
    pattern = Comparison
  )
) +
  
  # ------------------------------------------------------------
# Barplot with hatching
# ------------------------------------------------------------

geom_bar_pattern(
  stat = "identity",
  width = 0.45,
  colour = "black",
  
  pattern_fill = "black",
  pattern_colour = "black",
  
  pattern_density = 0.08,
  pattern_spacing = 0.03,
  pattern_angle = 45
) +
  
  # ------------------------------------------------------------
# Coloured facet strips
# ------------------------------------------------------------

facet_wrap2(
  ~gene,
  nrow = 1,
  
  strip = strip_themed(
    
    background_x = elem_list_rect(
      
      fill = c(
        "#f4b6b6",   # GRIN2A pink
        "#c7d9ec",   # KLC1 blue
        "#f3e6a2",   # SP4 yellow
        "#f2c7a7"    # STAG1 peach
      ),
      
      colour = "black"
    ),
    
    text_x = elem_list_text(
      face = "bold",
      size = 20
    )
  )
) +
  
  # ------------------------------------------------------------
# Zero line
# ------------------------------------------------------------

geom_hline(
  yintercept = 0,
  linetype = "dashed",
  colour = "grey50"
) +
  
  # ------------------------------------------------------------
# Significance stars
# ------------------------------------------------------------

geom_text(
  aes(
    label = Stars,
    y = ifelse(
      log2FC > 0,
      log2FC + 0.45,
      log2FC - 0.45
    )
  ),
  size = 5
) +
  
  # ------------------------------------------------------------
# Fill colours
# ------------------------------------------------------------

scale_fill_manual(
  values = c(
    "TRUE" = "#ff2d20",
    "FALSE" = "grey70"
  ),
  
  labels = c(
    "FALSE" = "Not significant",
    "TRUE" = "Significant"
  ),
  
  name = "Statistical significance (FDR)"
) +
  
  # ------------------------------------------------------------
# Pattern legend
# ------------------------------------------------------------

scale_pattern_manual(
  values = c(
    "oRG" = "none",
    "RGpooled" = "stripe"
  ),
  
  labels = c(
    "oRG" =
      "oRG vs Neurons\n(solid bar)",
    
    "RGpooled" =
      "Pooled RG (RG+oRG) vs Neurons\n(hatched bar)"
  ),
  
  name = "Comparison type"
) +
  
  # ------------------------------------------------------------
# Axis limits
# ------------------------------------------------------------

coord_cartesian(
  ylim = c(-4, 5)
) +
  
  # ------------------------------------------------------------
# Labels
# ------------------------------------------------------------

labs(
  x = "Brain region",
  
  y = expression(
    log[2]~fold~change~
      "(RG or oRG vs Neurons)"
  )
) +
  
  # ------------------------------------------------------------
# Theme
# ------------------------------------------------------------

theme_bw(base_size = 15) +
  
  theme(
    
    panel.grid = element_blank(),
    
    axis.title = element_text(
      face = "bold",
      size = 16
    ),
    
    axis.text = element_text(
      size = 13
    ),
    
    legend.position = "top",
    
    legend.box = "vertical",
    
    legend.direction = "horizontal",
    
    legend.title = element_text(
      face = "bold",
      size = 13
    ),
    
    legend.text = element_text(
      size = 12
    )
  )

# ============================================================
# Save figure
# ============================================================

output_file <- file.path(
  FIG_DIR,
  "SCZ_gene_enrichment_barplot.png"
)

ggsave(
  filename = output_file,
  plot = p,
  width = 12,
  height = 5,
  dpi = 600
)

# ============================================================
# Display plot
# ============================================================

print(p)

# ============================================================
# Session info
# ============================================================

sessionInfo()