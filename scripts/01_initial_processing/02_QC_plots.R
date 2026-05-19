############################################################
# 02_QC_plots.R
#
# PURPOSE:
# Generate publication-quality QC plots
############################################################

source("scripts/00_setup.R")

library(tidyverse)

############################################################
# Load QC table
############################################################

qc <- read_csv(file.path(DATA_PROC, "QC/QC_clean_table.csv"))

############################################################
# Factor ordering (important for clean plots)
############################################################

qc$region <- factor(qc$region, levels = c("FC","GE","Hipp","Thal","Cer"))

############################################################
# 1. Cells per sample
############################################################

p1 <- ggplot(qc, aes(x = sample, y = `Estimated Number of Cells`, fill = region)) +
  geom_col() +
  theme_classic() +
  labs(title = "Cells per sample", y = "Number of cells", x = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

############################################################
# 2. Genes per cell
############################################################

p2 <- ggplot(qc, aes(x = sample, y = `Median Genes per Cell`, fill = region)) +
  geom_col() +
  theme_classic() +
  labs(title = "Median genes per cell", y = "Genes", x = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

############################################################
# 3. Sequencing saturation
############################################################

p3 <- ggplot(qc, aes(x = sample, y = `Sequencing Saturation`, fill = region)) +
  geom_col() +
  theme_classic() +
  labs(title = "Sequencing saturation (%)", y = "%", x = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

############################################################
# 4. Fraction reads in cells
############################################################

p4 <- ggplot(qc, aes(x = sample, y = `Fraction Reads in Cells`, fill = region)) +
  geom_col() +
  theme_classic() +
  labs(title = "Fraction reads in cells (%)", y = "%", x = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

############################################################
# Save plots
############################################################

ggsave(file.path(FIG_DIR, "QC_cells_per_sample.pdf"), p1, width = 8, height = 4)
ggsave(file.path(FIG_DIR, "QC_genes_per_cell.pdf"), p2, width = 8, height = 4)
ggsave(file.path(FIG_DIR, "QC_saturation.pdf"), p3, width = 8, height = 4)
ggsave(file.path(FIG_DIR, "QC_fraction_reads.pdf"), p4, width = 8, height = 4)

############################################################
# DONE
############################################################

cat("QC plots saved in figures/\n")