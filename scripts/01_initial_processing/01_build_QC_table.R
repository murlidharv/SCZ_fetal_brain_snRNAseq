############################################################
# 01_build_QC_table.R  (FINAL WORKING VERSION)
############################################################

source("scripts/00_setup.R")

library(tidyverse)

############################################################
# 1. Paths
############################################################

INPUT_DIR  <- DATA_RAW
OUTPUT_DIR <- file.path(DATA_PROC, "QC")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

############################################################
# 2. List files
############################################################

files <- list.files(
  path = INPUT_DIR,
  pattern = "_metrics_summary\\.csv$",
  full.names = TRUE
)

############################################################
# 3. SAFE reader (NO pivot_longer)
############################################################

read_metrics <- function(file_path) {
  
  df <- read_csv(file_path, col_types = cols(.default = "c"))
  
  sample_name <- basename(file_path) %>%
    str_remove("_metrics_summary\\.csv")
  
  df$sample <- sample_name
  
  return(df)
}

############################################################
# 4. Combine
############################################################

qc_wide <- purrr::map_dfr(files, read_metrics)

############################################################
# 5. Extract metadata
############################################################

qc_wide <- qc_wide %>%
  separate(sample, into = c("fetus", "region", "replicate"), sep = "_", remove = FALSE)

############################################################
# 6. Save full table
############################################################

write_csv(
  qc_wide,
  file.path(OUTPUT_DIR, "QC_full_table.csv")
)

############################################################
# 7. Clean numeric columns
############################################################

clean_numeric <- function(x) {
  x %>%
    str_remove_all(",") %>%
    str_remove("%") %>%
    as.numeric()
}

qc_clean <- qc_wide %>%
  mutate(
    `Estimated Number of Cells` = clean_numeric(`Estimated Number of Cells`),
    `Mean Reads per Cell` = clean_numeric(`Mean Reads per Cell`),
    `Median Genes per Cell` = clean_numeric(`Median Genes per Cell`),
    `Total Genes Detected` = clean_numeric(`Total Genes Detected`),
    `Number of Reads` = clean_numeric(`Number of Reads`),
    `Valid Barcodes` = clean_numeric(`Valid Barcodes`),
    `Sequencing Saturation` = clean_numeric(`Sequencing Saturation`),
    `Q30 Bases in Barcode` = clean_numeric(`Q30 Bases in Barcode`),
    `Q30 Bases in RNA Read` = clean_numeric(`Q30 Bases in RNA Read`),
    `Reads Mapped Confidently to Genome` = clean_numeric(`Reads Mapped Confidently to Genome`),
    `Reads Mapped Confidently to Transcriptome` = clean_numeric(`Reads Mapped Confidently to Transcriptome`),
    `Fraction Reads in Cells` = clean_numeric(`Fraction Reads in Cells`)
  ) %>%
  select(sample, fetus, region, replicate, everything())

############################################################
# 8. Save clean table
############################################################

write_csv(
  qc_clean,
  file.path(OUTPUT_DIR, "QC_clean_table.csv")
)

############################################################
# DONE
############################################################

cat("\n✅ QC TABLE CREATED SUCCESSFULLY\n")
cat("Samples:", nrow(qc_clean), "\n")