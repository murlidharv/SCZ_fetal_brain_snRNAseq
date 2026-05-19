############################################################
# FINAL QC TABLE (Publication ready)
############################################################

source("scripts/00_setup.R")
library(tidyverse)

qc <- read_csv(file.path(DATA_PROC, "QC/QC_clean_table.csv"))

qc_pub <- qc %>%
  transmute(
    Fetus = fetus,
    Region = region,
    Rep = replicate,
    
    Cells = scales::comma(round(`Estimated Number of Cells`, 0)),
    `Reads/Cell` = scales::comma(round(`Mean Reads per Cell`, 0)),
    `Genes/Cell` = round(`Median Genes per Cell`, 0),
    
    `Saturation (%)` = round(`Sequencing Saturation`, 1),
    `Reads in Cells (%)` = round(`Fraction Reads in Cells`, 1)
  ) %>%
  arrange(Fetus, Region)

write_csv(qc_pub,
          file.path(DATA_PROC, "QC/QC_table_FINAL.csv"))