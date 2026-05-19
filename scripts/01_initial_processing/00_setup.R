# ============================================================
# 00_setup.R
# Project-wide setup for MSc_DL_Murali
# BEAR / Slurm / RStudio safe
# ============================================================

# ------------------------------------------------------------
# 1. Project root
#    (scripts/ is one level below project root)
# ------------------------------------------------------------
PROJECT_ROOT <- normalizePath(".")

message("Project root set to: ", PROJECT_ROOT)

# ------------------------------------------------------------
# 2. Standard directory layout
# ------------------------------------------------------------
DATA_RAW  <- file.path(PROJECT_ROOT, "data", "raw")
DATA_PROC <- file.path(PROJECT_ROOT, "data", "processed")
FIG_DIR   <- file.path(PROJECT_ROOT, "figures")
LOG_DIR   <- file.path(PROJECT_ROOT, "logs")

# ------------------------------------------------------------
# 3. Create directories if missing (safe on BEAR)
# ------------------------------------------------------------
dir.create(DATA_RAW,  recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_PROC, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR,   recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR,   recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 4. Reproducibility & global options
# ------------------------------------------------------------
set.seed(1234)

options(
  stringsAsFactors = FALSE,
  scipen = 999,
  bitmapType = "cairo"
)

# ------------------------------------------------------------
# 5. Minimal sanity check (fail early if broken)
# ------------------------------------------------------------
stopifnot(
  dir.exists(DATA_PROC),
  dir.exists(FIG_DIR)
)

# ------------------------------------------------------------
# 6. Startup message
# ------------------------------------------------------------
message(">>> Setup complete")
message("    DATA_PROC: ", DATA_PROC)
message("    FIG_DIR:   ", FIG_DIR)
