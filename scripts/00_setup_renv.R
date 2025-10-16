# =========================
# 00_setup_renv.R
# Initializes a per-project environment with renv and installs base deps.
# Run once from the project root.
# =========================

if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")

# 1) Initialize renv (creates renv/ and renv.lock; activates in .Rprofile)
if (is.null(renv::project())) renv::init(bare = TRUE)  # 'bare=TRUE' avoids snapshotting your global libs

# 2) Install base packages used in preprocessing & plotting
base_pkgs <- c(
  "tidyverse", "data.table", "dplyr", "ggplot2",
  "pracma", "clock", "slider"
)
renv::install(base_pkgs)

# 3) Snapshot exact versions to renv.lock
renv::snapshot(prompt = FALSE)

message("✅ renv ready. Use renv::restore() on fresh clones to reproduce the environment.")
