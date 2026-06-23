# ------------------------------------------------------------------------------
# File:     scripts/setup_r.R
# Purpose:  One-time bootstrap. Installs the required CRAN stack and snapshots
#           versions to renv.lock so collaborators get the same environment.
#
# Usage (from project root):
#           Rscript scripts/setup_r.R
# ------------------------------------------------------------------------------

if (getRversion() < "4.3.0") stop("Requires R >= 4.3.0; you have ", R.version.string)

REQUIRED_PKGS <- c(
  # --- core infrastructure ---------------------------------------------------
  "renv",                         # environment locking (must be first)
  "tidyverse", "haven",           # data wrangling + .dta import
  "here", "fs", "glue",           # paths + small utilities
  "log4r",                        # structured logging

  # --- regression core -------------------------------------------------------
  "fixest",                       # high-dim FE regression + clustered SEs (OLS, FE, IV, Sun-Abraham)
  "sandwich", "lmtest",           # alternative HC robust SE machinery
  "estimatr",                     # lm_robust / iv_robust (HC2/CR2 by default)
  "AER", "ivmodel",               # IV diagnostics + Anderson-Rubin CIs
  "fwildclusterboot",             # wild-cluster bootstrap for G < ~30
  "survey",                       # complex survey weights / svyglm

  # --- staggered / heterogeneity-robust DiD ----------------------------------
  "did",                          # Callaway-Sant'Anna (att_gt)
  "did2s",                        # Borusyak-Jaravel-Spiess two-stage DiD
  "DIDmultiplegt",                # de Chaisemartin-D'Haultfoeuille
  "staggered",                    # Roth-Sant'Anna efficient estimator
  "HonestDiD",                    # parallel-trends sensitivity bounds

  # --- double / debiased machine learning ------------------------------------
  "ddml",                         # DDML for partial linear / IV models (Ahrens et al.)
  "glmnet", "ranger", "xgboost",  # default first-stage learners for ddml

  # --- survival analysis -----------------------------------------------------
  "survival",                     # Cox PH (coxph), Surv(), survfit() — also ships with R

  # --- publication output ----------------------------------------------------
  "modelsummary", "kableExtra",   # tables (.tex, .csv, .html)
  "ggplot2", "patchwork",         # figures + multi-panel composition
  "scales", "broom"               # axis breaks + tidy(model)
)

cat("== R Research Pipeline: setup ==\n")
cat("R version:  ", R.version.string, "\n")
cat("Library:    ", .libPaths()[[1]], "\n\n")

missing <- REQUIRED_PKGS[!vapply(REQUIRED_PKGS, requireNamespace,
                                 logical(1), quietly = TRUE)]

if (length(missing) > 0) {
  cat("Installing:", paste(missing, collapse = ", "), "\n")
  install.packages(missing, repos = c(CRAN = "https://cloud.r-project.org"))
} else {
  cat("All required packages already installed.\n")
}

# Initialize renv if not already initialized; otherwise snapshot the current
# library state so the lockfile stays in sync.
if (!file.exists("renv.lock")) {
  cat("\nInitializing renv (creating renv.lock + renv/activate.R) ...\n")
  renv::init(bare = TRUE, restart = FALSE)
} else {
  cat("\nrenv.lock present; snapshotting current library state ...\n")
}
renv::snapshot(prompt = FALSE)

cat("\n== Setup complete ==\n")
cat("Next:  bash scripts/run_pipeline.sh\n")
