# ==============================================================================
# ===== (2026) Drosophilidae 16S Virus Infection: Dependencies Setup ===========
# ==============================================================================

# ------------------------------------------------------------------------------
# ----- 0. Initialisation ------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 0.1. Description -------------------------------------------------------

# The following script installs all dependencies used in this study.

# Note: This script will install the most recent versions of each dependency.
#       The exact versions of each tool used in this study are provided below.
#       It may be necessary to install exact versions if significant changes
#       have been made to a library since this study was published.

# ------------------------------------------------------------------------------
# ----- 1. Install Packages ----------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 1.1. Specify Packages --------------------------------------------------

packages <- data.frame(
  package = c("tidyverse", "MCMCglmm", "ape",   "minpack.lm", "here", "viridisLite", "scales", "patchwork", "ggtree"),
  version = c("2.0.0",     "2.36",     "5.8.1", "1.2.4",     "1.0.2", "0.4.2",       "1.4.0",  "1.3.0",     "3.16.0"),
  source = c("cran",       "cran",     "cran",  "cran",       "cran",  "cran",       "cran",   "cran",      "bioc"))

# ----- 1.2. Install Functions -------------------------------------------------

install_cran <- function(package) {
    install.packages(package, dependencies = TRUE)
}

install_bioc <- function(package) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}
  BiocManager::install(package, ask = FALSE, update = FALSE)
}

install_github <- function(repo) {
  if (!requireNamespace("remotes", quietly = TRUE)) {install.packages("remotes")}
  remotes::install_github(repo, dependencies = TRUE, upgrade = "never")
}


# ----- 1.3. Install Packages --------------------------------------------------

for (i in c(1:nrow(packages))) {
  
  package <- packages$package[i]
  source <- packages$source[i]
  
  if (requireNamespace(package, quietly = TRUE)) {
    print(sprintf("%s already installed, skipping.", package))
    next
  }
  
  print(sprintf("Installing %s", package))
  
  if (source == "cran") {
    install_cran(package)
  } else if (source == "bioc") {
    install_bioc(package)
  } else {
    install_github(package)
  }
  
}
