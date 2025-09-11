# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Packages to install from CRAN
cran_packages <- c("metafor", "mada", "readr", "dplyr", "tidyr", "ggplot2")

# Check for and install missing CRAN packages
to_install <- setdiff(cran_packages, rownames(installed.packages()))
if(length(to_install)) {
  install.packages(to_install)
}

# Install GitHub package for dmetar if not already installed
if (!requireNamespace("dmetar", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  remotes::install_github("MathiasHarrer/dmetar")
}

# Optional: renv snapshot for reproducibility
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}
# Uncomment the next line to initialize renv lockfile interactively
# renv::init()

# Check versions of installed packages
lapply(c(cran_packages, "dmetar"), function(pkg) {
  if(requireNamespace(pkg, quietly = TRUE)) {
    cat(pkg, "version:", as.character(packageVersion(pkg)), "\n")
  } else {
    cat(pkg, "is not installed\n")
  }
})

