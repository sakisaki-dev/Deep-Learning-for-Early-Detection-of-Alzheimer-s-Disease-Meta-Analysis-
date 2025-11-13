options(repos = c(CRAN = "https://cloud.r-project.org"))

cran_packages <- c("metafor", "mada", "readr", "dplyr", "tidyr", "ggplot2")

to_install <- setdiff(cran_packages, rownames(installed.packages()))
if(length(to_install)) {
  install.packages(to_install)
}

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

# Check versions of installed package
lapply(c(cran_packages, "dmetar"), function(pkg) {
  if(requireNamespace(pkg, quietly = TRUE)) {
    cat(pkg, "version:", as.character(packageVersion(pkg)), "\n")
  } else {
    cat(pkg, "is not installed\n")
  }
})

