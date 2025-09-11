packages <- c("metafor","mada","dmetar","readr","dplyr","tidyr","ggplot2")
to_install <- setdiff(packages, rownames(installed.packages()))
if(length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
# Snapshot with renv if you want
if (!requireNamespace("renv", quietly=TRUE)) install.packages("renv")
# renv::init() # run interactively if you want an R lockfile
