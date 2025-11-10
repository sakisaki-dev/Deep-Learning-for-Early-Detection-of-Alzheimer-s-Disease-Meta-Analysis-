# ---- PRISMA2020 flow diagram (static PNG + interactive HTML) ----
# Install + load
if (!requireNamespace("PRISMA2020", quietly = TRUE)) {
  install.packages("PRISMA2020", repos = "https://cloud.r-project.org")
}
library(PRISMA2020)

# Output dir (start fresh like your other figures)
figdir <- "alz-meta/R/generated_plots"
if (dir.exists(figdir)) unlink(figdir, recursive = TRUE, force = TRUE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

# Load the template CSV that ships with the package
csvFile <- system.file("extdata", "PRISMA.csv", package = "PRISMA2020")
tmpl    <- read.csv(csvFile, stringsAsFactors = FALSE)

# Peek available IDs if needed:
# print(unique(tmpl[[ if("data" %in% names(tmpl)) "data" else names(tmpl)[1] ]]))

# Helper to set a node count robustly
set_n <- function(df, id, n){
  idcol <- if ("data" %in% names(df)) "data" else if ("id" %in% names(df)) "id" else names(df)[1]
  if (!"n" %in% names(df)) df$n <- NA_real_
  idx <- which(df[[idcol]] == id)
  if (length(idx) == 0) {
    message("WARNING: ID not found in template: ", id)
    return(df)
  }
  suppressWarnings(df$n[idx] <- as.numeric(n))
  df
}

# ---- Your counts (from your notes) ----
# Databases: 454 (242 PubMed + 212 IEEE); Websites: 1 (arXiv)
# Duplicates removed -> 165 screened; 90 full-text assessed; 30 included.
tmpl <- set_n(tmpl, "database_results",         454)
tmpl <- set_n(tmpl, "register_results",         0)
tmpl <- set_n(tmpl, "website_results",          1)

tmpl <- set_n(tmpl, "duplicate_records",        290)
tmpl <- set_n(tmpl, "automation_tool_exclusions", 0)
tmpl <- set_n(tmpl, "other_exclusions",         0)

tmpl <- set_n(tmpl, "records_screened",         165)
tmpl <- set_n(tmpl, "records_excluded",         75)

tmpl <- set_n(tmpl, "reports_sought",           90)
tmpl <- set_n(tmpl, "reports_not_retrieved",    0)
tmpl <- set_n(tmpl, "reports_assessed",         90)

# Optional: split exclusions into reasons if you tracked them
tmpl <- set_n(tmpl, "reports_excluded",         60)  # total at full text
tmpl <- set_n(tmpl, "reports_excluded_other",   0)   # keep 0 if unused

# Included
tmpl <- set_n(tmpl, "new_studies",              30)
tmpl <- set_n(tmpl, "new_reports",              30)
tmpl <- set_n(tmpl, "total_studies",            30)
tmpl <- set_n(tmpl, "total_reports",            30)

# Convert to the list structure PRISMA2020 expects
data_list <- PRISMA_data(tmpl)

# Plot (static)
plot_static <- PRISMA_flowdiagram(
  data_list,
  interactive = FALSE,   # PNG/PDF/SVG
  previous    = FALSE,   # hide "previous studies" arm if not applicable
  other       = TRUE,    # show "other methods" (your arXiv record)
  fontsize    = 12
)
PRISMA_save(plot_static,
            filename = file.path(figdir, "prisma_flowdiagram.png"),
            filetype = "png",
            overwrite = TRUE)

# Plot (interactive HTML with tooltips/links)
plot_html <- PRISMA_flowdiagram(
  data_list,
  interactive = TRUE,
  previous    = FALSE,
  other       = TRUE,
  fontsize    = 12
)
PRISMA_save(plot_html,
            filename = file.path(figdir, "prisma_flowdiagram.html"),
            filetype = "html",
            overwrite = TRUE)

message("Saved: ", normalizePath(file.path(figdir, "prisma_flowdiagram.png")))
message("Saved: ", normalizePath(file.path(figdir, "prisma_flowdiagram.html")))
