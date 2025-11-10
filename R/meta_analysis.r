# ===============================
# Deep Learning for AD — Analysis (Zoomed Scales)
# ===============================

# ---- 0) Packages ----
req <- c("tidyverse","readr","janitor","stringr","scales","viridis",
         "metafor","GGally")
new <- setdiff(req, rownames(installed.packages()))
if(length(new)) install.packages(new, repos="https://cloud.r-project.org")
suppressPackageStartupMessages({
  library(tidyverse); library(readr); library(janitor); library(stringr)
  library(scales); library(viridis); library(metafor); library(GGally)
})

# ---- 1) IO ----
in_csv <- "alz-meta/data_clean/data_extracted_clean.csv"
figdir <- "alz-meta/R/generated_plots"

# Always start fresh (your request)
if (dir.exists(figdir)) unlink(figdir, recursive = TRUE, force = TRUE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

# ---- 2) Load + clean ----
raw <- read_csv(in_csv, show_col_types = FALSE)
df  <- raw %>% clean_names()

# Helpers

# Backward/forward compatible SI labeler
safe_label_si <- function(accuracy = 1) {
  if (utils::packageVersion("scales") >= "1.2.0") {
    scales::label_number(accuracy = accuracy, scale_cut = scales::cut_si(" "))
  } else if ("label_number_si" %in% getNamespaceExports("scales")) {
    scales::label_number_si(accuracy = accuracy)
  } else {
    scales::label_number(accuracy = accuracy)
  }
}

nr2na <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("NR","NaN","nan","NA","", "N/A")] <- NA
  x
}
num_safely <- function(x) suppressWarnings(readr::parse_number(nr2na(x)))
yn_norm <- function(x){
  z <- tolower(nr2na(x))
  dplyr::case_when(
    z %in% c("yes","y","true","1") ~ "Yes",
    z %in% c("no","n","false","0") ~ "No",
    z %in% c("partial","partially","some") ~ "Partial",
    TRUE ~ ifelse(is.na(z),"", str_to_title(z))
  )
}
mod_norm <- function(x){
  z <- tolower(nr2na(x)); z <- str_replace_all(z, fixed(" "), "")
  dplyr::case_when(
    str_detect(z,"multimodal|multi") ~ "Multimodal",
    str_detect(z,"eeg") ~ "EEG",
    str_detect(z,"pet") ~ "PET",
    str_detect(z,"fmri") ~ "fMRI",
    str_detect(z,"retina|fundus") ~ "Retinal/Fundus",
    str_detect(z,"mri") ~ "MRI",
    TRUE ~ str_to_title(nr2na(x))
  )
}
model_type_norm <- function(x){
  z <- str_squish(tolower(nr2na(x)))
  dplyr::case_when(
    str_detect(z,"ensemble") ~ "Ensemble",
    str_detect(z,"transformer") ~ "Transformer",
    str_detect(z,"gnn|graph") ~ "GNN",
    str_detect(z,"cnn|conv") ~ "CNN",
    str_detect(z,"deep") ~ "Deep Learning",
    str_detect(z,"svm|rf|qda|lda|knn|nb|logistic|ridge|ridg") ~ "Classical ML",
    TRUE ~ str_to_title(nr2na(x))
  )
}

df <- df %>%
  mutate(
    modality            = mod_norm(modality),
    model_type          = model_type_norm(model_type),
    external_validation = yn_norm(external_validation),
    sample_size_total   = num_safely(sample_size_total),
    auc_roc             = num_safely(auc_roc),
    sensitivity         = num_safely(sensitivity_recall),
    specificity         = num_safely(specificity),
    f1_precision_num    = num_safely(f1_precision),
    year                = str_extract(author_year_institution, "(19|20)\\d{2}") %>% as.numeric()
  )

message("N rows: ", nrow(df))
message("AUC available in: ", sum(!is.na(df$auc_roc)), " studies")

# ---- 2b) Smart zoom helpers (so clustered values are visible) ----
compute_zoom <- function(x, pad = 0.005, min_span = 0.06, hard = c(0.80, 1.00)) {
  x <- x[is.finite(x)]
  x <- x[!is.na(x)]
  if (!length(x)) return(hard)
  lo <- max(hard[1], quantile(x, 0.05, na.rm = TRUE) - pad)
  hi <- min(hard[2], quantile(x, 0.95, na.rm = TRUE) + pad)
  if ((hi - lo) < min_span) {
    mid <- median(x, na.rm = TRUE)
    lo <- max(hard[1], mid - min_span/2)
    hi <- min(hard[2], mid + min_span/2)
  }
  c(lo, hi)
}
auc_zoom <- compute_zoom(df$auc_roc)  # dynamic AUC axis
ss_zoom  <- compute_zoom(df$sample_size_total / max(df$sample_size_total, na.rm=TRUE)) # not used directly

# ---- 3) Meta-analysis of AUC (random-effects, logit) ----
dat_meta <- df %>%
  filter(!is.na(auc_roc), !is.na(sample_size_total), sample_size_total > 3) %>%
  mutate(
    auc_clip    = pmin(pmax(auc_roc, 0.001), 0.999),  # avoid 0/1
    yi          = qlogis(auc_clip),                   # logit(AUC)
    vi          = 1 / (sample_size_total * auc_clip * (1 - auc_clip)),
    first_author= str_trim(str_extract(author_year_institution, "^[^,]+")),
    short_label = dplyr::coalesce(
      if_else(!is.na(first_author) & !is.na(year),
              paste0(first_author, " ", year, " (", study_id, ")"),
              NA_character_),
      study_id, first_author, ""
    ),
    short_label = stringr::str_trunc(short_label, 28)
  ) %>% arrange(auc_clip)

res <- tryCatch(rma(yi = yi, vi = vi, data = dat_meta, method = "REML"),
                error = function(e) NULL)

if (!is.null(res)) {
  # (A) Metafor forest (AUC axis zoomed)
  png(file.path(figdir, "forest_auc.png"), width = 1400, height = 1600, res = 220)
  op <- par(mar = c(4,4,1,1))
  forest(res,
         slab    = dat_meta$short_label,
         transf  = transf.ilogit,           # show AUC
         xlab    = "AUC",
         refline = 0.5,
         alim    = auc_zoom,                # << zoom the x-axis
         at      = seq(auc_zoom[1], auc_zoom[2], by = 0.02),
         cex     = 0.72)
  par(op); dev.off()
  
  # (B) Compact ggplot forest (zoomed)
  studies <- dat_meta %>%
    mutate(
      auc_ci_lo = plogis(yi - 1.96*sqrt(vi)),
      auc_ci_hi = plogis(yi + 1.96*sqrt(vi)),
      short_label = forcats::fct_reorder(short_label, auc_clip)
    )
  p_forest <- ggplot(studies, aes(x = auc_clip, y = short_label)) +
    geom_errorbarh(aes(xmin = auc_ci_lo, xmax = auc_ci_hi), height = 0) +
    geom_point(size = 2) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey60") +
    coord_cartesian(xlim = auc_zoom) +                     # << zoom
    scale_x_continuous(labels = number_format(accuracy = 0.001)) +
    theme_minimal(base_size = 12) +
    labs(title = "Diagnostic Accuracy (AUC) — Random-Effects Model",
         x = "AUC", y = NULL)
  ggsave(file.path(figdir, "forest_auc_compact.png"),
         p_forest, width = 9, height = 10, dpi = 220)
  
  # Funnel (still on logit scale)
  png(file.path(figdir, "funnel_plot.png"), width = 1200, height = 1000, res = 180)
  funnel(res, xlab = "Logit(AUC)", ylab = "Standard Error")
  dev.off()
}

# ---- 4) AUC by Model Type (boxplot, zoomed) ----
p1 <- df %>%
  filter(!is.na(auc_roc), !is.na(model_type)) %>%
  mutate(model_type = fct_reorder(model_type, auc_roc, .fun = median, na.rm = TRUE)) %>%
  ggplot(aes(x = model_type, y = auc_roc, fill = model_type)) +
  geom_boxplot(outlier.alpha = 0.4, width = 0.7) +
  coord_cartesian(ylim = auc_zoom) +                      # << zoom
  scale_y_continuous(breaks = seq(auc_zoom[1], auc_zoom[2], by = 0.01),
                     labels = number_format(accuracy = 0.001)) +
  theme_minimal(base_size = 12) + theme(legend.position = "none") +
  labs(title = "AUC by Model Type (Zoomed)", x = "Model Type", y = "AUC")
ggsave(file.path(figdir, "auc_by_model_type.png"), p1, width = 10, height = 6, dpi = 200)

# ---- 5) AUC by Modality (violin + box, zoomed) ----
p2 <- df %>%
  filter(!is.na(auc_roc), !is.na(modality)) %>%
  mutate(modality = fct_reorder(modality, auc_roc, .fun = median, na.rm = TRUE)) %>%
  ggplot(aes(x = modality, y = auc_roc, fill = modality)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.18, outlier.alpha = 0.3) +
  coord_cartesian(ylim = auc_zoom) +                      # << zoom
  scale_y_continuous(breaks = seq(auc_zoom[1], auc_zoom[2], by = 0.01),
                     labels = number_format(accuracy = 0.001)) +
  theme_minimal(base_size = 12) + theme(legend.position = "none") +
  labs(title = "Performance by Imaging Modality (Zoomed)", x = "Modality", y = "AUC")
ggsave(file.path(figdir, "auc_by_modality.png"), p2, width = 10, height = 6, dpi = 200)

# ---- 6) Sample Size vs AUC (log10 X, zoomed Y) ----
p3 <- df %>%
  filter(!is.na(auc_roc), !is.na(sample_size_total), sample_size_total > 0) %>%
  ggplot(aes(x = sample_size_total, y = auc_roc, color = modality)) +
  geom_point(alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  scale_x_log10(labels = safe_label_si(accuracy = 1)) + # << log scale
  coord_cartesian(ylim = auc_zoom) +                      # << zoom
  theme_minimal(base_size = 12) +
  labs(title = "Relationship Between Sample Size (log10) and AUC",
       x = "Sample Size (log10)", y = "AUC", color = "Modality")
ggsave(file.path(figdir, "sample_vs_auc.png"), p3, width = 10, height = 6, dpi = 200)

# ---- 7) Heatmap: mean AUC by Model x Modality (centered on global mean) ----
global_mean <- mean(df$auc_roc, na.rm = TRUE)
heat_data <- df %>%
  filter(!is.na(auc_roc), !is.na(model_type), !is.na(modality)) %>%
  group_by(model_type, modality) %>%
  summarise(mean_auc = mean(auc_roc, na.rm = TRUE), n = n(), .groups = "drop")

p4 <- ggplot(heat_data, aes(x = model_type, y = modality, fill = mean_auc)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.3f\n(n=%d)", mean_auc, n)), size = 3) +
  scale_fill_viridis(limits = c(auc_zoom[1], auc_zoom[2]),    # << same zoomed scale
                     option = "C", direction = 1) +
  theme_minimal(base_size = 12) +
  labs(title = paste0("Average AUC by Model Type and Modality (Mean≈", sprintf("%.3f", global_mean), ")"),
       fill = "Mean AUC", x = "Model Type", y = "Modality")
ggsave(file.path(figdir, "heatmap_model_modality.png"), p4, width = 11, height = 7, dpi = 200)

# ---- 9) Trend of AUC over Time (zoomed) ----
p5 <- df %>%
  filter(!is.na(auc_roc), !is.na(year)) %>%
  ggplot(aes(x = year, y = auc_roc)) +
  geom_point(alpha = 0.85) +
  geom_smooth(se = FALSE) +
  coord_cartesian(ylim = auc_zoom) +                      # << zoom
  scale_y_continuous(labels = number_format(accuracy = 0.001)) +
  theme_minimal(base_size = 12) +
  labs(title = "AUC Over Time (Zoomed)", x = "Publication Year", y = "AUC")
ggsave(file.path(figdir, "auc_trend_over_time.png"), p5, width = 10, height = 6, dpi = 200)

# ---- 10) Sensitivity vs Specificity (both zoomed) ----
sens_zoom <- compute_zoom(df$sensitivity, pad = 0.01, min_span = 0.08, hard = c(0.70, 1.00))
spec_zoom <- compute_zoom(df$specificity, pad = 0.01, min_span = 0.08, hard = c(0.70, 1.00))

p6 <- df %>%
  filter(!is.na(sensitivity), !is.na(specificity)) %>%
  ggplot(aes(x = specificity, y = sensitivity, color = model_type)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  coord_cartesian(xlim = spec_zoom, ylim = sens_zoom) +   # << zoom both axes
  theme_minimal(base_size = 12) +
  labs(title = "Sensitivity vs Specificity (Zoomed)", x = "Specificity", y = "Sensitivity",
       color = "Model Type")
ggsave(file.path(figdir, "sens_spec_plot.png"), p6, width = 10, height = 6, dpi = 200)

# ---- 11) Quality plot (only if present) ----
if("quality_rating" %in% names(df) && any(!is.na(df$quality_rating) & df$quality_rating!="")){
  p7 <- df %>%
    filter(!is.na(quality_rating), quality_rating!="") %>%
    ggplot(aes(x = quality_rating, fill = quality_rating)) +
    geom_bar() + theme_minimal(base_size = 12) + theme(legend.position = "none") +
    labs(title = "Quality Assessment of Included Studies", x = "Quality Rating", y = "Count")
  ggsave(file.path(figdir, "quality_barplot.png"), p7, width = 8, height = 5, dpi = 200)
}

# ---- 12) Correlation matrix ----
corr_df <- df %>%
  transmute(AUC = auc_roc, Sensitivity = sensitivity, Specificity = specificity, F1 = f1_precision_num)
png(file.path(figdir, "correlation_matrix.png"), width = 1200, height = 1000, res = 170)
print(GGally::ggpairs(corr_df))
dev.off()


# ---- Console summary ----
if(!is.null(res)){
  cat("\nRandom-effects meta-analysis (logit AUC):\n")
  print(res)
  cat("\nPooled AUC (back-transformed):\n")
  cat(round(transf.ilogit(coef(res)), 4), " (95% CI ",
      round(transf.ilogit(res$ci.lb),4), "–",
      round(transf.ilogit(res$ci.ub),4), ")\n", sep = "")
}
cat("\nSaved figures to: ", normalizePath(figdir), "\n")
