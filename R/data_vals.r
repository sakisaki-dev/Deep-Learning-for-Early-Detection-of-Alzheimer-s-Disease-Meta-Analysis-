# ---- 1) Load + clean ----
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

# ---- 3) Summary + meta-analysis (safe) ----

library(metafor)
library(dplyr)

# Keep only valid AUCs strictly between 0 and 1
meta_df <- df %>%
  filter(!is.na(auc_roc), auc_roc > 0 & auc_roc < 1, !is.na(sample_size_total)) %>%
  mutate(
    logit_auc = log(auc_roc / (1 - auc_roc)),
    var_auc   = (1 / (auc_roc * (1 - auc_roc))) / sample_size_total
  )

# Check for NA or non-finite variances
meta_df <- meta_df %>% filter(is.finite(logit_auc), is.finite(var_auc))

if (nrow(meta_df) < 3) stop("Not enough valid studies for meta-analysis.")

# --- Random-effects meta-analysis on logit(AUC)
res_auc <- rma(yi = logit_auc, vi = var_auc, data = meta_df, method = "REML")

# Back-transform pooled AUC + CI
pooled_auc <- exp(res_auc$b) / (1 + exp(res_auc$b))
ci_auc_low <- exp(res_auc$ci.lb) / (1 + exp(res_auc$ci.lb))
ci_auc_high <- exp(res_auc$ci.ub) / (1 + exp(res_auc$ci.ub))

# --- Egger’s test
egger_test <- regtest(res_auc, model = "rma")

# --- Descriptive summaries (mean ± sd)
desc_stats <- df %>%
  summarise(
    mean_auc  = mean(auc_roc, na.rm = TRUE),
    sd_auc    = sd(auc_roc, na.rm = TRUE),
    mean_sens = mean(sensitivity, na.rm = TRUE),
    mean_spec = mean(specificity, na.rm = TRUE),
    mean_f1   = mean(f1_precision_num, na.rm = TRUE),
    n_auc     = sum(!is.na(auc_roc)),
    n_sens    = sum(!is.na(sensitivity)),
    n_spec    = sum(!is.na(specificity)),
    n_f1      = sum(!is.na(f1_precision_num))
  )

# ---- 4) Output ----
cat("\n--- Random-Effects Meta-Analysis (logit(AUC)) ---\n")
print(res_auc)

cat("\nPooled AUC (back-transformed): ",
    round(pooled_auc, 3), " [", round(ci_auc_low, 3), ", ", round(ci_auc_high, 3), "]\n")

cat("\n--- Egger’s Test ---\n")
print(egger_test)

cat("\n--- Descriptive Summary ---\n")
print(desc_stats)


