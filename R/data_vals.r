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




# =========================
# S3 — Sensitivity Analyses
# (drop-in; uses your existing `df`)
# =========================

suppressPackageStartupMessages(library(metafor))
suppressPackageStartupMessages(library(dplyr))

# If you already have figdir defined, great; otherwise set a safe default
if (!exists("figdir")) figdir <- "alz-meta/R/generated_plots"
if (!dir.exists(figdir)) dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
out_csv <- file.path(figdir, "S3_sensitivity_summary.csv")
out_tex <- file.path(figdir, "Table_S3_1_sensitivity.tex")

# ---- Base meta-analysis dataset (matches your pipeline) ----
base_dat <- df %>%
  filter(!is.na(auc_roc), auc_roc > 0 & auc_roc < 1,
         !is.na(sample_size_total), sample_size_total > 3) %>%
  mutate(
    p_clip   = pmin(pmax(auc_roc, 0.001), 0.999),
    yi_logit = qlogis(p_clip),                                      # effect on logit(AUC)
    vi_logit = 1 / (sample_size_total * p_clip * (1 - p_clip)),     # delta-method var (logit)
    vi_raw   = (p_clip * (1 - p_clip)) / sample_size_total          # delta-method var (raw AUC)
  )

stopifnot(nrow(base_dat) >= 3)

# ---- Helper to summarize metafor fits and back-transform ----
summarize_fit <- function(res, scale = c("logit","raw")) {
  scale <- match.arg(scale)
  k    <- res$k
  Q    <- unname(res$QE)
  tau2 <- ifelse(is.null(res$tau2), 0, unname(res$tau2))
  I2   <- if (is.finite(Q) && k > 1 && Q > 0) max(0, (Q - (k - 1)) / Q) * 100 else NA_real_
  
  if (scale == "logit") {
    pr <- try(predict(res, transf = transf.ilogit), silent = TRUE)
    if (inherits(pr, "try-error")) {
      est <- plogis(res$b); lo <- plogis(res$ci.lb); hi <- plogis(res$ci.ub)
      se  <- res$se
      pi_lo <- plogis(as.numeric(res$b) - 1.96 * sqrt(se^2 + tau2))
      pi_hi <- plogis(as.numeric(res$b) + 1.96 * sqrt(se^2 + tau2))
    } else {
      est <- unname(pr$pred); lo <- unname(pr$ci.lb); hi <- unname(pr$ci.ub)
      if (is.null(pr$pi.lb) || is.na(pr$pi.lb)) {
        se <- res$se
        pi_lo <- plogis(as.numeric(res$b) - 1.96 * sqrt(se^2 + tau2))
        pi_hi <- plogis(as.numeric(res$b) + 1.96 * sqrt(se^2 + tau2))
      } else {
        pi_lo <- unname(pr$pi.lb); pi_hi <- unname(pr$pi.ub)
      }
    }
  } else { # raw scale
    pr <- try(predict(res), silent = TRUE)
    est <- if (inherits(pr, "try-error")) as.numeric(res$b) else unname(pr$pred)
    lo  <- if (inherits(pr, "try-error")) as.numeric(res$ci.lb) else unname(pr$ci.lb)
    hi  <- if (inherits(pr, "try-error")) as.numeric(res$ci.ub) else unname(pr$ci.ub)
    se  <- res$se
    pi_lo <- as.numeric(res$b) - 1.96 * sqrt(se^2 + tau2)
    pi_hi <- as.numeric(res$b) + 1.96 * sqrt(se^2 + tau2)
  }
  
  tibble(
    k = k,
    pooled_auc = as.numeric(est),
    ci_lo = as.numeric(lo),
    ci_hi = as.numeric(hi),
    pred_lo = as.numeric(pi_lo),
    pred_hi = as.numeric(pi_hi),
    tau2 = tau2,
    Q = Q,
    I2 = I2
  )
}

# ---- Run the four analyses ----

# 1) Fixed-effects (logit)
res_FE  <- rma(yi = yi_logit, vi = vi_logit, data = base_dat, method = "FE")
sum_FE  <- summarize_fit(res_FE, scale = "logit") %>% mutate(analysis = "Fixed-effects (logit)")

# 2) REML (logit)
res_RE  <- rma(yi = yi_logit, vi = vi_logit, data = base_dat, method = "REML")
sum_RE  <- summarize_fit(res_RE, scale = "logit") %>% mutate(analysis = "REML (logit)")

# 3) Excluding AUC >= 0.99 (REML, logit)
dat_no99 <- base_dat %>% filter(p_clip < 0.99)
if (nrow(dat_no99) >= 3) {
  res_no99 <- rma(yi = yi_logit, vi = vi_logit, data = dat_no99, method = "REML")
  sum_no99 <- summarize_fit(res_no99, scale = "logit") %>% mutate(analysis = "REML (logit), exclude AUC≥0.99")
} else {
  sum_no99 <- tibble(analysis = "REML (logit), exclude AUC≥0.99", k = nrow(dat_no99),
                     pooled_auc = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
                     pred_lo = NA_real_, pred_hi = NA_real_, tau2 = NA_real_, Q = NA_real_, I2 = NA_real_)
}

# 4) External-validation only (REML, logit)
dat_ext <- base_dat %>% filter(external_validation == "Yes")
if (nrow(dat_ext) >= 3) {
  res_ext <- rma(yi = yi_logit, vi = vi_logit, data = dat_ext, method = "REML")
  sum_ext <- summarize_fit(res_ext, scale = "logit") %>% mutate(analysis = "REML (logit), external-only (Yes)")
} else {
  sum_ext <- tibble(analysis = "REML (logit), external-only (Yes)", k = nrow(dat_ext),
                    pooled_auc = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
                    pred_lo = NA_real_, pred_hi = NA_real_, tau2 = NA_real_, Q = NA_real_, I2 = NA_real_)
}

# 5) RAW AUC scale (REML)
res_raw <- rma(yi = p_clip, vi = vi_raw, data = base_dat, method = "REML")
sum_raw <- summarize_fit(res_raw, scale = "raw") %>% mutate(analysis = "REML (raw AUC scale)")

# ---- Collect, save, and print LaTeX-ready lines ----
S3_table <- bind_rows(sum_FE, sum_RE, sum_no99, sum_ext, sum_raw) %>%
  select(analysis, k, pooled_auc, ci_lo, ci_hi, pred_lo, pred_hi, tau2, Q, I2)

write.csv(S3_table, out_csv, row.names = FALSE)
cat("\nSaved S3 sensitivity summary to:", normalizePath(out_csv), "\n\n")

rnd3 <- function(x) ifelse(is.na(x), "NA", sprintf("%.3f", x))

FE_line  <- paste0("\\textbf{Fixed-effects vs REML:} pooled AUC = \\emph{",
                   rnd3(sum_FE$pooled_auc), "} (FE) vs \\emph{", rnd3(sum_RE$pooled_auc), "} (REML).")
no99_line<- paste0("\\textbf{Excluding AUC$\\ge$0.99:} pooled AUC = \\emph{", rnd3(sum_no99$pooled_auc), "}.")
ext_line <- paste0("\\textbf{External-validation only:} pooled AUC = \\emph{", rnd3(sum_ext$pooled_auc), "}.")
raw_line <- paste0("\\textbf{Raw AUC scale:} pooled AUC = \\emph{", rnd3(sum_raw$pooled_auc), "}.")

cat(FE_line, "\n", no99_line, "\n", ext_line, "\n", raw_line, "\n", sep = "")

# ---- 4) Output ----
cat("\n--- Random-Effects Meta-Analysis (logit(AUC)) ---\n")
print(res_auc)

cat("\nPooled AUC (back-transformed): ",
    round(pooled_auc, 3), " [", round(ci_auc_low, 3), ", ", round(ci_auc_high, 3), "]\n")

cat("\n--- Egger’s Test ---\n")
print(egger_test)

cat("\n--- Descriptive Summary ---\n")
print(desc_stats)


