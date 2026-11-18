# Deep Learning and Machine Learning for Early Detection of Alzheimer's Disease
## A Systematic Review and Meta-Analysis

[![medRxiv](https://img.shields.io/badge/medRxiv-2026.05.21.26353815-blue)](https://medrxiv.org/cgi/content/short/2026.05.21.26353815v1)

**Preprint:** https://medrxiv.org/cgi/content/short/2026.05.21.26353815v1
**OSF registration:** https://osf.io/h3p6x
**Author:** Saketh Machiraju, University of California, Santa Cruz

## Overview

A PRISMA 2020 systematic review and random-effects meta-analysis of how well machine learning and deep learning models diagnose Alzheimer's disease and mild cognitive impairment. Studies were searched in PubMed, IEEE Xplore and arXiv, screened in Rayyan, and the extracted AUC, sensitivity and specificity values were analysed in R with `metafor`.

## Study selection

From `R/generated_plots/prisma_flowdiagram.png` and `docs/search_queries.txt`:

* 455 records identified (242 PubMed, 212 IEEE Xplore, 1 arXiv)
* 346 records screened, 90 full-text reports assessed
* 60 reports excluded (45 with no AUC-ROC, 12 repeated method and dataset, 3 out of scope)
* 30 studies included

## Method

`R/meta_analysis.r` and `R/data_vals.r`:

* Each study's AUC is logit-transformed, with variance approximated as `1 / (n * AUC * (1 - AUC))`.
* Pooled with a random-effects model (`metafor::rma`, REML) and back-transformed to the AUC scale.
* Egger's regression test (`regtest`) for funnel-plot asymmetry.
* Sensitivity analyses: fixed-effects, excluding AUC >= 0.99, external-validation studies only, and pooling on the raw AUC scale.
* Figures: forest plots, funnel plot, AUC by model type and modality, model-by-modality heatmap, sample size vs AUC, AUC over time, sensitivity vs specificity, and a correlation matrix.

## Results

See the preprint for the full results. The committed outputs in this repo show:

* Random-effects pooled AUC of 0.963 (95% CI 0.940 to 0.978) across the 22 studies that reported an AUC strictly between 0 and 1 and a sample size (sensitivity summary, commit `4818cb6`).
* Very high heterogeneity: I² = 99.2%, with a 95% prediction interval of 0.692 to 0.997.
* Pooled AUC of 0.947 after excluding AUC >= 0.99 (k = 19) and 0.977 for externally validated studies only (k = 6).
* Deep learning on MRI is the most common pairing: mean AUC 0.956 across 11 studies (`heatmap_model_modality.png`).

## Repository layout

| Path | Contents |
| :--- | :--- |
| `R/meta_analysis.r` | Data cleaning, meta-analysis and all figures |
| `R/data_vals.r` | Pooled estimate, Egger's test, descriptive stats and the sensitivity-analysis table. Run after `meta_analysis.r` in the same R session, since it reuses `in_csv` and loaded packages |
| `R/install_packages.R` | Installs `metafor`, `mada`, `readr`, `dplyr`, `tidyr`, `ggplot2` and `dmetar` |
| `R/generated_plots/` | Committed figures, including the PRISMA flow diagram |
| `analysis/articles-selected.csv` | Rayyan export of the 90 full-text reports assessed |
| `docs/protocols.md` | Review protocol (PICO, eligibility, synthesis plan) |
| `docs/search_queries.txt` | Database search strings, filters and hit counts |

## Reproducing

The extracted study-level table (`data_clean/data_extracted_clean.csv`) is not in the repository (`data_clean/` is gitignored), so the analysis cannot be rerun from this repo alone. With that file in place:

```bash
# meta_analysis.r reads and writes paths under alz-meta/, so run R from the
# directory that contains the repo checked out as alz-meta/
Rscript alz-meta/R/install_packages.R
R -e 'source("alz-meta/R/meta_analysis.r"); source("alz-meta/R/data_vals.r")'
```

`meta_analysis.r` also installs `tidyverse`, `janitor`, `scales`, `viridis` and `GGally` if missing, and deletes and regenerates `R/generated_plots/` on each run. `environment.yml` describes a Python environment, but no Python code is in the repo.

## Limitations

* Heterogeneity is extreme (I² above 99%), so the pooled AUC summarises a very mixed set of datasets, validation designs and tasks.
* Only 6 of the pooled studies reported external validation.
* The AUC variance is a sample-size approximation, not the variance reported by each study.

## Citation

Machiraju, S. (2026). Deep Learning and Machine Learning for Early Detection of Alzheimer's Disease: A Systematic Review and Meta-Analysis. *medRxiv*. https://doi.org/10.64898/2026.05.21.26353815

## License

MIT (see `LICENSE`).
