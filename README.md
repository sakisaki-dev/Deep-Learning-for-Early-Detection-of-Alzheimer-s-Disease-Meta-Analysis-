# Deep Learning and Machine Learning for Early Detection of Alzheimer's Disease
## A Systematic Review and Meta-Analysis

[![medRxiv](https://img.shields.io/badge/medRxiv-2026.05.21.26353815-blue)](https://medrxiv.org/cgi/content/short/2026.05.21.26353815v1)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

**Preprint:** https://medrxiv.org/cgi/content/short/2026.05.21.26353815v1  
**OSF Protocol:** https://osf.io/h3p6x  
**Author:** Saketh Machiraju · University of California, Santa Cruz

---

## Overview
Systematic review and meta-analysis of 30 studies (2015–2025) evaluating 
machine learning and deep learning diagnostic performance for Alzheimer's 
disease and mild cognitive impairment. Pooled AUC of 0.962 across MRI, 
PET, EEG, and multimodal pipelines using PRISMA 2020 guidelines and 
random-effects meta-analysis in R.

---

## Key Findings
- Pooled AUC: **0.962** (95% CI: 0.939–0.977)
- Mean sensitivity: **0.914** · Mean specificity: **0.913**
- DL+MRI most validated pairing (AUC 0.956, n=11 studies)
- Publication bias confirmed via Egger's test (p=0.0003)

---

## Repo Structure
- `data_clean/` → processed extraction tables
- `analysis/` → R and Python analysis scripts
- `figures/` → forest plots, funnel plots, heatmaps
- `docs/` → protocol and OSF registration

---

## Reproducibility
- R dependencies: `metafor`, `tidyverse`, `GGally`, `viridis`
- Python dependencies: `requirements.txt`
- R setup: `R/install_packages.R`

---

## Citation
If you use this work, please cite:

Machiraju, S. (2026). Deep Learning and Machine Learning for Early 
Detection of Alzheimer's Disease: A Systematic Review and Meta-Analysis. 
*medRxiv*. https://doi.org/10.64898/2026.05.21.26353815

---

## License
CC BY 4.0 — free to share and adapt with attribution
