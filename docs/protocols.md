# Protocol: ML/DL models for Alzheimer's diagnosis — systematic review & meta-analysis

## Title
Deep Learning for Early Detection of Alzheimer’s Disease: A Systematic Review and 
Meta-Analysis of Diagnostic Accuracy

## Author
- Saketh Machiraju UCSC Undergraduate 

## Background
Alzheimer’s disease (AD) is the most common cause of dementia, affecting millions of individuals worldwide and placing a significant burden on patients, caregivers, and healthcare systems. Early and accurate diagnosis is essential for effective management, yet current clinical assessments often face challenges in sensitivity, specificity, and scalability. Advances in machine learning (ML) and deep learning (DL) offer the potential to transform AD diagnosis by leveraging neuroimaging, biomarker data, and clinical records to detect disease patterns that are difficult for humans to discern.

In recent years, a rapidly growing body of research has applied ML/DL models—ranging from traditional classifiers such as support vector machines and random forests to convolutional and transformer-based neural networks—for tasks such as distinguishing AD from healthy controls, predicting conversion from mild cognitive impairment, and integrating multimodal data sources. Many of these studies report promising performance metrics, particularly high AUC-ROC values, but results remain inconsistent due to differences in datasets, validation strategies, and reporting practices.

Given this heterogeneity, a systematic review and meta-analysis is needed to rigorously evaluate the diagnostic accuracy of ML/DL models for Alzheimer’s disease. By synthesizing existing evidence across imaging and biomarker modalities, model types, and study designs, this project aims to clarify the current state of performance, assess methodological rigor, and identify gaps that need to be addressed before clinical implementation.

This work builds on prior experience with biomedical imaging, computer vision, and applied machine learning projects, and will emphasize reproducibility through preregistration, open-source code, and transparent reporting


## Objectives
Primary: Summarize and meta-analyze the diagnostic performance of ML/DL models for AD.

## PICO
- Population: Adults with suspected AD (inc. MCI)
- Index: ML/DL diagnostic models (supervised models, CNNs, RF, SVM, Transformers)
- Comparator: Clinical diagnosis or baseline models where reported
- Outcomes: AUC-ROC, sensitivity, specificity, F1, calibration

## Eligibility criteria
- Inclusion: Peer-reviewed/ preprints that report diagnostic performance of ML/DL on human AD datasets; report at least one performance metric with sample sizes.
- Exclusion: Animal studies, purely algorithmic papers without evaluation on human data, review papers.

## Search strategy
- Databases: PubMed, IEEE Xplore, arXiv, Scopus, Web of Science.
- Search terms (example): ("Alzheimer*" OR "MCI") AND ("machine learning" OR "deep learning" OR "CNN" OR "random forest" OR "SVM") AND ("accuracy" OR "AUC" OR "sensitivity" OR "specificity")

## Data extraction
- Study metadata, population, sample sizes, model type, input modality (MRI/PET/CSF), training/validation method, reported metrics and 2×2 confusion table if available.

## Risk of bias
- Use QUADAS-2 or a diagnostic model risk tool; record applicability concerns.

## Data synthesis / ML specifics
- For AUCs: random-effects meta-analysis (metafor)
- For sensitivity/specificity: bivariate/mada models
- Subgroup analyses: modality (MRI vs PET), external validation vs internal only, DL vs classical ML.

## Performance metrics
- Primary: AUC-ROC
- Secondary: sensitivity, specificity, F1, calibration

## Reproducibility & data sharing
- Code: GitHub repository
- Data: list of public datasets used; raw extracted tables included in /data_raw



