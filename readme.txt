GNHSF: Population-based Metaproteomics Reveals Key Microbial Functions in Metabolic Diseases and Aging

This repository contains the main analysis codes used to generate figures in the GNHSF study.


├── code
│   ├── figure_R                      # This directory contains the R scripts used for generating the figures in the research paper
│   │   ├── arrange_link.R            # This script establishes the correspondence relationships among peptides, proteins, and taxa, generating the protein-taxa mapping files required by all subsequent analysis scripts.
│   │   ├── fig1_figs2.R              # This script generates: Fig 1D: Average identification counts per sample at each taxonomic level; Fig S2D: Distribution histogram of sample identification counts; Fig S2E: Breakdown of microbial proteins versus human proteins
│   │   ├── fig2.R                    # This script processes and visualizes results from Generalized Linear Model (GLM) analysis. Note: Run fig2_part_glm.R first to obtain complete GLM results. Includes: Fig S7B: GLM associations grouped by clinical categories; Fig 2A: Summary of top 6 associations; Fig 2B-D: Heatmap visualization of the most significant associations at different taxonomic levels
│   │   ├── fig3.R                    # This script analyzes metaproteomic features associated with aging. Note: Run fig3_part_glmm.R first to calculate within-subject associations using Generalized Linear Mixed Models (GLMM). Includes: Fig 3A-F: Aging-associated metaproteomic features; Fig S8A: Functional and taxonomic annotations of age-associated microbial protein groups
│   │   ├── fig4.R                    # This script identifies and visualizes metaproteomic features commonly associated with metabolic diseases. Includes: Fig 4A-C: Shared metaproteomic signatures across metabolic diseases
│   │   ├── fig5.R                    # This script performs medication-weighted GLM calculations and generates related visualizations. Includes: Fig 5B-G: Medication-responsive metaproteomic features in metabolic diseases Fig S11B: Medication-specific proteins and their corresponding species in T2D
│   │   ├── fig6_figs8_figs9_figs11.R # This script analyzes and visualizes T2D-associated features. Note: Perform GLM analysis on the FH cohort and run machine learning code to export ML-related features before executing relevant sections. Includes: Fig 6A: Network visualization of T2D-associated species; Fig 6C: Network visualization of T2D-associated metaproteomic features; Fig S8B: Comparison between metaproteomics and metagenomics data; Fig S9A-B: T2D-related species and their produced microbial protein groups; Fig S11A: GLM associations of M. elsdenii proteins with T2D and T2D medication
│   │   ├── fig7.R                     # This script visualizes in vivo and in vitro biological validation data. Includes: Fig 7: All panels for biological validation experiments
│   │   ├── figs1.R                    # This script generates supplementary figure 1. Note: Run figs1_part_getBC.R first to calculate Bray-Curtis distance matrices for all replicate types. Includes: Fig S1A-B, Fig 1C: Correlation coefficients and Bray-Curtis distances for all replicates; Fig S1C: PCoA of all 2,514 samples
│   │   ├── figs2_mapping.R            # This script calculates the proportion of each sample annotated to taxa or functions for Fig S2A-C: Annotation coverage statistics
│   │   ├── figs3_count.R              # This script generates Fig S3A-H: Count statistics of top features
│   │   ├── figs4.R                    # This script generates all panels in Fig S4: Complete supplementary figure 4
│   │   ├── figs5_figs6.R              # This script calculates and visualizes all core features of the GNHSF metaproteomic dataset. Includes Fig S5 & Fig S6: All panels showing core metaproteomic features
│   │   ├── figs7.R                    # This script performs Fig S7A: PERMANOVA analysis
│   └── ML_py                          # This directory contains the python scripts used for machine learning in the research paper
│       ├── evaluate.py                # Calculate AUC and other metrics from predicted probabilities and true labels
│       ├── test_model_extra.py        # Evaluate model performance on the external test set
│       ├── test_model_inter.py        # Evaluate model performance on the internal test set
│       ├── train_model.py             # Train models on different proteomics datasets
│       └── validation_analysis.py
├── figures                            # This directory contains the python scripts and results of ROC-AUC and PR-AUC of machine learning models
│   ├── auc_curves                     # ROC-AUC curve plots of internal and external tests
│   │   ├── external_roc.pdf
│   │   └── internal_roc.pdf
│   ├── boxplot                        # Box plots of ROC-AUC and PR-AUC
│   │   ├── seed_auc_boxplot_external.pdf
│   │   ├── seed_auc_boxplot_internal.pdf
│   │   ├── seed_prauc_boxplot_external.pdf
│   │   └── seed_prauc_boxplot_internal.pdf
│   ├── plot_auc_boxplots.py          # Generate boxplots illustrating the distribution of ROC-AUC scores across 20 random seeds.
│   └── plot_pr_auc_boxplots.py       # Create boxplots that depict the distribution of PR-AUC (Precision-Recall Area Under the Curve) scores across 20 random seeds
└── readme.txt
