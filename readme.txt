GNHSF: Population-based Metaproteomics Reveals Key Microbial Functions in Metabolic Diseases and Aging


├── code
│   ├── figure_R                      # This directory contains the R scripts used for generating the figures in the research paper
│   │   ├── arrange_link.R
│   │   ├── fig1_figs2.R
│   │   ├── fig1.R
│   │   ├── fig2_part_glm.R
│   │   ├── fig2.R
│   │   ├── fig3_part_glmm.R
│   │   ├── fig3.R
│   │   ├── fig4.R
│   │   ├── fig5.R
│   │   ├── fig6_figs8_figs9_figs11.R
│   │   ├── fig6.R
│   │   ├── fig7.R
│   │   ├── figs1_part_getBC.R
│   │   ├── figs1.R
│   │   ├── figs2_mapping.R
│   │   ├── figs3_count.R
│   │   ├── figs4.R
│   │   ├── figs5_figs6.R
│   │   ├── figs7.R
│   │   ├── glmm.R
│   │   └── glm.R
│   └── ML_py                          # This directory contains the python scripts used for machine learning in the research paper
│       ├── evaluate.py                # Calculate AUC and other metrics from predicted probabilities and true labels
│       ├── test_model_extra.py        # Evaluate model performance on the external test set
│       ├── test_model_inter.py        # Evaluate model performance on the internal test set
│       ├── train_model.py             # Train models on different proteomics datasets
│       └── validation_analysis.py
├── figures                            # This directory contains the python scripts and results of ROC-AUC and PR-AUC of machine learning models
│   ├── auc_curves
│   │   ├── external_roc.pdf
│   │   └── internal_roc.pdf
│   ├── boxplot
│   │   ├── seed_auc_boxplot_external.pdf
│   │   ├── seed_auc_boxplot_internal.pdf
│   │   ├── seed_prauc_boxplot_external.pdf
│   │   └── seed_prauc_boxplot_internal.pdf
│   ├── plot_auc_boxplots.py
│   └── plot_pr_auc_boxplots.py
└── readme.txt

