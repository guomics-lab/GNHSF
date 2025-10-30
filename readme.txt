.
├── code/                     # Modeling scripts and utilities
│   ├── train_model.py        # Train models on different proteomics datasets
│   ├── test_model_extra.py   # Evaluate model performance on the external test set
│   ├── test_model_inter.py   # Evaluate model performance on the internal test set
│   └── evaluate.py           # Calculate AUC and other metrics from predicted probabilities and true labels
│
├── figures/                  # Scripts for generating visualizations (e.g., boxplots and AUC curves)
│
└── result/                   # Output directory containing modeling results
    └── glm_xgb_gain_freq_42_output_XXX/  # Results for a specific proteomics dataset configuration
        ├── external_preds/               # Prediction outputs on the external test set across different random seeds
        ├── internal_preds/               # Prediction outputs on the internal test set for different cross-validation folds
        ├── models/                       # Trained model checkpoints saved during training
        └── selected_features/            # Features selected by the model during the feature selection process