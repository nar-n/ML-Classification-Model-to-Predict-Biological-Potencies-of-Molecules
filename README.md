# Bioactivity Prediction with Random Forest Classifier

This project focuses on developing machine learning models to predict the bioactivity of molecules sourced from PubChem, classified as either "active" or "inactive" against a specified target. 

## Overview:
- The dataset undergoes preprocessing to ensure quality, involving the conversion of SMILES to Canonical SMILES, and the removal of duplicates and misleading activity labels.
- Molecular descriptors, encompassing physicochemical properties and fingerprints, are computed from the preprocessed dataset.
- The dataset is divided into training and testing sets, and Random Forest classifiers are trained using various combinations of molecular fingerprints and descriptors.
- Performance evaluation is conducted using classification metrics like precision, recall, and F1-score.
- Four models are created and evaluated, differing in the types of descriptors employed, achieving accuracies ranging from 93% to 94%.

# MLModel Activity Project

This project processes molecular SMILES data, standardises it, and removes duplicates to prepare it for machine learning tasks.

## Project Structure

```
MLModel_Activity/
├── data/
│   └── NP_004081-0.0.csv          # Original CSV file with molecular data
├── src/
│   ├── __init__.py                # Initialize src package
│   ├── data_loader.py             # Loads and preprocesses data
│   ├── smiles_processing.py       # SMILES processing and standardisation
│   ├── duplicate_handler.py       # Duplicate detection and handling
│   └── main.py                    # Main script for executing the pipeline
├── requirements.txt               # Project dependencies
└── README.md                      # Project overview and instructions
```

## Instructions

1. **Setup**:
   Install dependencies using:
   ```bash
   pip install -r requirements.txt
   ```

2. **Run**:
   Execute the pipeline by running:
   ```bash
   python src/main.py
   ```

## Description of Modules

- `data_loader.py`: Loads molecular data from a CSV file.
- `smiles_processing.py`: Standardises SMILES strings to canonical format.
- `duplicate_handler.py`: Identifies and displays duplicate SMILES entries.
- `main.py`: Main script to run the entire process.
