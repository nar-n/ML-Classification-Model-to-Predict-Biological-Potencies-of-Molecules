# Bioactivity Prediction with Machine Learning

This project focuses on developing machine learning models to predict the bioactivity of molecules sourced from PubChem, classified as either "active" or "inactive" against a specified target. 

## Overview:
- The dataset undergoes preprocessing to ensure quality, involving the conversion of SMILES to Canonical SMILES, and the removal of duplicates and misleading activity labels.
- Molecular descriptors, encompassing physicochemical properties and fingerprints, are computed from the preprocessed dataset.
- The dataset is divided into training and testing sets, and Random Forest classifiers are trained using various combinations of molecular fingerprints and descriptors.
- Performance evaluation is conducted using classification metrics like precision, recall, and F1-score.
- Four models are created and evaluated, differing in the types of descriptors employed, achieving accuracies ranging from 93% to 94%.
