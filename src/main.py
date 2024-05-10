from src.data_loader import load_data
from src.smiles_processing import add_canonical_smiles_column
from src.duplicate_handler import find_duplicates, display_duplicates

def main():
    # Load data from CSV file
    filepath = 'data/NP_004081-0.0.csv'
    df = load_data(filepath)

    # Add canonical SMILES column for standardisation
    df = add_canonical_smiles_column(df)

    # Find and display duplicates
    duplicates_indices = find_duplicates(df)
    display_duplicates(df)

if __name__ == "__main__":
    main()
