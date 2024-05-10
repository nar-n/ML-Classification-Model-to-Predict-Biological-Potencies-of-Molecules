import pandas as pd
from rdkit import Chem

def load_data(filepath):
    """
    Loads the dataset from a CSV file and converts SMILES strings 
    into molecule objects that can be further analysed.
    
    Args:
        filepath (str): Path to the CSV file.

    Returns:
        pd.DataFrame: DataFrame containing molecules and related data.
    """
    df = pd.read_csv(filepath)
    df['Molecule'] = [Chem.MolFromSmiles(s) for s in df['SMILES']]
    print("Dataset dimensions:", df.shape)
    print("Check if 'Molecule' column has any NA values:", df[['Molecule']].isna().value_counts())
    print("List of SMILES strings that could not be converted:", df[df['Molecule'].isna() == True].SMILES.tolist())
    return df
