from rdkit import Chem

def generate_canonical_smiles(smiles):
    """
    Converts a SMILES string into its canonical (standardised) format,
    which helps in identifying the molecule regardless of structure variations.
    
    Args:
        smiles (str): Original SMILES string.
        
    Returns:
        str: Canonical SMILES string or "Invalid SMILES" if not valid.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
    return "Invalid SMILES"

def add_canonical_smiles_column(df):
    """
    Adds a new column with canonical SMILES representations to the DataFrame.
    
    Args:
        df (pd.DataFrame): DataFrame with SMILES strings.
        
    Returns:
        pd.DataFrame: DataFrame updated with a 'Canonical_SMILES' column.
    """
    df["Canonical_SMILES"] = df["SMILES"].apply(generate_canonical_smiles)
    df = df.reset_index(drop=True)
    return df
