def find_duplicates(df, column_name="Canonical_SMILES"):
    """
    Identifies duplicates in the specified column, usually in the canonical SMILES format.
    
    Args:
        df (pd.DataFrame): DataFrame with the canonical SMILES column.
        column_name (str): Column to check for duplicates.
        
    Returns:
        list: List of indices where duplicates are found.
    """
    duplicates_indices = df[df[column_name].duplicated()].index.tolist()
    print(f"Total duplicates found in '{column_name}':", len(duplicates_indices))
    return duplicates_indices

def display_duplicates(df, column_name="Canonical_SMILES"):
    """
    Prints out rows that contain duplicate values in the specified column.
    
    Args:
        df (pd.DataFrame): DataFrame with potential duplicates.
        column_name (str): Column to display duplicates from.
        
    Returns:
        pd.DataFrame: DataFrame containing only duplicate rows.
    """
    duplicates_df = df[df[column_name].duplicated() == True]
    print("Duplicate entries found:\n", duplicates_df)
    return duplicates_df
