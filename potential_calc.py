# potential_calc

# This file contains functions to calculate potential from contacts of the training set and given an imput TCR. 
# 1) calculate_potential(contacts_df, peptide=False)
# 2) get_potential(row, df_residues, col_id)

#Import libraries
import pandas as pd
import numpy as np

def calculate_potential(df, peptide=False):
    """
    Calculate potential for all residue pairs based on observed and expected frequencies.
    
    Args:
        df (pd.DataFrame): DataFrame containing residue contact information.
        
    Returns:
        pd.DataFrame: DataFrame with all possible residue pairs and their potential.
    """

    residues = list('ACDEFGHIKLMNPQRSTVWY')  # 20 standard amino acids

    if peptide==True:
        residues_p = list('ADEFGHIKLMNPQRSTVWY') #not C
        all_pairs = pd.DataFrame([(a, b) for a in residues_p for b in residues], columns=['residue_from', 'residue_to'])
    else:
        all_pairs = pd.DataFrame([(a, b) for a in residues for b in residues], columns=['residue_from', 'residue_to'])
        
    # Count the number of observed contacts
    contact_counts = df.groupby(['residue_from', 'residue_to']).size().reset_index(name='count')
    
    # Add pseudocount to each pair
    contact_counts['count'] += 1  # Add pseudocount
    # Merge with all_pairs to ensure every possible residue pair is included
    contact_counts = all_pairs.merge(contact_counts, on=['residue_from', 'residue_to'], how='left')
    # Replace NaN with pseudocount for missing pairs
    contact_counts['count'] = contact_counts['count'].fillna(1)

    # Calculate observed frequencies pobs(a, b)
    total_contacts = contact_counts['count'].sum()
    contact_counts['pobs'] = contact_counts['count'] / total_contacts
    
    # Calculate pobs(a) and pobs(b)
    pobs_a = contact_counts.groupby('residue_from')['pobs'].sum().reset_index(name='pobs_a')
    pobs_b = contact_counts.groupby('residue_to')['pobs'].sum().reset_index(name='pobs_b')
    
    # Merge to get pobs(a) and pobs(b) for each pair
    contact_counts = contact_counts.merge(pobs_a, left_on='residue_from', right_on='residue_from')
    contact_counts = contact_counts.merge(pobs_b, left_on='residue_to', right_on='residue_to')
    
    # Calculate expected frequencies pexp(a, b)
    contact_counts['pexp'] = contact_counts['pobs_a'] * contact_counts['pobs_b']
    
    # Calculate potential
    contact_counts['potential'] = -np.log(contact_counts['pobs'] / contact_counts['pexp'])
    
    return contact_counts[['residue_from', 'residue_to', 'potential']]
    
def get_potential(row, df_residues, col1_id, col2_id):
    """
    Retrieve the potential value for a given residue pair from a DataFrame.
    
    Args:
        row (pd.Series): A row from the DataFrame containing contact information.
        df_residues (pd.DataFrame): DataFrame containing potential values.
        col_id (str): The column name that identifies the residue in the row.
    
    Returns:
        float or None: The potential value if a match is found; otherwise, None.
    """
    try:
        # Get the residue values from the row
        residue_from = row[col1_id]
        residue_to = row[col2_id]
        
        # Find the matching rows in the residues DataFrame
        match = df_residues[(df_residues['residue_from'] == residue_from) & 
                            (df_residues['residue_to'] == residue_to)]
        
        if not match.empty:
            return match['potential'].values[0]
        return None
    
    except KeyError as e:
        print(f"Error: Key {e} not found in row or DataFrame.")
        return None


    
