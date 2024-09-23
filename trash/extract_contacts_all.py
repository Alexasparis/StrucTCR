# Extract contacts:

# This file contains functions to extract contacting residues and filter the resulting dataframe.
# 1) extract_contacts (pdb_file_list, distance)
# 2) filter_contacts (contact_df, minimum_n_contacting_atoms)

#Import libraries
import os
import pandas as pd
import numpy as np
from Bio import PDB


residue_mapping = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
    'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
    'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
    'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

def extract_contacts(pdb_files, distance=5):
    """
    Extract contacting atoms between residues based on a certain distance (default = 5),
    including reciprocal contacts.
    
    Args:
        pdb_files (list of str): List of file paths to PDB files.
        distance (float): Distance threshold for considering contacts.
    
    Returns:
        pd.DataFrame: DataFrame containing contacts with columns ['pdb_id', 'chain_from', 'chain_to', 
                                                                 'residue_from', 'residue_to', 
                                                                 'resid_from', 'resid_to', 
                                                                 'atom_from', 'atom_to', 'dist'].
    """
    contacts = []
    for pdb_file in pdb_files:
        pdb_id = os.path.basename(pdb_file).split('.')[0]
        print(f"Extracting contacts from {pdb_id}")
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure(pdb_id, pdb_file)
        model = structure[0]

        for chain_from in model.get_chains():
            for chain_to in model.get_chains():
                if chain_from != chain_to:
                    residues_from = list(chain_from.get_residues())
                    residues_to = list(chain_to.get_residues())
                    
                    for residue_from in residues_from:
                        for residue_to in residues_to:
                            atoms_from = list(residue_from.get_atoms())
                            atoms_to = list(residue_to.get_atoms())
                            
                            for atom_from in atoms_from:
                                for atom_to in atoms_to:
                                    dist = np.linalg.norm(atom_from.coord - atom_to.coord)
                                    if dist <= distance:  # Threshold for contact
                                        # Map residue names to single-letter codes
                                        res_from_single = residue_mapping.get(residue_from.get_resname(), 'X')  # 'X' for unknown residues
                                        res_to_single = residue_mapping.get(residue_to.get_resname(), 'X')  # 'X' for unknown residues
                                        
                                        # Original contact
                                        contacts.append([
                                            pdb_id, 
                                            chain_from.get_id(), 
                                            chain_to.get_id(), 
                                            res_from_single,  # Single-letter code
                                            res_to_single,  # Single-letter code
                                            residue_from.get_id()[1],  # Residue index
                                            residue_to.get_id()[1],  # Residue index
                                            atom_from.get_id(), 
                                            atom_to.get_id(), 
                                            dist
                                        ])
                                        
                                        # Reciprocal contact
                                        contacts.append([
                                            pdb_id, 
                                            chain_to.get_id(),  # Swap chains
                                            chain_from.get_id(),  # Swap chains
                                            res_to_single,  # Swap residues
                                            res_from_single,  # Swap residues
                                            residue_to.get_id()[1],  # Swap residue indices
                                            residue_from.get_id()[1],  # Swap residue indices
                                            atom_to.get_id(),  # Swap atoms
                                            atom_from.get_id(),  # Swap atoms
                                            dist  # Same distance
                                        ])
    
    return pd.DataFrame(contacts, columns=['pdb_id', 'chain_from', 'chain_to', 'residue_from', 'residue_to', 'resid_from', 'resid_to', 'atom_from', 'atom_to', 'dist'])


def filter_contacts(contacts, tcra_chain, tcrb_chain, peptide_chain, mhc_chain, threshold=2):
    """
    Filter contacting residues to get only those with more than a certain number of atoms contacting (default=2).
    Splits the dataframe to get two dataframes with the contacting residues between TCR-peptide and TCR-MHC respectively.

    Args:
        contacts (df): DataFrame containing contacts between residues. Format ['pdb_id', 'chain_from', 'chain_to', 'residue_from', 'residue_to', 'resid_from', 'resid_to', 'atom_from', 'atom_to', 'dist']
        tcra_chain (str): Chain ID for TCRA chain.
        tcrb_chain (str): Chain ID for TCRB chain.
        peptide_chain (str): Chain ID for peptide chain.
        mhc_chain (str): Chain ID for MHC chain.
        threshold (int): Minimum number of atoms contacting (default = 2).
        
    Returns:
        tuple: DataFrames containing filtered contacts for TCR-peptide and TCR-MHC.    
    """
    
    # Get occurrences with more than the threshold number of atoms contacting
    contacts_unique = contacts[['pdb_id', 'chain_from', 'chain_to', 'residue_from', 'residue_to', 'resid_from', 'resid_to']]
    counts = contacts_unique.groupby(['pdb_id', 'chain_from', 'chain_to', 'residue_from', 'residue_to', 'resid_from', 'resid_to']).size()
    duplicates = counts[counts >= threshold].index
    filtered_contacts = contacts_unique.set_index(['pdb_id', 'chain_from', 'chain_to', 'residue_from', 'residue_to', 'resid_from', 'resid_to'])
    contacts_filtered = filtered_contacts.loc[duplicates].reset_index()
    contacts_filtered_unique = contacts_filtered.drop_duplicates()
    
    # Filter out rows where residue_from or residue_to is 'X'
    contacts_filtered_unique = contacts_filtered_unique[
        (contacts_filtered_unique['residue_from'] != 'X') &
        (contacts_filtered_unique['residue_to'] != 'X')
    ]
    
    # TCR/peptide contacts
    contacts_TCR_p = contacts_filtered_unique[
        (contacts_filtered_unique['chain_from'].isin([tcra_chain, tcrb_chain])) & 
        (contacts_filtered_unique['chain_to'] == peptide_chain)
    ]

    # TCR/MHC contacts
    contacts_TCR_MHC = contacts_filtered_unique[
        (contacts_filtered_unique['chain_from'].isin([tcra_chain, tcrb_chain])) & 
        (contacts_filtered_unique['chain_to'] == mhc_chain)
    ]
    
    return contacts_TCR_p, contacts_TCR_MHC
