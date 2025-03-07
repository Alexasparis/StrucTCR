# Extract contacts: This file contais functions to extract and filter contacts

# This file contains functions to extract contacting residues and filter the resulting dataframe.
# You can apply this functions with the script ./script/contact_maps_pdb.py and with the notebook ./notebooks/extract_contacts.ipynb.
# 1) extract_contacts (pdb_file_list, chain_dict, distance=5): extract contacts between TCR-P and TCR-MHC
# 2) filter_contacts(contacts_df, tcra_chain_id, tcrb_chain_id, peptide_chain_id, mhc_chain_id, threshold_n_atoms, remove_X=True): Splits df into TCR-p and TCR-MHC

#Import libraries
import os
import pandas as pd
import numpy as np
from Bio import PDB
from utils import residue_mapping

def extract_contacts(pdb_files, chain_dict, distance=5):
    """
    Extract contacts between residues of TCR chains and peptide/TCR chains and MHC in a list of PDB files.

    Args:
        pdb_files (list): List of PDB files. Can be a list with a single PDB file.
        chain_dict (dict): Dictionary mapping PDB IDs to chain information. Obtained with parse_general_file function from utils.py
        distance (float): Distance threshold for contacting residues (default=5).

    Returns:
        df: DataFrame containing contacts between residues. Format ['pdb_id', 'chain_from', 'chain_to', 'residue_from', 'residue_to', 'resid_from', 'resid_to', 'atom_from', 'atom_to', 'dist']
    """
    contacts = []
    
    if residue_mapping is None:
        print("Error: residue_mapping is None")
        return pd.DataFrame() 
    
    for pdb_file in pdb_files:
        try:
            pdb_id = os.path.basename(pdb_file).split('.')[0]
            print(f"Extracting contacts from {pdb_id}")
            parser = PDB.PDBParser(QUIET=True)
            structure = parser.get_structure(pdb_id, pdb_file)
            model = structure[0]

            chains = chain_dict.get(pdb_id)
            if not chains:
                continue
            
            chain_pairs = [
                (chains['tcra_chain'], chains['mhc_chain']),
                (chains['tcrb_chain'], chains['mhc_chain']),
                (chains['tcra_chain'], chains['peptide_chain']),
                (chains['tcrb_chain'], chains['peptide_chain'])]
            
            for chain_from_id, chain_to_id in chain_pairs:
                if chain_from_id and chain_to_id: 
                    try:
                        chain_from = model[chain_from_id]
                        chain_to = model[chain_to_id]
                    except KeyError:
                        print(f"Chain not found in {pdb_id}: {chain_from_id} or {chain_to_id}")
                        continue
                    
                    residues_from = list(chain_from.get_residues())
                    residues_to = list(chain_to.get_residues())
                    
                    for residue_from in residues_from:
                        for residue_to in residues_to:
                            if residue_from.id[0] == " " and residue_to.id[0] == " ":
                                atoms_from = list(residue_from.get_atoms())
                                atoms_to = list(residue_to.get_atoms())
                                
                                for atom_from in atoms_from:
                                    if atom_from.get_parent().id[0] == " " and atom_from.get_name() != "HETATM":
                                        for atom_to in atoms_to:
                                            if atom_to.get_parent().id[0] == " " and atom_to.get_name() != "HETATM":
                                                dist = np.linalg.norm(atom_from.coord - atom_to.coord)
                                                if dist <= distance:
                                                    res_from_single = residue_mapping.get(residue_from.get_resname(), residue_from.get_resname())
                                                    res_to_single = residue_mapping.get(residue_to.get_resname(), residue_to.get_resname())
                                                    
                                                    contacts.append([
                                                        pdb_id, 
                                                        chain_from.get_id(), 
                                                        chain_to.get_id(), 
                                                        res_from_single, 
                                                        res_to_single, 
                                                        residue_from.get_id()[1], 
                                                        residue_to.get_id()[1], 
                                                        atom_from.get_id(), 
                                                        atom_to.get_id(), 
                                                        dist])
        except Exception as e:
            print(f"Error extracting contacts in {pdb_id}: {e}")
            continue

    return pd.DataFrame(contacts, columns=['pdb_id', 'chain_from', 'chain_to', 'residue_from', 'residue_to', 'resid_from', 'resid_to', 'atom_from', 'atom_to', 'dist'])



def filter_contacts(contacts, tcra_chain, tcrb_chain, peptide_chain, mhc_chain, threshold=1, remove_X=True):
    """
    Filter contacting residues to get only those with more than a certain number of atoms contacting (default=2).
    Splits the dataframe to get two dataframes with the contacting residues between TCR-peptide and TCR-MHC respectively.
    
    Args:
        contacts (df): DataFrame containing contacts between residues. Format ['pdb_id', 'chain_from', 'chain_to', 'residue_from', 'residue_to', 'resid_from', 'resid_to', 'atom_from', 'atom_to', 'dist']
        tcra_chain (str): Chain ID for TCRA chain.
        tcrb_chain (str): Chain ID for TCRB chain.
        peptide_chain (str): Chain ID for peptide chain.
        mhc_chain (str): Chain ID for MHC chain.
        threshold (int): Minimum number of atoms contacting (default = 1).
        remove_X (bool): Whether to filter out residues with 'X'. Default is True.
        
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
    
    # Remove rows with 'X' if remove_X is True
    if remove_X:
        contacts_filtered_unique = contacts_filtered_unique[
            (contacts_filtered_unique['residue_from'] != 'X') &
            (contacts_filtered_unique['residue_to'] != 'X')]
    
    # TCR/peptide contacts
    contacts_TCR_p = contacts_filtered_unique[
        (contacts_filtered_unique['chain_from'].isin([tcra_chain, tcrb_chain])) & 
        (contacts_filtered_unique['chain_to'] == peptide_chain)]

    # TCR/MHC contacts
    contacts_TCR_MHC = contacts_filtered_unique[
        (contacts_filtered_unique['chain_from'].isin([tcra_chain, tcrb_chain])) & 
        (contacts_filtered_unique['chain_to'] == mhc_chain)]
    
    return contacts_TCR_p, contacts_TCR_MHC

def filter_contacts_weighted(contacts, tcra_chain, tcrb_chain, peptide_chain, mhc_chain, remove_X=True):
    """
    Filter the contacts between residues for TCR-peptide and TCR-MHC interactions without grouping.
    
    Args:
        contacts (df): DataFrame containing contacts between residues. Format ['pdb_id', 'chain_from', 'chain_to', 'residue_from', 'residue_to', 'resid_from', 'resid_to', 'atom_from', 'atom_to', 'dist']
        tcra_chain (str): Chain ID for TCRA chain.
        tcrb_chain (str): Chain ID for TCRB chain.
        peptide_chain (str): Chain ID for peptide chain.
        mhc_chain (str): Chain ID for MHC chain.
        remove_X (bool): Whether to filter out residues with 'X'. Default is True.
        
    Returns:
        tuple: DataFrames containing raw contact occurrences for TCR-peptide and TCR-MHC.
    """
    # Remove rows with 'X' if remove_X is True
    if remove_X:
        contacts = contacts[
            (contacts['residue_from'] != 'X') &
            (contacts['residue_to'] != 'X')]
    
    # TCR/peptide contacts
    contacts_TCR_p = contacts[
        (contacts['chain_from'].isin([tcra_chain, tcrb_chain])) & 
        (contacts['chain_to'] == peptide_chain)]

    # TCR/MHC contacts
    contacts_TCR_MHC = contacts[
        (contacts['chain_from'].isin([tcra_chain, tcrb_chain])) & 
        (contacts['chain_to'] == mhc_chain)]
    
    return contacts_TCR_p, contacts_TCR_MHC

