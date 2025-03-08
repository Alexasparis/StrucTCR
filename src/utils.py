# utils.py: This file contains utility functions that are shared by the different scripts in the project.

# Finctions:
# 1) residue_mapping: Dictionary mapping three-letter residue names to one-letter codes.
# 2) extract_sequences (pdb_file): Extract sequences for all chains from a PDB file.
# 3) extract_specific_sequences (pdb_file, chain_types_dict): Extracts specific sequences (TCRA, TCRB, and peptide) from a PDB file based on chain types.
# 4) parse_general_file (general_file): Parses the general file and creates a dictionary mapping PDB IDs to specific chain information such as 'tcra_chain', 'tcrb_chain', 'peptide_chain', and 'mhc_chain'.
# 5) calculate_sequence_distance (seq1, seq2): calculates lenshtein ditance between sequences.

import os
import pandas as pd
from Bio import PDB
from Levenshtein import distance as levenshtein_distance

residue_mapping = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
    'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
    'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
    'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}

def extract_sequences(pdb_file):
    """
    Extract sequences for all chains from a PDB file.

    Args:
        pdb_file (str): Path to the PDB file.

    Returns:
        dict: Dictionary with chain IDs as keys and sequences as values.
    """
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)
    sequences = {}

    for model in structure:
        for chain in model.get_chains():
            chain_id = chain.get_id()
            sequence = []
            for residue in chain:
                if PDB.is_aa(residue):  # Check if the residue is an amino acid
                    res_name = residue.get_resname()
                    sequence.append(residue_mapping.get(res_name, 'X'))  # 'X' for unknown residues
            # Join one-letter codes to form the sequence
            sequences[chain_id] = ''.join(sequence)
    
    return sequences

def extract_specific_sequences(pdb_file, chain_types_dict):
    """
    Extracts specific sequences (TCRA, TCRB, and peptide) from a PDB file based on chain types.
    
    :param pdb_file: Path to the PDB file.
    :param chain_types_dict: Dictionary mapping PDB IDs to chain information.
    :param extract_sequences_fn: Function to extract sequences from the PDB file.
    :return: Tuple of sequences (TCRA, TCRB, peptide).
    """
    pdb_id = os.path.basename(pdb_file).split('.')[0]  # Extract PDB ID from the file path
    seqTCRA, seqTCRB, seqPep = '', '', ''
    
    # Call external function to extract sequences from the PDB file
    chain_seq_dict = extract_sequences(pdb_file)
    
    # Check if the pdb_id is in the chain_types_dict
    if pdb_id not in chain_types_dict:
        print(f"{pdb_id} not found in the chain types dictionary.")
        return seqTCRA, seqTCRB, seqPep

    # Get chain information for the current pdb_id
    chains = chain_types_dict.get(pdb_id, {'tcra_chain': 'D',
                                            'tcrb_chain': 'E',
                                            'peptide_chain': 'C',
                                            'b2m_chain': 'B',
                                            'mhc_chain': 'A'})
    
    # Extract sequences based on the chain information
    if chains['tcra_chain'] is not None:
        seqTCRA = chain_seq_dict.get(chains['tcra_chain'], '')
    if chains['tcrb_chain'] is not None:
        seqTCRB = chain_seq_dict.get(chains['tcrb_chain'], '')
    if chains['peptide_chain'] is not None:
        seqPep = chain_seq_dict.get(chains['peptide_chain'], '')
    
    return seqTCRA, seqTCRB, seqPep

def parse_general_file(general_file):
    """
    Parses the general file and creates a dictionary mapping PDB IDs to specific chain information
    such as 'tcra_chain', 'tcrb_chain', 'peptide_chain', and 'mhc_chain'.
    
    :param general_file: Path to the general file.
    :return: A dictionary where keys are PDB IDs and values are dictionaries with chain information
    """
    # Read the general file into a pandas DataFrame
    df = pd.read_csv(general_file, sep='\t')
    pdb_dict = {}

    # Group by PDB ID and iterate through each group
    for pdb_id, group in df.groupby('pdb.id'):
        chains = {
            'tcra_chain': None,
            'tcrb_chain': None,
            'peptide_chain': None,
            'mhc_chain': None,
            'b2m_chain': None}

        # Iterate through the rows in the group to assign chain types
        for _, row in group.iterrows():
            chain_id = row['chain.id']
            chain_type = row['chain.type']
            chain_component = row['chain.component']
            chain_supertype = row['chain.supertype']

            # Assign chain IDs based on their component and type
            if chain_component == 'TCR' and chain_type == 'TRA':
                chains['tcra_chain'] = chain_id
            elif chain_component == 'TCR' and chain_type == 'TRB':
                chains['tcrb_chain'] = chain_id
            elif chain_component == 'PEPTIDE':
                chains['peptide_chain'] = chain_id
            elif chain_component == 'MHC' and chain_supertype == 'MHCI' and chain_type == 'MHCa':
                chains['mhc_chain'] = chain_id
            elif chain_component == 'MHC' and chain_supertype == 'MHCI' and chain_type == 'MHCb':
                chains['b2m_chain'] = chain_id
        
        # Add the chain information for this pdb_id to the dictionary
        pdb_dict[pdb_id] = chains

    return pdb_dict

def calculate_sequence_distance(seq1, seq2):
    """
    Calculates the Levenshtein distance between two sequences.
    
    :param seq1: First sequence.
    :param seq2: Second sequence.
    :return: Levenshtein distance between seq1 and seq2.
    """
    return levenshtein_distance(seq1, seq2)
