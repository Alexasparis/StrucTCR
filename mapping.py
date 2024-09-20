# Mapping

# This file contains functions to map TCR residues to imgt numbering and to map neoantigen residues to reference epitope. 
# 1) extract_sequences (pdb_file)
# 2) extract_residues_and_resids (pdb_file, chain_id)
# 3) run_anarci (sequence, chain_id)
# 4) parse_anarci_output (anarci_output)
# 5) map_imgt_to_original (imgt_seq_tuple, pdb_resids_tuple)
# 6) list_to_dict(mapping_list)
# 7) get_imgt_mapping_dict (pdb_id, chain_id, mapping_dict)
# 8) map_resid (row, imgt_mapping_dict)
# 9) add_imgt_mappings(df, imgt_mapping_dict)
# 10) map_epitope_residue (row, epitope_sequence)
# 11) global_alignment (seq1, seq2)
# 12) renumber_seq2_based_on_alignment (aligned_seq1, aligned_seq2)
# 13) map_alignment_to_residues(aligned_seq_pdb, aligned_seq_query, residues_tupple)

#Import libraries:
import os
import re
import subprocess
import pandas as pd

from Bio.Align import PairwiseAligner
from Bio import PDB
from extract_contacts import residue_mapping


####### MAP TCR #######

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
    
def extract_residues_and_resids(pdb_file, chain_id):
    """
    Extract the residue IDs and residues (in one-letter code) from a specific chain in a PDB file.
    
    Args:
        pdb_file (str): Path to the PDB file.
        chain_id (str): Chain ID to extract residues from.
    
    Returns:
        list of tuples: List of tuples where each tuple contains (resid, residue_one_letter).
    """
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)

    residues = []
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                for residue in chain:
                    # Extract the residue ID
                    resid = residue.get_id()[1]
                    
                    # Get the 3-letter residue name and convert to 1-letter code
                    resname = residue.get_resname()
                    
                    # Convert 3-letter code to 1-letter code
                    residue_one_letter = PDB.Polypeptide.protein_letters_3to1.get(resname, 'X')  # Use 'X' for unknown residues

                    residues.append((resid, residue_one_letter))
    
    return residues

def run_anarci(sequence, chain_id):
    """
    Run ANARCI with the given sequence and chain ID.

    Args:
        sequence (str): Sequence in one-letter codes.
        chain_id (str): Chain ID for identification.

    Returns:
        str: Output from ANARCI.
    """
    # Define the path for the temporary file
    temp_fasta = 'temp.fasta'
    
    # Write the sequence to the temporary file
    with open(temp_fasta, 'w') as f:
        f.write(f">{chain_id}\n{sequence}\n")
    
    # Ensure the temporary file is removed after processing
    try:
        # Check if file was created successfully
        if not os.path.exists(temp_fasta):
            raise FileNotFoundError(f"{temp_fasta} was not created successfully.")
        
        # Run ANARCI with IMGT scheme and capture both stdout and stderr
        result = subprocess.run(['ANARCI', '-i', temp_fasta, '--scheme', 'imgt'], 
                                capture_output=True, text=True, check=True)
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"ANARCI failed with error:\n{e.stderr}")
        return None
    finally:
        # Remove the temporary file
        if os.path.exists(temp_fasta):
            os.remove(temp_fasta)
        
def parse_anarci_output(anarci_output):
    """
    Parse the output of ANARCI to extract IMGT numbering and residues.
    
    Args:
        anarci_output (str): Output from ANARCI as a string.
    
    Returns:
        list of tuples: A list where each tuple contains (IMGT_number, residue).
    """
    # Regular expression to capture the relevant lines from the ANARCI output
    pattern = r'^([A-Z])\s+(\d+)\s+([A-Z\-])'
    matches = re.findall(pattern, anarci_output, re.MULTILINE)
    
    # Convert matches into a list of tuples (IMGT_number, residue)
    imgt_numbered_seq = []
    for match in matches:
        try:
            chain_letter, imgt_num, residue = match
            imgt_numbered_seq.append((int(imgt_num), residue))
        except ValueError as e:
            print(f"Error processing match: {match}. Error: {e}")
    
    return imgt_numbered_seq

def map_imgt_to_original(imgt_numbered_seq, pdb_resids):
    """
    Map the original numbering of a sequence from the PDB 'resids' to the IMGT numbering.
    
    Args:
        imgt_numbered_seq (list of tuples): The IMGT numbered sequence as tuples (IMGT_number, residue).
        pdb_resids (list of tuples): The original residue numbers from the PDB file as tuples (resid, residue_one_letter).
    
    Returns:
        list of tuples: A list where each tuple contains (original_resid, IMGT_number, residue).
    """
    mapping = []
    pdb_resid_index = 0  # Index to keep track of the PDB resids

    for imgt_pos, residue in imgt_numbered_seq:
        if residue != "-":  # If it's not a gap, map to the original PDB resid
            if pdb_resid_index < len(pdb_resids):
                original_resid, _ = pdb_resids[pdb_resid_index]
                mapping.append((original_resid, imgt_pos, residue))
                pdb_resid_index += 1
            else:
                # This condition should not happen unless there's a mismatch between PDB resids and the sequence length
                mapping.append((None, imgt_pos, residue))
        else:
            # Gaps in IMGT numbering don't correspond to any original PDB resid
            mapping.append((None, imgt_pos, residue))
    
    return mapping

def list_to_dict(mapping_list):
    """
    Converts a list of tuples to a dictionary for easier lookup.
    
    Args:
        mapping_list (list): A list of tuples where each tuple contains (original_residue, mapped_residue, residue_code).
    
    Returns:
        dict: A dictionary mapping original_residue -> mapped_residue.
    """
    return {orig: mapped for orig, mapped, _ in mapping_list if orig is not None}

def get_imgt_mapping_dict(pdb_id, chain, imgt_mappings):
    """
    Retrieves the correct mapping dictionary for a given PDB ID and chain.
    
    Args:
        pdb_id (str): PDB identifier.
        chain (str): Chain used.
        imgt_mappings (dict): Mapping dictionary with PDB ID and chain-specific mappings.
    
    Returns:
        dict: The appropriate mapping dictionary for the given PDB and chain.
    """
    mapping_list = imgt_mappings.get(pdb_id, {}).get(chain, [])
    
    # If it's a list of tuples, convert it to a dictionary
    if isinstance(mapping_list, list):
        return list_to_dict(mapping_list)
    
    return mapping_list  # Already a dictionary

def map_resid(row, imgt_mappings):
    """
    Maps residue IDs to IMGT numbering based on the PDB ID and chain.
    
    Args:
        row (pd.Series): A row from the DataFrame containing residue and chain information.
        imgt_mappings (dict): General IMGT mapping dictionary.
    
    Returns:
        pd.Series: A series with 'imgt_from' mapped values.
    """
    imgt_from = '-'
    
    # Retrieve the IMGT mapping for chain_from 
    imgt_mapping_from = get_imgt_mapping_dict(row['pdb_id'], row['chain_from'], imgt_mappings)
    
    # Apply the mapping
    imgt_from = imgt_mapping_from.get(row['resid_from'], row['resid_from'])
        
    return pd.Series([imgt_from], index=['imgt_from'])

def add_imgt_mappings(df, imgt_mappings):
    """
    Adds IMGT mappings to the DataFrame.

    Args:
        df (pd.DataFrame): Original DataFrame with residue information.
        imgt_mappings (dict): General IMGT mapping dictionary.

    Returns:
        pd.DataFrame: Updated DataFrame with 'imgt_from' columns added.
    """
    # Apply the mapping function to each row of the DataFrame
    mappings = df.apply(lambda row: map_resid(row, imgt_mappings), axis=1)
    
    # Add the 'imgt_from' columns to the original DataFrame
    df[['imgt_from']] = mappings
    
    return df

####### MAP EPITOPE #######

def map_epitope_residue(row, epitope_sequence):
    index = row['resid_to'] - 1
    if 0 <= index < len(epitope_sequence):
        return epitope_sequence[index]
    else:
        return None

####### MAP MHCI #######

def global_alignment(seq1, seq2, match=1, mismatch=-1, gap_open=-10, gap_extend=-0.5):
    """
    Perform global alignment between two sequences using the Needleman-Wunsch algorithm with custom scoring parameters.
    
    Parameters:
    seq1 (str): The first sequence.
    seq2 (str): The second sequence.
    match (int): Score for a match.
    mismatch (int): Penalty for a mismatch.
    gap_open (float): Penalty for opening a gap.
    gap_extend (float): Penalty for extending a gap.
    
    Returns:
    tuple: Contains aligned sequences and the alignment score.
    """
    # Initialize PairwiseAligner with custom scoring parameters
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = match
    aligner.mismatch_score = mismatch
    aligner.open_gap_score = gap_open
    aligner.extend_gap_score = gap_extend
    
    # Perform the alignment
    alignments = aligner.align(seq1, seq2)
    
    # Choose the best alignment (highest score)
    best_alignment = alignments[0]
    
    # Extract aligned sequences and alignment score from the best alignment
    aligned_seq1 = best_alignment[0]
    aligned_seq2 = best_alignment[1]
    score = best_alignment.score
    
    return aligned_seq1, aligned_seq2, score
    
def renumber_seq2_based_on_alignment(seq1, seq2):
    """
    Renumber seq2 based on the alignment with seq1, starting from the first aligned residue.
    
    Args:
        seq1 (str): The first sequence (aligned).
        seq2 (str): The second sequence (aligned).
    
    Returns:
        list of tuples: List of tuples where each tuple contains (original_number, new_number, residue).
    """
    # Find the index of the first aligned residue in seq1
    first_aligned_index_seq1 = next((i for i, char in enumerate(seq1) if char != '-'), None)
    
    if first_aligned_index_seq1 is None:
        raise ValueError("No aligned residue found in seq1.")
    
    # Find all residues in seq2 aligned with seq1
    renumbered_residues = []
    new_number = 1
    
    for i in range(len(seq2)):
        if seq1[i] != '-':  # If seq1 has an aligned residue
            if seq2[i] != '-':  # Only include residues from seq2 that are not gaps
                renumbered_residues.append((i + 1, new_number, seq2[i]))
                new_number += 1
    
    return renumbered_residues
    
def map_alignment_to_residues(aligned_seq_pdb, aligned_seq_query, residues_M):
    """
    Map the aligned residues in the PDB sequence to the aligned query sequence.
    
    Args:
        aligned_seq_pdb (str): The aligned PDB sequence.
        aligned_seq_query (str): The aligned query sequence.
        residues_M (list of tuples): List of (resid, residue) tuples from the PDB.
    
    Returns:
        list of tuples: Mapped residues with their original IDs.
    """
    mapped_residues = []
    pdb_index = 0  # Index for residues_M
    for i in range(len(aligned_seq_pdb)):
        if aligned_seq_pdb[i] != '-':  # If it's not a gap in the PDB sequence
            pdb_resid = residues_M[pdb_index][0]  # Get the original resid
            pdb_res = residues_M[pdb_index][1]  # Get the original residue
            pdb_index += 1  # Move to the next residue in residues_M
        else:
            pdb_resid = '-'  # This represents a gap in the PDB sequence
        
        query_res = aligned_seq_query[i]
        
        # Append the mapping (resid from PDB, residue from PDB, residue from query sequence)
        mapped_residues.append((pdb_resid, pdb_res if aligned_seq_pdb[i] != '-' else '-', query_res))
    
    return mapped_residues

