# Mapping

# This file contains functions to map TCR residues to imgt numbering and to map neoantigen residues to reference epitope. 
# 1) extract_residues_and_resids (pdb_file, chain_id)
# 2) run_anarci (sequence, chain_id)
# 3) parse_anarci_output (anarci_output)
# 4) parse_CDR3(anarci_output)
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
import re
import subprocess
import pandas as pd
import concurrent
from Bio.Align import PairwiseAligner
from Bio import PDB


####### MAP TCR #######
    
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
    Execute ANARCI to assign IMGT numbering to a TCR sequence.
    """
    try:
        command=f"ANARCI -i {sequence} --scheme imgt"
        result = subprocess.run(command, shell=True, capture_output=True, text=True, check=True)
        return result.stdout
    except subprocess.CalledProcessError as e:
        return f"Error with {chain_id}: {e.stderr}"

def run_anarci_parallel(sequences, chain_ids):
    results = {}
    with concurrent.futures.ThreadPoolExecutor() as executor:
        future_to_chain = {executor.submit(run_anarci, seq, chain): chain for seq, chain in zip(sequences, chain_ids)}
        for future in concurrent.futures.as_completed(future_to_chain):
            chain_id = future_to_chain[future]
            try:
                result = future.result()  
                results[chain_id] = result
            except Exception as exc:
                print(f'Error processing {chain_id}: {exc}')
    return results
        
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

def parse_CDR3(anarci_output):
    cdr3, seq = "", ""
    for i in anarci_output.split("\n"):
        if i and i[0] != "#" and i[0] != "/":
            parts = i.rstrip().split()
            if len(parts) >= 3:
                _, num, *residues = parts
                num = int(num)  
                res = residues[-1]  
                if res != "-":
                    seq += res
                    if 104 <= num <= 118: #105-117
                        cdr3 += res
            else:
                print(f"Línea inesperada: {parts}")
    
    return cdr3, seq

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
    pdb_resid_index = 0  # Index for PDB residues
    for imgt_pos, residue in imgt_numbered_seq:
        if residue != "-":  # Only process non-gap residues in IMGT
            for original_resid, residue1 in pdb_resids[pdb_resid_index:]:
                if residue1 == residue:
                    mapping.append((original_resid, imgt_pos, residue))
                    pdb_resid_index += 1
                    break
                else:
                    pdb_resid_index += 1
            else:
                mapping.append((None, imgt_pos, residue))
        else:
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
