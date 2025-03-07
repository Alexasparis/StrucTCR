# mapping.py: This file contains functions to map TCR residues to imgt numbering and to map neoantigen residues to reference epitope. 

# Functions:
# 1) extract_residues_and_resids (pdb_file, chain_id)
# 2) run_anarci (sequence)
# 3) parse_anarci_output (anarci_output)
# 4) parse_CDR3(anarci_output)
# 5) map_imgt_to_original (imgt_seq_tuple, pdb_resids_tuple)
# 6) get_imgt_mapping_dict (pdb_id, chain_id, mapping_dict)
# 7) map_resid (row, imgt_mapping_dict)
# 8) add_imgt_mappings(df, imgt_mapping_dict)
# 9) map_epitope_residue (row, epitope_sequence)
# 10) global_alignment (seq1, seq2)
# 11) renumber_seq2_based_on_alignment (aligned_seq1, aligned_seq2)
# 12) map_alignment_to_residues(aligned_seq_pdb, aligned_seq_query, residues_tupple)

#Import libraries:
import subprocess
import pandas as pd
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

def run_anarci(sequence):
    try:
        if '*' in sequence:
            print("Warning: Stop codon (*) found in the sequence. Removing it...")
            sequence = sequence.replace('*', '')
        elif 'X' in sequence:
            print("Warning: Unknown amino acid (X) found in the sequence. Removing it...")
            sequence = sequence.replace('X', '')
            
        command = f"ANARCI -i {sequence} --scheme imgt"
        result = subprocess.run(command, shell=True, capture_output=True, text=True, check=True)
        
        return result.stdout
    
    except subprocess.CalledProcessError as e:
        print(f"Error: {e.stderr}")
        return e.stdout

def parse_anarci_output(anarci_output):
    """
    Parse the output of ANARCI to extract IMGT numbering and residues, ensuring uniqueness of IMGT numbers.
    
    Args:
        anarci_output (str): Output from ANARCI as a string.
    
    Returns:
        list of tuples: A list where each tuple contains (IMGT_number, residue).
    """
    lines = anarci_output.split('\n')
    imgt_numbered_seq = []
    
    for line in lines:
        if line.startswith('#') or not line.strip():
            continue
        parts = line.split()
        if len(parts) >= 3:
            chain_letter = parts[0]
            imgt_num = int(parts[1])  # IMGT number
            residue = parts[-1]        # Residue
            imgt_numbered_seq.append((imgt_num, residue))
    return imgt_numbered_seq

def parse_CDR3(parsed_anarci_output):
    """
    Parse the CDR1, CDR2, and CDR3 regions from the IMGT-numbered sequence output.
    
    Args:
        parsed_anarci_output (list of tuples): The IMGT-numbered sequence as a list of tuples.
        
    Returns:
        seq_cdr3, seq_var_full: The CDR3 sequence and the full variable region sequence.
    """
    cdr3_seq = ''.join([residue for imgt_num, residue in parsed_anarci_output if 104 <= imgt_num <= 118 and residue != "-"])
    full_var_seq = ''.join([residue for imgt_num, residue in parsed_anarci_output if residue != "-"])
    return cdr3_seq, full_var_seq

def map_imgt_to_original(imgt_numbered_seq, pdb_resids):
    """
    Map the original numbering of a sequence from the PDB 'resids' to the IMGT numbering.
    
    Args:
        imgt_numbered_seq (list of tuples): The IMGT numbered sequence as tuples (IMGT_number, residue).
        pdb_resids (list of tuples): The original residue numbers from the PDB file as tuples (resid, residue_one_letter).
    
    Returns:
        list of tuples: A list where each tuple contains (original_resid, IMGT_number, residue).
    """
    full_var_seq = ''.join([residue for imgt_num, residue in imgt_numbered_seq if residue != "-"])
    full_pdb_seq = ''.join([residue for resid, residue in pdb_resids if residue != "-"])
    var_start = full_pdb_seq.find(full_var_seq)
    
    if var_start == -1:
        raise ValueError("Full variable region not found in the protein sequence.")
    
    var_end = var_start + len(full_var_seq) - 1
    try:
        tupples_of_indexes_var_start_to_var_end = []
        for i in range(var_start,var_end+1):
            tupples_of_indexes_var_start_to_var_end.append(pdb_resids[i])
        
        tupples_imgt_valid = [t for t in imgt_numbered_seq if t[1] != '-']
        
        if len(tupples_of_indexes_var_start_to_var_end) != len(tupples_imgt_valid):
            raise ValueError("The variable region and IMGT sequence lengths do not match. Alignment may be incorrect.")
        
        mapping = [(tupples_of_indexes_var_start_to_var_end[i][0], tupples_imgt_valid[i][0], tupples_imgt_valid[i][1])
                   for i in range(len(tupples_imgt_valid))]

    except Exception as e:
        raise ValueError(f"An error occurred while processing the sequences: {str(e)}")

    return mapping

def get_imgt_mapping_dict(pdb_id, chain, imgt_mappings):
    """
    Retrieves the correct mapping dictionary for a given PDB ID and chain.
    
    Args:
        pdb_id (str): PDB identifier.
        chain (str): Chain used.
        imgt_mappings (dict): Mapping dictionary with PDB ID and chain-specific mappings.
    
    Returns:
        dict: The appropriate mapping dictionary for the given PDB and chain {1:2, 2:3, 3:6 ...}.
    """
    mapping_list = imgt_mappings.get(pdb_id, {}).get(chain, [])
    
    if isinstance(mapping_list, list):
        return {orig: mapped for orig, mapped, _ in mapping_list if orig is not None}
    
    return mapping_list  

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
    imgt_mapping_from = get_imgt_mapping_dict(row['pdb_id'], row['chain_from'], imgt_mappings)
    row_resid_from = row['resid_from']
    imgt_from = imgt_mapping_from.get(row_resid_from, row_resid_from)
        
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
    mappings = df.apply(lambda row: map_resid(row, imgt_mappings), axis=1)

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
        list of tuples: Mapped residues with their original IDs. (pdb_resid, pdb_res, query_res)
    """
    mapped_residues = []
    pdb_index = 0 
    for i in range(len(aligned_seq_pdb)):
        if aligned_seq_pdb[i] != '-':  
            pdb_resid = residues_M[pdb_index][0] 
            pdb_res = residues_M[pdb_index][1]  
            pdb_index += 1  
        else:
            pdb_resid = '-'  
             
        query_res = aligned_seq_query[i]
        mapped_residues.append((pdb_resid, pdb_res if aligned_seq_pdb[i] != '-' else '-', query_res))

    return mapped_residues
