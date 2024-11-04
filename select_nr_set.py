#!/usr/bin/env python3
# Select non redundant structures:

# Example of usage:
# python3 select_nr_set.py --general ./structures_annotation/general.txt --pdb_folder pdb_files --output summary_clustering.pdb --nr_folder pdb_nr --distance 6

# This file contains functions to cluster pdb_files according to CDR3a+CDR3b+peptide Levenshtein distance.
#  
# 1) parse_CDR3_extended (anarci_output)
# 2) parse_general_file (general.txt)
# 3) extract_specific_sequences (pdb_file, chain_types_dict)
# 4) calculate_sequence_distance (seq1, seq2)
# 5) get_distance_sum(cdr3a1, cdr3a2, cdr3b1, cdr3b2, peptide1, peptide2)
# 6) cluster_pdbs(distance_matrix, distance_threshold=6)
# 7) vector_to_df(clusters, pdb_ids, col_name1="pdb_id", col_name2="cluster_id")
# 8) get_non_redundant_structures(df_clusters)
# 9) copy_non_redundant_pdbs(pdb_folder, pdb_nonred_folder, df_non_redundant)
# 10) main():

# Load libraries
import os
import shutil
import numpy as np
import pandas as pd
import argparse

from strsimpy.damerau import Damerau
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, cut_tree

from mapping import run_anarci, extract_sequences

def parse_CDR3_extended(anarci_output):
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
                    if 104 <= num <= 118:
                        cdr3 += res
            else:
                print(f"Unexpected line: {parts}")
    
    return cdr3, seq
    
def parse_general_file(general_file):
    """
    Parses the general file and creates a dictionary mapping PDB IDs to chain information.
    
    :param general_file: Path to the general file.
    :return: A dictionary where keys are PDB IDs and values are dictionaries of chain information.
    """
    df = pd.read_csv(general_file, sep='\t')
    pdb_dict = {}
    for _, row in df.iterrows():
        pdb_id = row['pdb.id'].split('.')[0]
        chain_id = row['chain.id']
        chain_type = row['chain.type']
        chain_component = row['chain.component']
        chain_supertype = row['chain.supertype']
        chain_info = {
            'chain.id': chain_id,
            'chain.type': chain_type,
            'chain.component': chain_component,
            'chain.supertype': chain_supertype
        }
        if pdb_id not in pdb_dict:
            pdb_dict[pdb_id] = {}
        pdb_dict[pdb_id][chain_id] = chain_info
    return pdb_dict

def extract_specific_sequences(pdb_file, chain_types_dict):
    """
    Extracts specific sequences (TCRA, TCRB, and peptide) from a PDB file based on chain types,
    excluding any 'X' characters from sequences.

    :param pdb_file: Path to the PDB file.
    :param chain_types_dict: Dictionary mapping PDB IDs to chain information.
    :return: Tuple of sequences (TCRA, TCRB, peptide) without 'X' characters.
    """
    pdb_id = pdb_file.split('/')[-1].split('.')[0]
    seqTCRA, seqTCRB, seqPep = '', '', ''
    
    chain_seq_dict = extract_sequences(pdb_file)
    
    if pdb_id not in chain_types_dict:
        print(f"{pdb_id} not found in the chain types dictionary.")
        return seqTCRA, seqTCRB, seqPep

    skip_x_removal = pdb_id in {"6zky", "6zkz"}
    
    for chain_id, chain_info in chain_types_dict[pdb_id].items():
        chain_type = chain_info['chain.type']
        
        sequence_str = chain_seq_dict.get(chain_id, '')
        if not skip_x_removal:
            sequence_str = sequence_str.replace('X', '')
            
        if chain_type == 'TRA':
            seqTCRA = sequence_str
        elif chain_type == 'TRB':
            seqTCRB = sequence_str
        elif chain_type == 'PEPTIDE':
            seqPep = sequence_str

    return seqTCRA, seqTCRB, seqPep

def calculate_sequence_distance(seq1, seq2):
    """
    Calculates the Damerau-Levenshtein distance between two sequences.
    
    :param seq1: First sequence.
    :param seq2: Second sequence.
    :return: Damerau-Levenshtein distance between seq1 and seq2, as an integer.
    """
    damerau = Damerau()  
    distance = damerau.distance(seq1, seq2)
    return int(distance) if distance.is_integer() else distance

def get_distance_sum(cdr3a1, cdr3a2, cdr3b1, cdr3b2, peptide1, peptide2):
    """
    Computes the sum of distances between CDR3 and peptide sequences.
    
    :param cdr3a1: CDR3 sequence of the first structure.
    :param cdr3a2: CDR3 sequence of the second structure.
    :param cdr3b1: CDR3 sequence of the first structure.
    :param cdr3b2: CDR3 sequence of the second structure.
    :param peptide1: Peptide sequence of the first structure.
    :param peptide2: Peptide sequence of the second structure.
    :return: Sum of distances between the provided sequences.
    """
    dist_cdr3a = calculate_sequence_distance(cdr3a1, cdr3a2)
    dist_cdr3b = calculate_sequence_distance(cdr3b1, cdr3b2)
    dist_peptide = calculate_sequence_distance(peptide1, peptide2)
    dist_sum = dist_cdr3a + dist_cdr3b + dist_peptide
    return dist_sum

def cluster_pdbs(distance_matrix, distance_threshold=6):
    """
    Performs clustering of PDB structures based on a distance matrix.
    
    :param distance_matrix: Numpy array of pairwise distances between structures.
    :param distance_threshold: Distance threshold for clustering.
    :return: Array of cluster IDs for each PDB structure.
    """
    condensed_dist_matrix = squareform(distance_matrix)
    linkage_matrix = linkage(condensed_dist_matrix, method='complete')
    clusters = cut_tree(linkage_matrix, height=distance_threshold).flatten()
    return clusters

def vector_to_df(clusters, pdb_ids, col_name1="pdb_id", col_name2="cluster_id"):
    """
    Converts cluster information and PDB IDs into a DataFrame.
    
    :param clusters: Array of cluster IDs.
    :param pdb_ids: List of PDB IDs.
    :param col_name1: Column name for PDB IDs.
    :param col_name2: Column name for cluster IDs.
    :return: DataFrame with PDB IDs and corresponding cluster IDs.
    """
    df = pd.DataFrame({col_name1: pdb_ids, col_name2: clusters})
    return df

def get_non_redundant_structures(df_clusters):
    """
    Identifies non-redundant structures by first sorting alphabetically by pdb_id,
    then selecting the first pdb_id for each unique cluster.

    :param df_clusters: DataFrame with cluster information.
    :return: DataFrame of non-redundant structures.
    """

    df_sorted = df_clusters.sort_values(by='pdb_id')
    
    df_non_redundant = df_sorted.drop_duplicates(subset='cluster_id', keep='first').reset_index(drop=True)
    
    return df_non_redundant

def copy_non_redundant_pdbs(pdb_folder, pdb_nonred_folder, df_non_redundant):
    """
    Copies non-redundant PDB files to a new directory.
    
    :param pdb_folder: Path to the folder containing the original PDB files.
    :param pdb_nonred_folder: Path to the folder where non-redundant PDB files will be copied.
    :param df_non_redundant: DataFrame with non-redundant PDB IDs.
    """
    if not os.path.exists(pdb_nonred_folder):
        os.makedirs(pdb_nonred_folder)

    for pdb_id in df_non_redundant['pdb_id']:
        src_file = os.path.join(pdb_folder, f"{pdb_id}.pdb")
        dst_file = os.path.join(pdb_nonred_folder, f"{pdb_id}.pdb")
        if os.path.exists(src_file):
            shutil.copy(src_file, dst_file)
        else:
            print(f"File {src_file} does not exist and was not copied.")


def generate_distance_dataframe(sequences):
    """
    Generates a DataFrame with pairs of structures, their Damerau-Levenshtein distances between cdr3a, cdr3b, and peptide,
    and the sum of these distances.
    
    :param sequences: List of dictionaries, each containing pdb_id, cdr3a, cdr3b, and peptide.
    :return: DataFrame with pairs of IDs and corresponding distances.
    """
    data = []
    
    # Calculate distances between all combinations of sequence pairs, including redundancies and themselves
    for i in range(len(sequences)):
        for j in range(len(sequences)):
            # Extract IDs and sequences for each pair
            pdb_id1, pdb_id2 = sequences[i]['pdb_id'], sequences[j]['pdb_id']
            cdr3a1, cdr3a2 = sequences[i]['cdr3a'], sequences[j]['cdr3a']
            cdr3b1, cdr3b2 = sequences[i]['cdr3b'], sequences[j]['cdr3b']
            peptide1, peptide2 = sequences[i]['peptide'], sequences[j]['peptide']
            
            # Calculate Damerau-Levenshtein distances using the defined function
            dist_cdr3a = calculate_sequence_distance(cdr3a1, cdr3a2)
            dist_cdr3b = calculate_sequence_distance(cdr3b1, cdr3b2)
            dist_peptide = calculate_sequence_distance(peptide1, peptide2)
            dist_sum = dist_cdr3a + dist_cdr3b + dist_peptide
            
            # Add pair data to the dataset
            data.append({
                "pdb.id.1": pdb_id1, 
                "pdb.id.2": pdb_id2, 
                "cdr3a.1": cdr3a1, 
                "cdr3b.1": cdr3b1, 
                "peptide.1": peptide1, 
                "cdr3a.2": cdr3a2, 
                "cdr3b.2": cdr3b2, 
                "peptide.2": peptide2, 
                "dist.cdr3a": dist_cdr3a, 
                "dist.cdr3b": dist_cdr3b, 
                "dist.peptide": dist_peptide, 
                "dist.sum": dist_sum
            })
    
    return pd.DataFrame(data)

def main():
    parser = argparse.ArgumentParser(description='Process PDB structures and perform clustering.')
    parser.add_argument("-g","--general", type=str, required=True, help='Path to the general file.')
    parser.add_argument("-p","--pdb_folder", type=str, required=True, help='Path to the folder containing PDB files.')
    parser.add_argument("-o","--output", type=str, default='./structures_annotation/summary_PDB_clustering.csv', help='Path to the output CSV file.')
    parser.add_argument("-nr","--nr_folder", type=str, default='pdb_nr', help='Path to the folder for non-redundant PDB files.')
    parser.add_argument("-d", "--distance", type=float, default=6, help='Distance threshold for clustering.')

    args = parser.parse_args()

    print("Parsing the general file...")
    chain_types_dict = parse_general_file(args.general)
    
    print("Listing PDB files...")
    pdb_files = sorted([os.path.join(args.pdb_folder, f) for f in os.listdir(args.pdb_folder) if f.endswith('.pdb')])
    print("Extracting sequences and calculating distances...")
    sequences = []
    
    for pdb_file in pdb_files:
        pdb_id = os.path.basename(pdb_file).split('.')[0]
        
        print(f"Processing {pdb_id}...")
        seqTCRA, seqTCRB, seqPep = extract_specific_sequences(pdb_file, chain_types_dict)
        cdr3a, cdr3b = '', ''
        if seqTCRA:
            anarci_output_a = run_anarci(seqTCRA, 'TRA')
            if anarci_output_a:
                cdr3a, _ = parse_CDR3_extended(anarci_output_a)
        if seqTCRB:
            anarci_output_b = run_anarci(seqTCRB, 'TRB')
            if anarci_output_b:
                cdr3b, _ = parse_CDR3_extended(anarci_output_b)

        sequences.append({
            'pdb_id': pdb_id,
            'cdr3a': cdr3a,
            'cdr3b': cdr3b,
            'peptide': seqPep
        })

    print("Generating distance DataFrame...")
    df_distances = generate_distance_dataframe(sequences)
    print(df_distances)
    df_distances.to_csv('./distances.csv', index=False)

    print("Calculating distance matrix...")
    pdb_ids = [seq['pdb_id'] for seq in sequences]  # Get pdb_ids
    
    n = len(pdb_ids)  
    distance_matrix = np.zeros((n, n))
    id_to_index = {seq_id: idx for idx, seq_id in enumerate(pdb_ids)}   

    for index, row in df_distances.iterrows():
        id1 = row['pdb.id.1']
        id2 = row['pdb.id.2']
        distance = row['dist.sum']
        
        i = id_to_index[id1]
        j = id_to_index[id2]
        
        distance_matrix[i, j] = distance
        distance_matrix[j, i] = distance  
    distance_matrix = np.round(distance_matrix).astype(int)
    distance_df = pd.DataFrame(distance_matrix, index=pdb_ids, columns=pdb_ids)
    distance_df.to_csv('./distance_matrix.csv', index=True)

    print("Performing clustering...")
    clusters = cluster_pdbs(distance_matrix, distance_threshold=args.distance)

    print("Creating DataFrame with clusters...")
    df_clusters = vector_to_df(clusters, pdb_ids)
    
    print("Identifying non-redundant structures...")
    df_non_redundant = get_non_redundant_structures(df_clusters)
    pdb_nonred = df_non_redundant['pdb_id'].tolist()
    
    df_sequences = pd.DataFrame(sequences)
    df_sequences['cluster_id'] = df_clusters['cluster_id']
    df_sequences['nonred'] = df_sequences['pdb_id'].apply(lambda x: x in pdb_nonred)
    df_sequences = df_sequences[['pdb_id', 'peptide', 'cdr3a', 'cdr3b', 'cluster_id', 'nonred']]
    df_sequences = df_sequences.sort_values(by='cluster_id')
    df_sequences.to_csv(args.output, index=False)
    print(f"Summary DataFrame saved to '{args.output}'.")

    print("Copying non-redundant PDB files...")
    copy_non_redundant_pdbs(args.pdb_folder, args.nr_folder, df_non_redundant)
    print(f"Non-redundant PDB files copied to '{args.nr_folder}'.")

if __name__ == "__main__":
    main()

