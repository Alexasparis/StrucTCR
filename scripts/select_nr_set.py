#!/usr/bin/env python3
# Select non redundant structures based on a Levenshtein distance threshold for CDR3a, CDR3b, and peptide sequences:

# Example of usage:
# python3 select_nr_set.py --pdb_folder pdb_files --output summary_clustering.pdb --nr_folder pdb_nr --distance 6 --method complete

# This file contains functions to cluster pdb_files according to CDR3a+CDR3b+peptide Levenshtein distance.

# 1) calculate_sequence_distance (seq1, seq2)
# 2) get_distance_sum(cdr3a1, cdr3a2, cdr3b1, cdr3b2, peptide1, peptide2)
# 3) cluster_pdbs(distance_matrix, distance_threshold=6)
# 4) vector_to_df(clusters, pdb_ids, col_name1="pdb_id", col_name2="cluster_id")
# 5) get_non_redundant_structures(df_clusters)
# 6) copy_non_redundant_pdbs(pdb_folder, pdb_nonred_folder, df_non_redundant)
# 7) main():

# Load libraries
import os
import shutil
import numpy as np
import pandas as pd
import argparse

from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, fcluster
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))
from utils import extract_specific_sequences, parse_general_file, extract_sequences, calculate_sequence_distance
from mapping import run_anarci, parse_CDR3
from config import STRUCTURES_ANNOTATION_DIR

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

def cluster_pdbs(distance_matrix, distance_threshold=6, method='complete'):
    """
    Performs clustering of PDB structures based on a distance matrix.
    
    :param distance_matrix: Numpy array of pairwise distances between structures.
    :param distance_threshold: Distance threshold for clustering.
    :return: Array of cluster IDs for each PDB structure.
    """
    condensed_dist_matrix = squareform(distance_matrix)
    linkage_matrix = linkage(condensed_dist_matrix, method=method)
    clusters = fcluster(linkage_matrix, distance_threshold, criterion='distance')
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
    Identifies non-redundant structures by taking the first structure of each cluster, sorted alphabetically.
    
    :param df_clusters: DataFrame with cluster information.
    :return: DataFrame of non-redundant structures.
    """
    df_clusters_sorted = df_clusters.sort_values(by=['cluster_id', 'pdb_id'])
    df_non_redundant = df_clusters_sorted.groupby('cluster_id').first().reset_index()
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

def main():
    parser = argparse.ArgumentParser(description='Process PDB structures and perform clustering.')
    parser.add_argument("-p","--pdb_folder", type=str, required=True, help='Path to the folder containing PDB files.')
    parser.add_argument("-o","--output", type=str, default='./structures_annotation/summary_PDB_clustering.csv', help='Path to the output CSV file.')
    parser.add_argument("-nr","--nr_folder", type=str, default='pdb_nr', help='Path to the folder for non-redundant PDB files.')
    parser.add_argument("-d", "--distance", type=float, default=6, help='Distance threshold for clustering.')
    parser.add_argument("-m", "--method", type=str, default='complete', help='Clustering method.')

    args = parser.parse_args()

    print("Parsing the general file...")
    general_path = os.path.join(STRUCTURES_ANNOTATION_DIR, "general.txt")
    chain_dict = parse_general_file(general_path)
    
    print("Listing PDB files...")
    pdb_files = sorted([os.path.join(args.pdb_folder, f) for f in os.listdir(args.pdb_folder) if f.endswith('.pdb')])
    
    print("Extracting sequences and calculating distances...")
    sequences = []
    
    for pdb_file in pdb_files:
        pdb_id = os.path.basename(pdb_file).split('.')[0]
        
        # Verifica si el pdb_id está en chain_dict
        if pdb_id in chain_dict:
            print(f"Processing {pdb_id}...")
            seqTCRA, seqTCRB, seqPep = extract_specific_sequences(pdb_file, chain_dict, extract_sequences)
        else:
            print(f"Warning: {pdb_id} not found in chain_dict. Using default chains.")
            chain_dict = {f'{pdb_id}': {'tcra_chain': 'D',
                                        'tcrb_chain': 'E',
                                        'peptide_chain': 'C',
                                        'b2m_chain': 'B',
                                        'mhc_chain': 'A'}}
            seqTCRA, seqTCRB, seqPep = extract_specific_sequences(pdb_file, chain_dict, extract_sequences)

        if not seqTCRA or not seqTCRB or not seqPep:
            print(f"Warning: Missing sequences for {pdb_id}. Skipping...")
            continue
        
        seqTCRA = seqTCRA.replace("X", "")
        seqTCRB = seqTCRB.replace("X", "")
        seqPep = seqPep.replace("X", "")
        
        cdr3a, cdr3b = '', ''
        if seqTCRA:
            anarci_output_a = run_anarci(seqTCRA, 'TRA')
            if anarci_output_a:
                cdr3a, _ = parse_CDR3(anarci_output_a)
        if seqTCRB:
            anarci_output_b = run_anarci(seqTCRB, 'TRB')
            if anarci_output_b:
                cdr3b, _ = parse_CDR3(anarci_output_b)

        sequences.append({
            'pdb_id': pdb_id,
            'cdr3a': cdr3a,
            'cdr3b': cdr3b,
            'peptide': seqPep})
    
    print("Calculating distance matrix...")
    n = len(sequences)
    distance_matrix = np.zeros((n, n))
    
    for i in range(n):
        for j in range(i + 1, n):
            dist_sum = get_distance_sum(
                sequences[i]['cdr3a'], sequences[j]['cdr3a'],
                sequences[i]['cdr3b'], sequences[j]['cdr3b'],
                sequences[i]['peptide'], sequences[j]['peptide'])
            distance_matrix[i, j] = dist_sum
            distance_matrix[j, i] = dist_sum

    print("Performing clustering...")
    clusters = cluster_pdbs(distance_matrix, distance_threshold=args.distance, method=args.method)

    print("Creating DataFrame with clusters...")
    pdb_ids = [seq['pdb_id'] for seq in sequences]
    df_clusters = vector_to_df(clusters, pdb_ids)
    
    print("Identifying non-redundant structures...")
    df_non_redundant = get_non_redundant_structures(df_clusters)
    pdb_nonred = df_non_redundant['pdb_id'].tolist()
    
    # Create a DataFrame for all sequences with 'nonred' column
    df_sequences = pd.DataFrame(sequences)
    df_sequences['cluster_id'] = df_clusters['cluster_id']
    df_sequences['nonred'] = df_sequences['pdb_id'].apply(lambda x: x in pdb_nonred)

    # Reorder columns and sort by 'cluster_id'
    df_sequences = df_sequences[['pdb_id', 'peptide', 'cdr3a', 'cdr3b', 'cluster_id', 'nonred']]
    df_sequences = df_sequences.sort_values(by='cluster_id')
    
    # Save to CSV
    df_sequences.to_csv(args.output, index=False)
    
    print(f"Summary DataFrame saved to '{args.output}'.")

    # Optional: Copy non-redundant PDB files to a new directory
    print("Copying non-redundant PDB files...")
    copy_non_redundant_pdbs(args.pdb_folder, args.nr_folder, df_non_redundant)
    print(f"Non-redundant PDB files copied to '{args.nr_folder}'.")

if __name__ == "__main__":
    main()