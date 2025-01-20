# Find contact map:

# This file contains functions to find the most similar TCR given an input TCR using tcrdist3. 
# 1) get_germlines (sequence)

# TCRdist
# 2) create_dataframe (alpha_outputs_list, beta_outputs_list, pdb_id_list)
# 3) find_closest_tcr(df, alpha_seq, beta_seq, tcr_name, one_value=True)

# Levenshtein distance
# 4) parse_cdrs(alfa_anarci_output, beta_anarci_output, epitope)
# 5) compute_pairwise_distances (cdrs_alpha1, cdrs_alpha2, cdrs_beta1, cdrs_beta2, epitope1, epitope2)
# 6) compute_combined_distance(cdrs_alpha1, cdrs_beta1, epitope1, cdrs_alpha2, cdrs_beta2, epitope2)
# 7) create_distance_matrix(input_df, tcr_alpha_sequence, tcr_beta_sequence, epitope)
# 8) get_min_combined_distance(df_results)

#Import libraries
from anarci import anarci
import pandas as pd
from tcrdist.repertoire import TCRrep
import numpy as np
import pickle
from utils import extract_specific_sequences, extract_sequences, parse_general_file
from mapping import run_anarci, parse_CDR3, global_alignment
from select_nr_set import calculate_sequence_distance

seq_dict=parse_general_file("./structures_annotation/general.txt")
    
def get_germlines(seq:str):
    '''
    Get the VJ germlines from TCRa or TCRb sequences
    '''
    input_seq = [('seq', seq)]
    try:
        results = anarci(input_seq, scheme="imgt", output=False, assign_germline=True)
        v_gene = results[1][0][0]['germlines']['v_gene'][0][1]
        j_gene = results[1][0][0]['germlines']['j_gene'][0][1]
    except:
        v_gene = 'NA'
        j_gene = 'NA'
    return v_gene, j_gene

### TCRdist 
    
def create_dataframe(alpha_outputs, beta_outputs, pdb_ids):
    """
    Creates a DataFrame with columns for CDR3, VJ germlines, and PDB ID for multiple alpha and beta chain outputs.
    """
    data = {
        'pdb_id': [],
        'cdr3_b_aa': [],
        'v_b_gene': [],
        'j_b_gene': [],
        'cdr3_a_aa': [],
        'v_a_gene': [],
        'j_a_gene': []
    }

    for anarci_output_alpha, anarci_output_beta, pdb_id in zip(alpha_outputs, beta_outputs, pdb_ids):
        # Parse ANARCI outputs
        cdr3_a, seq_a = parse_CDR3(anarci_output_alpha)
        cdr3_b, seq_b = parse_CDR3(anarci_output_beta)
        
        # Get germlines
        v_a, j_a = get_germlines(seq_a)
        v_b, j_b = get_germlines(seq_b)
        
        # Append to the data dictionary
        data['pdb_id'].append(pdb_id)
        data['cdr3_b_aa'].append(cdr3_b)
        data['v_b_gene'].append(v_b)
        data['j_b_gene'].append(j_b)
        data['cdr3_a_aa'].append(cdr3_a)
        data['v_a_gene'].append(v_a)
        data['j_a_gene'].append(j_a)

    # Create DataFrame
    df = pd.DataFrame(data)
    return df

def find_closest_tcr(df, alpha_seq, beta_seq, epitope_seq, tcr_name, one_value=True):
    """
    Finds the closest TCR to the given sequences.

    Parameters:
    - df (pd.DataFrame): The DataFrame containing existing TCR information.
    - alpha_seq (str): The TCR alpha chain sequence of the new TCR.
    - beta_seq (str): The TCR beta chain sequence of the new TCR.
    - epitope_seq (str): The epitope sequence for comparison.
    - tcr_name (str): The TCR name for the new entry.

    Returns:
    - list: The list of pdb_ids of the closest TCRs, ensuring `tcr_name` and `pdb_id` don't match.
    """

    pdb_id = tcr_name

    # Process the alpha chain
    anarci_output_alpha = run_anarci(alpha_seq, "TRA")
    cdr3_alpha, _ = parse_CDR3(anarci_output_alpha)
    v_gene_alpha, j_gene_alpha = get_germlines(alpha_seq)
    
    # Process the beta chain
    anarci_output_beta = run_anarci(beta_seq, "TRB")
    cdr3_beta, _ = parse_CDR3(anarci_output_beta)
    v_gene_beta, j_gene_beta = get_germlines(beta_seq)
    
    # Create a new row as DataFrame
    new_row = pd.DataFrame({
        'pdb_id': [pdb_id],
        'cdr3_a_aa': [cdr3_alpha],
        'v_a_gene': [v_gene_alpha],
        'j_a_gene': [j_gene_alpha],
        'cdr3_b_aa': [cdr3_beta],
        'v_b_gene': [v_gene_beta],
        'j_b_gene': [j_gene_beta],
        'count': [1]
    })

    # Create TCRrep for the new TCR row
    tr_current = TCRrep(cell_df=new_row, 
            organism='human', 
            chains=['alpha', 'beta'], 
            compute_distances=True,
            db_file='alphabeta_gammadelta_db.tsv')
    
    # To search in different datasets
    #length=len(epitope_seq)
    #pdb_ids_samelength=[]
    #for pdb_id in df['pdb_id']:
    #    pdb_path=f"./pdb_files/{pdb_id}.pdb"
    #    a, b, e = extract_specific_sequences(pdb_path, seq_dict, extract_sequences)
    #    if len(e)==length:
    #        pdb_ids_samelength.append(pdb_id)

    # To search in this dataset
    with open('./structures_annotation/epitope_dict.pkl', 'rb') as f:
        epitope_dict = pickle.load(f)
    epitope_length = len(epitope_seq)
    pdb_entries_samelength = epitope_dict.get(f"{epitope_length}", [])
    pdb_ids_samelength = [pdb_id for pdb_id, _ in pdb_entries_samelength]
    # Filter DataFrame for PDBs with the same epitope length
    df_samelength = df[df['pdb_id'].isin(pdb_ids_samelength)]
    # Avoid the current TCR
    df_samelength = df_samelength[df_samelength['pdb_id'] != pdb_id]

    # Compute TCRrep for the filtered PDBs
    tr = TCRrep(cell_df=df_samelength, 
            organism='human', 
            chains=['alpha', 'beta'], 
            compute_distances=True,
            db_file='alphabeta_gammadelta_db.tsv')
    
    # Compute distances between the new TCR and existing TCRs
    tr.compute_rect_distances(df=tr_current.clone_df, df2=tr.clone_df)

    # Sum the alpha and beta chain distances to get a global distance
    tr.alpha_beta = tr.rw_alpha + tr.rw_beta
    array = tr.alpha_beta
    pdbs = tr.clone_df['pdb_id'].values
    min_value = np.min(array)
    min_indices = np.where(array == min_value)[1]
    min_pdb_ids = pdbs[min_indices].tolist()
    
    #In other datasets
    #distances={}
    #for pdb_id in min_pdb_ids:
    #    pdb_path=f"./pdb_files/{pdb_id}.pdb"
    #    a, b, e = extract_specific_sequences(pdb_path, seq_dict, extract_sequences) #Usa el dict aquí apra no iterar
    #    dist=calculate_sequence_distance(epitope_seq, e)
    #    distances[pdb_id]=dist

    distances = {}
    # Calculate the epitope distance for each PDB with the minimum distance (this dataset)
    for pdb_id in min_pdb_ids:
        sequence = next((seq for pid, seq in epitope_dict.get(f"{epitope_length}", []) if pid == pdb_id), None)
        if sequence:
            dist = calculate_sequence_distance(epitope_seq, sequence)
            distances[pdb_id] = dist
    # Select the PDB with the minimum epitope distance
    min_pdb_ids = [k for k, v in distances.items() if v == min(distances.values())]
    # If there are more than one value, select the first one
    if one_value == True:
        min_pdb_ids = min_pdb_ids[0]
        return min_pdb_ids
    else:
        return min_pdb_ids
    
# Levenshtein distance

def parse_cdrs(alfa_anarci_output, beta_anarci_output, epitope):
    """
    Extracts CDR1, CDR2, CDR2.5, and CDR3 sequences from ANARCI output for both alpha and beta chains,
    and includes the epitope in the results.
    
    Parameters:
    - alfa_anarci_output (str): The ANARCI output containing the CDR information for the alpha chain.
    - beta_anarci_output (str): The ANARCI output containing the CDR information for the beta chain.
    - epitope (str): The epitope associated.
    
    Returns:
    - dict: A dictionary containing CDR sequences (CDR1, CDR2, CDR2.5, and CDR3), and epitope for both alpha and beta chains.
    """
    # Initialize dictionaries to store sequences
    cdrs_alpha = {'CDR1': '', 'CDR2': '', 'CDR2.5': '', 'CDR3': '', 'FullSeq': ''}
    cdrs_beta = {'CDR1': '', 'CDR2': '', 'CDR2.5': '', 'CDR3': '', 'FullSeq': ''}
    
    # Process alpha chain
    for line in alfa_anarci_output.split("\n"):
        if line and line[0] != "#" and line[0] != "/":
            parts = line.rstrip().split()
            if len(parts) >= 3:
                _, num, *residues = parts
                num = int(num)  # Position number
                res = residues[-1]  # Amino acid residue
                
                if res != "-":
                    cdrs_alpha['FullSeq'] += res
                    
                    # CDR1: Residue numbers 27-38
                    if 27 <= num <= 38:
                        cdrs_alpha['CDR1'] += res
                    
                    # CDR2: Residue numbers 56-65
                    if 56 <= num <= 65:
                        cdrs_alpha['CDR2'] += res
                    
                    # CDR2.5: Residue numbers 81-86
                    if 81 <= num <= 86:
                        cdrs_alpha['CDR2.5'] += res
                    
                    # CDR3: Residue numbers 105-117
                    if 105 <= num <= 117:
                        cdrs_alpha['CDR3'] += res
    
    # Process beta chain
    for line in beta_anarci_output.split("\n"):
        if line and line[0] != "#" and line[0] != "/":
            parts = line.rstrip().split()
            if len(parts) >= 3:
                _, num, *residues = parts
                num = int(num)  # Position number
                res = residues[-1]  # Amino acid residue
                
                if res != "-":
                    cdrs_beta['FullSeq'] += res
                    
                    # CDR1: Residue numbers 27-38
                    if 27 <= num <= 38:
                        cdrs_beta['CDR1'] += res
                    
                    # CDR2: Residue numbers 56-65
                    if 56 <= num <= 65:
                        cdrs_beta['CDR2'] += res
                    
                    # CDR2.5: Residue numbers 81-86
                    if 81 <= num <= 86:
                        cdrs_beta['CDR2.5'] += res
                    
                    # CDR3: Residue numbers 105-117
                    if 105 <= num <= 117:
                        cdrs_beta['CDR3'] += res

    # Return combined results with PDB ID and epitope
    return {
        'Alpha': cdrs_alpha,
        'Beta': cdrs_beta,
        'Epitope': epitope
    }

def compute_pairwise_distances(cdrs_alpha1, cdrs_alpha2, cdrs_beta1, cdrs_beta2, epitope1, epitope2):
    """Calculate pairwise distances between sequences of two TCRs."""
    distances_alpha = []
    distances_beta = []
    distances_epitope=[]
    # Calculate pairwise distances for CDRs
    for cdr_name in ['CDR1', 'CDR2', 'CDR2.5', 'CDR3']:
        # Alpha chain
        seq1_alpha = cdrs_alpha1[cdr_name]
        seq2_alpha = cdrs_alpha2[cdr_name]
        distance_alpha = calculate_sequence_distance(seq1_alpha, seq2_alpha)
        distances_alpha.append(distance_alpha)
        
        # Beta chain
        seq1_beta = cdrs_beta1[cdr_name]
        seq2_beta = cdrs_beta2[cdr_name]
        distance_beta = calculate_sequence_distance(seq1_beta, seq2_beta)
        distances_beta.append(distance_beta)
        
    
    # Calculate pairwise distance for epitopes
    distance_epitope = calculate_sequence_distance(epitope1, epitope2)
    distances_epitope.append(distance_epitope)
    
    return distances_alpha, distances_beta, distances_epitope

def compute_combined_distance(cdrs_alpha1, cdrs_beta1, epitope1, cdrs_alpha2, cdrs_beta2, epitope2):
    """Compute a combined metric of Levenshtein distances between sequences of two TCRs."""
    distances_alpha, distances_beta, distances_epitope = compute_pairwise_distances(cdrs_alpha1, cdrs_alpha2, cdrs_beta1, cdrs_beta2, epitope1, epitope2)
    combined_distance = sum(distances_alpha+distances_beta+distances_epitope)  # Sum up all the distances for a combined score
    return distances_alpha, distances_beta, distances_epitope, combined_distance

def create_distance_matrix(input_df, tcr_alpha_sequence, tcr_beta_sequence, epitope):
    results = []

    anarci_output_alpha = run_anarci(tcr_alpha_sequence, 'alpha')
    anarci_output_beta = run_anarci(tcr_beta_sequence, 'beta')
    
    cdrs = parse_cdrs(anarci_output_alpha, anarci_output_beta, epitope)
    cdrs_alpha = cdrs['Alpha']
    cdrs_beta = cdrs['Beta']
    epitope_input = cdrs['Epitope']
    
    
    for _, row in input_df.iterrows():
        pdb_id = row['pdb_id']
        
        cdrs_alpha_row = {
            'CDR1': row['cdr1_a'],
            'CDR2': row['cdr2_a'],
            'CDR2.5': row['cdr2.5_a'],
            'CDR3': row['cdr3_a']
        }
        
        cdrs_beta_row = {
            'CDR1': row['cdr1_b'],
            'CDR2': row['cdr2_b'],
            'CDR2.5': row['cdr2.5_b'],
            'CDR3': row['cdr3_b']
        }
        
        epitope_row = row['epitope']
        
        distances_alpha, distances_beta, distances_epitope, combined_distance = compute_combined_distance(
            cdrs_alpha, cdrs_beta, epitope_input,
            cdrs_alpha_row, cdrs_beta_row, epitope_row
        )
        
        results.append({
            'pdb_id': pdb_id,
            'distance_CDR1_alpha': distances_alpha[0],
            'distance_CDR2_alpha': distances_alpha[1],
            'distance_CDR2.5_alpha': distances_alpha[2],
            'distance_CDR3_alpha': distances_alpha[3],
            'distance_CDR1_beta': distances_beta[0],
            'distance_CDR2_beta': distances_beta[1],
            'distance_CDR2.5_beta': distances_beta[2],
            'distance_CDR3_beta': distances_beta[3],
            'distance_epitope': distances_epitope[0],
            'combined_distance': combined_distance
        })
    
    results_df = pd.DataFrame(results)
    return results_df

def get_min_combined_distance(df_results):
    min_idx = df_results['combined_distance'].idxmin()
    min_row = df_results.loc[min_idx]
    min_pdb_id = min_row['pdb_id']
    return min_pdb_id


def find_closest_tcr2(df, alpha_seq, beta_seq, epitope_seq, tcr_name, mhc_allele, mhc_seq=None, structure=True, one_value=True):
    pdb_id = tcr_name

    # Process alpha chain
    anarci_output_alpha = run_anarci(alpha_seq, "TRA")
    cdr3_alpha, _ = parse_CDR3(anarci_output_alpha)
    v_gene_alpha, j_gene_alpha = get_germlines(alpha_seq)

    # Process beta chain
    anarci_output_beta = run_anarci(beta_seq, "TRB")
    cdr3_beta, _ = parse_CDR3(anarci_output_beta)
    v_gene_beta, j_gene_beta = get_germlines(beta_seq)

    # Create a new TCR entry
    new_row = pd.DataFrame({
        'pdb_id': [pdb_id],
        'cdr3_a_aa': [cdr3_alpha],
        'v_a_gene': [v_gene_alpha],
        'j_a_gene': [j_gene_alpha],
        'cdr3_b_aa': [cdr3_beta],
        'v_b_gene': [v_gene_beta],
        'j_b_gene': [j_gene_beta],
        'count': [1]})

    # Compute TCRrep for the new TCR
    tr_current = TCRrep(cell_df=new_row, organism='human', chains=['alpha', 'beta'], 
                        compute_distances=True, db_file='alphabeta_gammadelta_db.tsv')
    
    # Load epitope dictionary and filter MHC alleles
    with open('./structures_annotation/epitope_dict.pkl', 'rb') as f:
        epitope_dict = pickle.load(f)
    epitope_length = len(epitope_seq)

    # Filter by MHC allele
    alleles_df = pd.read_csv('./structures_annotation/mhc_annotation.csv')
    pdb_ids_allele = alleles_df.loc[alleles_df['mhci_allele'] == mhc_allele, 'pdb_id'].values
    df_samelength = df[df['pdb_id'].isin(pdb_ids_allele) & (df['pdb_id'] != pdb_id)]

    if df_samelength.empty:
        print("No TCRs with the same MHC allele")
        print("Searching in the whole dataset")
        
        # Extract MHC sequence if needed
        if structure:
            seqs_query = extract_sequences(f"./pdb_files/{pdb_id}.pdb")
            mhc_seq_query_id = seq_dict[pdb_id]['mhc_chain']
            mhc_seq_query = seqs_query[mhc_seq_query_id]
        elif mhc_seq:
            mhc_seq_query = mhc_seq
        else:
            mhc_seqs_df = pd.read_csv('./input/input_MHCs2.csv')
            mhc_allele = mhc_allele + ":01:01"  
            mhc_seq_query = mhc_seqs_df[mhc_seqs_df['mhc_allele'] == mhc_allele]['mhc_seq'].values[0]
        
        
        # Perform global alignment with all MHC sequences
        scores = {}
        for pdb_file in os.listdir('./pdb_files'):
            if pdb_file.endswith('.pdb') and pdb_file != f"{pdb_id}.pdb":
                pdb_idmhc = pdb_file.split('.')[0]
                mhc_seqid = seq_dict[pdb_idmhc]['mhc_chain']
                pdb_path = f"./pdb_files/{pdb_file}"
                seqs = extract_sequences(pdb_path)
                mhc_seq = seqs[mhc_seqid]
                _, _, score = global_alignment(mhc_seq_query, mhc_seq) 
                scores[pdb_idmhc] = score

        # Select PDBs with the highest score
        max_score = max(scores.values())
        max_pdb_ids = [pdb_id for pdb_id, score in scores.items() if score == max_score]


        # Collect MHC alleles for the selected PDBs
        mhc_alleles_max = []
        for pdb_id in max_pdb_ids:
            mhc_seqid = seq_dict[pdb_id]['mhc_chain']
            pdb_path = f"./pdb_files/{pdb_id}.pdb"
            seqs = extract_sequences(pdb_path)
            mhc_seq = seqs[mhc_seqid]
            mhc_allele = alleles_df[alleles_df['pdb_id'] == pdb_id]['mhci_allele'].values
            mhc_alleles_max.append(mhc_allele)
        
        # Filter the DataFrame for the selected MHC alleles
        pdb_ids = []
        for mhc_allele in mhc_alleles_max:
            if isinstance(mhc_allele, (list, np.ndarray)):
                pdb_ids_allele = alleles_df[alleles_df['mhci_allele'].isin(mhc_allele)]['pdb_id'].values
            else:
                pdb_ids_allele = alleles_df[alleles_df['mhci_allele'] == mhc_allele]['pdb_id'].values
            pdb_ids.extend(pdb_ids_allele)

        df_samelength = df[df['pdb_id'].isin(pdb_ids)]
        df_samelength = df_samelength[df_samelength['pdb_id'] != pdb_id]
        
    # Further filter by epitope length
    samelength = [row['pdb_id'] for _, row in df_samelength.iterrows() if epitope_length == len(extract_specific_sequences(f"./pdb_files/{row['pdb_id']}.pdb", seq_dict, extract_sequences)[2])]
    df_samelength = df_samelength[df_samelength['pdb_id'].isin(samelength)]
    
    if df_samelength.empty:
        # To search in this dataset
        print("There are not pdb_files with same or similar MHC alleles and same epitope length doing reverse")
        with open('./structures_annotation/epitope_dict.pkl', 'rb') as f:
            epitope_dict = pickle.load(f)
        epitope_length = len(epitope_seq)
        pdb_entries_samelength = epitope_dict.get(f"{epitope_length}", [])
        pdb_ids_samelength = [pdb_id for pdb_id, _ in pdb_entries_samelength]  
        # Exclude current TCR_id from pdb_ids_samelength
    
        for pdb_id in pdb_ids_samelength:
            if pdb_id == tcr_name:
                pdb_ids_samelength.remove(pdb_id)


        # Extract MHC sequence 
        if structure:
            seqs_query = extract_sequences(f"./pdb_files/{pdb_id}.pdb")
            mhc_seq_query_id = seq_dict[pdb_id]['mhc_chain']
            mhc_seq_query = seqs_query[mhc_seq_query_id]
        elif mhc_seq:
            mhc_seq_query = mhc_seq
        else:
            mhc_seqs_df = pd.read_csv('./input/input_MHCs2.csv')
            mhc_allele = mhc_allele + ":01:01"  
            mhc_seq_query = mhc_seqs_df[mhc_seqs_df['mhc_allele'] == mhc_allele]['mhc_seq'].values[0]

        # Filter for MHC similarity.
        scores_len={}
        for len_pdb_id in pdb_ids_samelength:
            if len_pdb_id != pdb_id:
                path_len=f"./pdb_files/{len_pdb_id}.pdb"
                mhc_seqid_len = seq_dict[len_pdb_id]['mhc_chain']
                seqs_len = extract_sequences(path_len)
                mhc_seq_len = seqs_len[mhc_seqid_len]
                aligned_seq_len, aligned_seq_query_len, score_len = global_alignment(mhc_seq_query, mhc_seq_len) 
                scores_len[len_pdb_id] = score_len
        
        max_score_len = max(scores_len.values())
        max_pdb_ids_len = [len_pdb_id for len_pdb_id, score_len in scores_len.items() if score_len == max_score_len]

        mhc_alleles_max_len = []
        for len_pdb_id in max_pdb_ids_len:
            mhc_seqid_len = seq_dict[len_pdb_id]['mhc_chain']
            pdb_path_len = f"./pdb_files/{len_pdb_id}.pdb"
            seqs_len = extract_sequences(pdb_path_len)
            mhc_seq_len = seqs_len[mhc_seqid_len]
            mhc_allele_len = alleles_df[alleles_df['pdb_id'] == len_pdb_id]['mhci_allele'].values
            mhc_alleles_max_len.append(mhc_allele_len)
        
        pdb_ids_len = []
        for mhc_allele_len in mhc_alleles_max_len:
            # Si mhc_allele_len es un array/lista, utiliza .isin()
            if isinstance(mhc_allele_len, (list, np.ndarray)):
                pdb_ids_allele_len = alleles_df[alleles_df['mhci_allele'].isin(mhc_allele_len)]['pdb_id'].values
            else:
                # Si es un único valor, filtra directamente
                pdb_ids_allele_len = alleles_df[alleles_df['mhci_allele'] == mhc_allele_len]['pdb_id'].values
            pdb_ids_len.extend(pdb_ids_allele_len)
            
        pdb_ids_len = [pdb_id for pdb_id in pdb_ids_len if pdb_id in pdb_ids_samelength] 
        df_samelength = df[df['pdb_id'].isin(pdb_ids_len)]
        df_samelength = df_samelength[df_samelength['pdb_id'] != pdb_id]
    # Compute TCRrep for the filtered PDBs
    tr = TCRrep(cell_df=df_samelength, organism='human', chains=['alpha', 'beta'], 
                compute_distances=True, db_file='alphabeta_gammadelta_db.tsv')

    tr.compute_rect_distances(df=tr_current.clone_df, df2=tr.clone_df)
    tr.alpha_beta = tr.rw_alpha + tr.rw_beta
    min_value = np.min(tr.alpha_beta)
    min_indices = np.where(tr.alpha_beta == min_value)[1]
    min_pdb_ids = tr.clone_df['pdb_id'].values[min_indices].tolist()

    distances = {
        pdb_id: calculate_sequence_distance(epitope_seq, extract_specific_sequences(f"./pdb_files/{pdb_id}.pdb", seq_dict, extract_sequences)[2])
        for pdb_id in min_pdb_ids
    }

    min_pdb_ids = [k for k, v in distances.items() if v == min(distances.values())]
    return min_pdb_ids[0] if one_value else min_pdb_ids

