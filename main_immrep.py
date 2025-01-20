# Main function to rank TCRs given an input neoantigenic sequence.
# EXAMPLE: 
# python main_positions_fast.py -cm ./contact_maps -t ./input/input_test/1mi5_df.csv -m ./input/input_MHCs.csv -e "FLRGRAYGL" -a "B*08:01" -mp ./model/TCR_MHC_potential.csv -tp ./model/TCR-p -g ./structures_annotation/general.txt -s scores_example.csv 

import argparse
import pandas as pd
import os
import warnings  
import ast
import glob
from concurrent.futures import ProcessPoolExecutor
import pickle

from utils import *
from find_contact_map import *
from potential_calc import get_potential
from extract_contacts import filter_contacts
from mapping import run_anarci_parallel, add_imgt_mappings, map_epitope_residue, parse_anarci_output, extract_residues_and_resids, map_imgt_to_original, global_alignment, map_alignment_to_residues

warnings.filterwarnings("ignore")

def folder(arg):
    """ Custom type for argparse to handle folder paths. """
    if not os.path.isdir(arg):
        raise argparse.ArgumentTypeError(f"The folder {arg} does not exist.")
    return arg

def process_tcr(tcr_id, alpha_seq, beta_seq, epitope_seq, tcr_p_potential, tcr_mhc_potential, mhc_df, chain_dict, verbose, mhc_allele=None):
    """ Process a single TCR and compute scores. """
    result_string = ""
    result_string += f"\n{'-'*40}\n------ Processed TCR: {tcr_id} ------\n{'-'*40}\n"
    results = {
        "tcr_id": None,
        "mhc_allele": None,
        "score_tcr_p1":None,
        "score_tcr_p2": None,
        "score_tcr_p3": None,
        "score_tcr_p4": None,
        "score_tcr_p5": None,
        "score_tcr_p6": None,
        "score_tcr_p7": None,
        "score_tcr_p8": None,
        "score_tcr_p9": None,
        "score_tcr_all":None,
        "score_tcr_mhc": None}

    # Processing similar TCRs and extracting pdb file if exists
    pdb_id_similar = None
    chains = {} 
    alpha_seq=str(alpha_seq)
    beta_seq=str(beta_seq)
    epitope_seq=str(epitope_seq)
    tcr_id=str(tcr_id)
    try:
        tcr_dist=pd.read_csv("./structures_annotation/TCRdist_df.csv")
        mhc_original = mhc_allele.replace(":01:01", "")
        similar_tcr = find_closest_tcr2(tcr_dist, alpha_seq, beta_seq, epitope_seq, str(tcr_id), mhc_original, mhc_seq=None, structure=False, one_value=True)
        pdb_id_similar = similar_tcr       
    except Exception as e:
        result_string += f"\nNo similar TCRs found for {tcr_id}. {e}\n"

    # If pdb found, process contacts and mappings
    if pdb_id_similar:
        pdb_file_path =  f"./contact_maps/{pdb_id_similar}_contacts.csv"
        if os.path.isfile(pdb_file_path):
            contacts_df = pd.read_csv(pdb_file_path)
            chains = chain_dict.get(pdb_id_similar, {})
            if all(chains.values()):  
                contacts_TCR_p, contacts_TCR_MHC = filter_contacts(
                    contacts_df,
                    chains['tcra_chain'],
                    chains['tcrb_chain'],
                    chains['peptide_chain'],
                    chains['mhc_chain'],
                    threshold=1,
                    remove_X=False)

        # Mapping with IMGT numbering
        result_string += "\n Sequences reenumbered with IMGT convention\n"
        alpha_pdb, beta_pdb, epitope_pdb = extract_specific_sequences(f"./pdb_files/{pdb_id_similar}.pdb", chain_dict, extract_sequences)
        sequences = [alpha_pdb, beta_pdb]
        chain_ids = ['alpha', 'beta']
        anarci_results = run_anarci_parallel(sequences, chain_ids)
        parsed_anarci_D = parse_anarci_output(anarci_results['alpha'])
        parsed_anarci_E = parse_anarci_output(anarci_results['beta'])
        residues_D = extract_residues_and_resids(f"./pdb_files/{pdb_id_similar}.pdb", chains['tcra_chain'])
        residues_E = extract_residues_and_resids(f"./pdb_files/{pdb_id_similar}.pdb", chains['tcrb_chain'])
        mapping_D = map_imgt_to_original(parsed_anarci_D, residues_D)
        mapping_E = map_imgt_to_original(parsed_anarci_E, residues_E)
        imgt_mappings = {pdb_id_similar: {chains['tcra_chain']: mapping_D, chains['tcrb_chain']: mapping_E}}
        contacts_imgt_M = add_imgt_mappings(contacts_TCR_MHC, imgt_mappings)
        contacts_imgt_P = add_imgt_mappings(contacts_TCR_p, imgt_mappings)

        # Process input TCR
        sequences_tcr = [alpha_seq, beta_seq]
        chain_ids_tcr = ['alpha', 'beta']
        result_string += f"\nMapping TCR {tcr_id} into {pdb_id_similar}\n"
        anarci_results_tcr = run_anarci_parallel(sequences_tcr, chain_ids_tcr)
        parsed_tcr_D = parse_anarci_output(anarci_results_tcr['alpha'])
        parsed_tcr_E = parse_anarci_output(anarci_results_tcr['beta'])
        imgt_dict_D = dict(parsed_tcr_D)
        imgt_dict_E = dict(parsed_tcr_E)

        contacts_TCR_p[tcr_id] = contacts_TCR_p.apply(lambda row: imgt_dict_D.get(row['imgt_from'], None) if row['chain_from'] == chains['tcra_chain']
                            else imgt_dict_E.get(row['imgt_from'], None) if row['chain_from'] == chains['tcrb_chain']
                            else None, axis=1)
            
        contacts_TCR_MHC[tcr_id] = contacts_TCR_MHC.apply(lambda row: imgt_dict_D.get(row['imgt_from'], None) if row['chain_from'] == chains['tcra_chain']
                            else imgt_dict_E.get(row['imgt_from'], None) if row['chain_from'] == chains['tcrb_chain']
                            else None, axis=1)

        ##### PROCESSING EPITOPE #####
        result_string += f"\n-> Processed Epitope: {epitope_seq}\n"
        contacts_TCR_p['epitope'] = contacts_TCR_p.apply(lambda row: map_epitope_residue(row, epitope_seq), axis=1)

        # Add TCR-Potential (looping for P1-P9)
        result_string += "\nCalculated TCR-peptide potential\n"
        total_scores = {}
        
        # Dictionary of assignments for different peptide lengths
        assignments_dict = {
            4: {1: 1, 2: 2, 3: 3, 4: 9},
            8: {1: 1, 2: 2, 3: 3, 4: 5, 5: 5, 6: 7, 7: 8, 8: 9},
            9: {1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9},
            10: {1: 1, 2: 2, 3: 3, 4: 4, 5: 6, 6: 5, 7: 6, 8: 6, 9: 8, 10: 9},
            11: {1: 1, 2: 2, 3: 3, 4: 4, 5: 4, 6: 4, 7: 8, 8: 8, 9: 7, 10: 8, 11: 9},
            12: {1: 1, 2: 2, 3: 3, 4: 2, 5: 5, 6: 7, 7: 7, 8: 6, 9: 9, 10: 9, 11: 9, 12: 9},
            13: {1: 1, 2: 2, 3: 3, 4: 6, 5: 3, 6: 4, 7: 8, 8: 6, 9: 8, 10: 9, 11: 9, 12: 9, 13: 9}}

        # Select the appropriate assignments dictionary based on the peptide length
        peptide_length = len(epitope_seq)
        assignments = assignments_dict.get(peptide_length)

        # Raise an error if the peptide length is unsupported
        if not assignments:
            raise ValueError(f"Unsupported peptide length: {peptide_length}")

        # Calculate potential scores and totals for each position
        for i in range(1, 10):
            potential_column = f'potential_P{i}'
            
            # Apply the potential function based on the residue assignment
            contacts_TCR_p[potential_column] = contacts_TCR_p.apply(
                lambda row: get_potential(row, tcr_p_potential[f'tcr_p_potential_P{i}'], tcr_id, "epitope")
                if assignments.get(row['resid_to']) == i else None, axis=1)
            
            # Sum the potential scores for each position
            total_scores[potential_column] = contacts_TCR_p[potential_column].sum(skipna=True)
            
            # Add the total score to the result string
            result_string += f"Total score for P{i}: {total_scores[potential_column]}\n"

        # Calculate total score for all positions
        total_score_all = sum(total_scores.values())
        result_string += f"Total score for all positions: {total_score_all}\n"

        #### PROCESSING MHC ####
        if mhc_allele:
            reference_allele = mhc_allele
            mhc_seq = mhc_df[mhc_df['mhc_allele'] == reference_allele]['mhc_seq'].values[0]

            result_string += f"\n-> Processing MHC-I allele: {reference_allele}\n"
                
            seq_pdb = extract_sequences(f"./pdb_files/{pdb_id_similar}.pdb")
            aligned_seq_pdb, aligned_seq_query, score = global_alignment(seq_pdb[chains['mhc_chain']], mhc_seq)
            residues_M = extract_residues_and_resids(f"./pdb_files/{pdb_id_similar}.pdb", chains['mhc_chain']) 
            mapped_residues = map_alignment_to_residues(aligned_seq_pdb, aligned_seq_query, residues_M)
                
            df_tuples = pd.DataFrame(mapped_residues, columns=['resid', 'mhc_pdb', reference_allele])
            contacts_TCR_MHC_updated = pd.merge(
                    contacts_TCR_MHC, 
                    df_tuples[['resid', reference_allele]], 
                    left_on='resid_to', 
                    right_on='resid', 
                    how='left')
                
            contacts_TCR_MHC_updated = contacts_TCR_MHC_updated.drop(columns=['resid'])

            # Add TCR-MHC potential
            result_string += "\nCalculating TCR-MHC-I potential\n"
            contacts_TCR_MHC_updated['potential'] = contacts_TCR_MHC_updated.apply(
                    lambda row: get_potential(row, tcr_mhc_potential, tcr_id, reference_allele), axis=1)
                
            contacts_TCR_MHC_updated['potential'] = pd.to_numeric(contacts_TCR_MHC_updated['potential'], errors='coerce')
            total_score_mhc = contacts_TCR_MHC_updated['potential'].sum()
            result_string +=f"Total score for TCR-MHC {reference_allele}: {total_score_mhc}"
        else:
            result_string += "\nNo MHC allele provided. Skipping MHC processing.\n"
            total_score_mhc = None

        # Append results to the results container
        results["tcr_id"]=tcr_id
        results["mhc_allele"]=reference_allele if mhc_allele else None
        for i in range(1, 10):
            results[f"score_tcr_p{i}"]=total_scores.get(f"potential_P{i}", 0)
        results["score_tcr_all"]=total_score_all
        results["score_tcr_mhc"]=total_score_mhc if total_score_mhc else None

    else:
        result_string += f"\nNo similar TCR found for input TCR. Skipping {tcr_id}.\n"
    
    if verbose:
        print(result_string)

    return results

def main():
    parser = argparse.ArgumentParser(description='Process input TCRs and rank them based on statistic potential.')
    parser.add_argument("-t", "--tcr_file", type=str, required=True, help='CSV file with TCRs to process. Format: tcr_id,alpha_seq,beta_seq')
    parser.add_argument("-e", "--epitope_seq", type=str, required=True, help='Epitope df')
    parser.add_argument("-a", "--mhc_allele", type=str, required=False, help='Annotation df')
    parser.add_argument("-m", "--mhc_file", type=str, required=True, help='CSV file with MHCs to process. Format: mhc_allele,mhc_seq')
    parser.add_argument("-tp", "--tcr_potential_folder", type=folder, required=True, help="Directory containing TCR-peptide potential files. Format: TCR_p_potential_PX.csv. X=1:9")
    parser.add_argument("-mp", "--mhc_potential", type=str, required=True, default='./model/TCR_MHC_potential.csv', help='CSV file with TCR-MHC potential. Format: residue_from,residue_to,potential')
    parser.add_argument("-s", "--scores_filename", type=str, required=False, default='scores.csv', help='Output CSV file for scores (default: scores.csv)')
    parser.add_argument("-v", "--verbose", action="store_false", default=True, help='Verbose mode')

    args = parser.parse_args()


    # Create output directory if it doesn't exist
    output_dir = "output"
    os.makedirs(output_dir, exist_ok=True)
    
    # Load data
    tcr_df = pd.read_csv(args.tcr_file)
    mhc_df = pd.read_csv(args.mhc_file)
    epitope_df = pd.read_csv(args.epitope_seq)
    epitope_df['tcr_id'] = epitope_df['tcr_id'].astype(str)
    mhc_annotation = pd.read_csv(args.mhc_allele) 
    mhc_annotation['tcr_id'] = mhc_annotation['tcr_id'].astype(str)

    print("Loading models...")
    csv_files = glob.glob(os.path.join(args.tcr_potential_folder, "*_P*.csv"))
    tcr_p_potential = {}
    if not csv_files:
        print(f"No CSV files found in the directory: {args.tcr_potential_folder}")
    else:
        for file in csv_files:
            try:
                file_name = os.path.basename(file)
                potential_number = file_name.split('_P')[1].split('.csv')[0]
                potential_key = f'tcr_p_potential_P{potential_number}'
                tcr_p_potential[potential_key] = pd.read_csv(file)
            except Exception as e:
                print(f"Failed to load {file}: {e}")
                
    tcr_mhc_potential = pd.read_csv(args.mhc_potential)
    chain_dict = parse_general_file("./structures_annotation/general.txt")
    similarity_df = pd.read_csv("./structures_annotation/closest_tcr.csv")
    
    # Process each TCR in parallel
    results = []
    with ProcessPoolExecutor(max_workers=5) as executor:
        for tcr_id, alpha_seq, beta_seq in zip(tcr_df['tcr_id'], tcr_df['alpha_seq'], tcr_df['beta_seq']):
            #Exclude TCR_id 182 and 21960
            if tcr_id != "182" and tcr_id != "21960":
                try:
                    # Allele selection 
                    mhc_allele = mhc_annotation[mhc_annotation['tcr_id'] == str(tcr_id)]['mhc_allele'].values[0]
                    mhc_allele = mhc_allele + ":01:01"
                    print(f"Processing TCR {tcr_id} with MHC allele {mhc_allele}...")
                    epitope_seq = epitope_df[epitope_df['tcr_id'] == str(tcr_id)]['epitope'].values[0]
                    print(f"Processing TCR {tcr_id} with epitope {epitope_seq}...")
                    result = executor.submit(process_tcr, tcr_id, alpha_seq, beta_seq, epitope_seq, tcr_p_potential, tcr_mhc_potential, mhc_df, chain_dict, args.verbose, mhc_allele)
                    results.append(result)
                except Exception as e:
                    print(f"Error processing TCR {tcr_id}: {e}")
            else:
                print(f"Skipping TCR {tcr_id}...")
                
    # Collect and save results
    results_list = [r.result() for r in results]
    results_df = pd.DataFrame(results_list)
    output_file = os.path.join("output", args.scores_filename)
    results_df.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}")
    
if __name__ == "__main__":
    main()
