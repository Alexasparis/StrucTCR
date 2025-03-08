#!/usr/bin/env python3

# This script is used to give an score to the binding of a query TCR to a query pMHC complex.


# Example of execution: python main.py -i input/input_df.csv -e "LLFGYPVYV" -a "A*02:01" -tcrp ./model/TCR-p_all/ -tcrm ./model/TCR_MHC_all.csv -o ./output/scores2.csv -w 5
import argparse
import pandas as pd
import os  
import ast
from concurrent.futures import ProcessPoolExecutor
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

from utils import parse_general_file, extract_sequences, extract_specific_sequences
from potential_calc import get_potential
from extract_contacts import filter_contacts
from mapping import run_anarci, add_imgt_mappings, map_epitope_residue, parse_anarci_output, extract_residues_and_resids, map_imgt_to_original, global_alignment, map_alignment_to_residues
from find_contact_map import find_closest_tcr
from config import DATA_DIR,STRUCTURES_ANNOTATION_DIR, CONTACT_MAPS_DIR, PDB_DIR

def process_tcr(tcr_id, alpha_seq, beta_seq, epitope_seq, tcr_p_potential, tcr_mhc_potential, mhc_df, chain_dict, verbose, mhc_allele):
    """ Process a TCR and compute scores."""
    result_string = "" # Verbose
    result_string += f"\n{'-'*40}\n------ Processed TCR: {tcr_id} ------\n{'-'*40}\n"
    peptide_length = len(epitope_seq)
    results = {"tcr_id": None,
                **{f"score_tcr_p{i+1}": None for i in range(peptide_length)},
                "score_tcr_all": None,
                "score_tcr_mhc": None,}

    # Processing similar TCRs and extracting pdb file if exists
    similar_tcr = None
    chains = {} 
    try:
        tcr_dist=pd.read_csv(os.path.join(STRUCTURES_ANNOTATION_DIR, "TCRdist_df.csv"))
        similar_tcr = find_closest_tcr(tcr_dist, alpha_seq, beta_seq, epitope_seq, str(tcr_id), mhc_allele, mhc_seq=None, structure=False, one_value=True)
        result_string += f"\nThe most similar TCR to {tcr_id} is {similar_tcr}."       
    except Exception as e:
        result_string += f"\nNo similar TCRs found for {tcr_id}. {e}\n"

    # If pdb found, process contacts and mappings
    if similar_tcr:
        pdb_cm_path =  os.path.join(CONTACT_MAPS_DIR, f"{similar_tcr}_contacts.csv")
        pdb_file_path = os.path.join(PDB_DIR, f"{similar_tcr}.pdb")
        if os.path.isfile(pdb_cm_path):
            contacts_df = pd.read_csv(pdb_cm_path)
            chains = chain_dict.get(similar_tcr, {'tcra_chain': 'D',
                                                  'tcrb_chain': 'E',
                                                  'peptide_chain': 'C',
                                                  'b2m_chain': 'B',
                                                  'mhc_chain': 'A'})
            if all(chains.values()):  
                contacts_TCR_p, contacts_TCR_MHC = filter_contacts(
                    contacts_df,
                    chains['tcra_chain'],
                    chains['tcrb_chain'],
                    chains['peptide_chain'],
                    chains['mhc_chain'],
                    remove_X=True)

        # Mapping with IMGT convention the similar TCR
        result_string += "\n Sequences reenumbered with IMGT convention\n"
        alpha_pdb, beta_pdb, epitope_pdb = extract_specific_sequences(pdb_file_path, chain_dict)
        
        anarci_a = run_anarci (alpha_pdb)
        anarci_b = run_anarci (beta_pdb)
        
        parsed_anarci_a = parse_anarci_output(anarci_a)
        parsed_anarci_b = parse_anarci_output(anarci_b)
        
        residues_a = extract_residues_and_resids(pdb_file_path, chains['tcra_chain'])
        residues_b = extract_residues_and_resids(pdb_file_path, chains['tcrb_chain'])
        
        mapping_a = map_imgt_to_original(parsed_anarci_a, residues_a)
        mapping_b = map_imgt_to_original(parsed_anarci_b, residues_b)

        imgt_mappings = {similar_tcr: {chains['tcra_chain']: mapping_a, chains['tcrb_chain']: mapping_b}}
        contacts_TCR_p = add_imgt_mappings(contacts_TCR_p, imgt_mappings)
        contacts_TCR_MHC = add_imgt_mappings(contacts_TCR_MHC, imgt_mappings)
        
        # Process input TCR
        result_string += f"\nMapping TCR {tcr_id} into {similar_tcr}\n"
        anarci_input_a = run_anarci (alpha_seq)
        anarci_input_b = run_anarci (beta_seq)
        
        parsed_input_a = parse_anarci_output(anarci_input_a)
        parsed_input_b = parse_anarci_output(anarci_input_b)
        
        imgt_dict_a = dict(parsed_input_a)
        imgt_dict_b = dict(parsed_input_b)

        contacts_TCR_p[tcr_id] = contacts_TCR_p.apply(lambda row: imgt_dict_a.get(row['imgt_from'], None) if row['chain_from'] == chains['tcra_chain']
                            else imgt_dict_b.get(row['imgt_from'], None) if row['chain_from'] == chains['tcrb_chain']
                            else None, axis=1)
        
        contacts_TCR_MHC[tcr_id] = contacts_TCR_MHC.apply(lambda row: imgt_dict_a.get(row['imgt_from'], None) if row['chain_from'] == chains['tcra_chain']
                            else imgt_dict_b.get(row['imgt_from'], None) if row['chain_from'] == chains['tcrb_chain']
                            else None, axis=1)

        ##### PROCESSING EPITOPE #####
        result_string += f"\n-> Processed Epitope: {epitope_seq}\n"
        contacts_TCR_p['epitope'] = contacts_TCR_p.apply(lambda row: map_epitope_residue(row, epitope_seq), axis=1)

        # Add TCR-Potential (looping for P1-P9)
        result_string += "\nCalculated TCR-peptide potential\n"
        total_scores = {}

        # Calculate potential scores and totals for each position
        try:
            if peptide_length > 13 or peptide_length < 8:
                raise ValueError(f"Unsupported peptide length: {peptide_length}. Max supported length is 13.")

            if peptide_length <= 10:
                for i in range(1, peptide_length+1):
                    potential_column = f'potential_P{i}'
                    # Apply the potential function based on the residue assignment
                    contacts_TCR_p[potential_column] = contacts_TCR_p.apply(
                        lambda row: get_potential(row, tcr_p_potential[f'tcr_p_potential_P{i}'], tcr_id, "epitope")
                        if row['resid_to'] == i else None, axis=1)
                    
                    # Sum the potential scores for each position
                    total_scores[potential_column] = tcr_p_potential[potential_column].sum(skipna=True)
                    total_contacts_position = contacts_TCR_p[potential_column].notnull().sum()

                    # Add the total score to the result string
                    result_string += f"Total score for P{i}: {total_scores[potential_column]}"
                    result_string += f"Total contacts for P{i}: {total_contacts_position}\n"
            else:
                assignments_dict = {11: {1: 1, 2: 2, 3: 3, 4: 4, 5: 4, 6: 5, 7: 5, 8: 8, 9: 5, 10: 9, 11: 10},
                                12: {1: 1, 2: 2, 3: 3, 4: 2, 5: 5, 6: 7, 7: 7, 8: 6, 9: 9, 10: 9, 11: 9, 12: 10},
                                13: {1: 1, 2: 2, 3: 3, 4: 6, 5: 3, 6: 4, 7: 8, 8: 6, 9: 8, 10: 9, 11: 9, 12: 9, 13: 10}}
                for i in range(1, peptide_length+1):
                    potential_column = f'potential_P{i}'
                    
                    real_index = assignments_dict[peptide_length][i]

                    contacts_TCR_p[potential_column] = contacts_TCR_p.apply(
                        lambda row: get_potential(row, tcr_p_potential[f'tcr_p_potential_P{real_index}'], tcr_id, "epitope")
                        if row['resid_to'] == i else None, axis=1)
                    
                    # Sum the potential scores for each position
                    total_scores[potential_column] = tcr_p_potential[potential_column].sum(skipna=True)
                    total_contacts_position = contacts_TCR_p[potential_column].notnull().sum()

                    # Add the total score to the result string
                    result_string += f"Total score for P{i}: {total_scores[potential_column]}"
                    result_string += f"Total contacts for P{i}: {total_contacts_position}\n"

        except Exception as e:
            print("Error calculating TCR-p potential")

        # Calculate total score for all positions
        total_contacts_tcr = len(tcr_p_potential)
        total_score_all = sum(total_scores.values())/total_contacts_tcr
        result_string += f"Total score for all positions normalized by the number of total contacts ({total_contacts_tcr}): {total_score_all}\n"

        #### PROCESSING MHC ####
        if mhc_allele:
            try:
                mhc_seq_match = mhc_df[mhc_df['mhc_allele'] == mhc_allele]
                if not mhc_seq_match.empty:
                    mhc_seq = mhc_seq_match['mhc_seq'].values[0]
                else:
                    print(f"Warning: No MHC sequence found for allele {mhc_allele}")
                    mhc_seq = None
                    return None
                
                result_string += f"\n-> Processing MHC-I allele: {mhc_allele}\n"
                seq_pdb = extract_sequences(pdb_file_path)
                aligned_seq_pdb, aligned_seq_query, score = global_alignment(seq_pdb[chains['mhc_chain']], mhc_seq)
                residues_M = extract_residues_and_resids(pdb_file_path, chains['mhc_chain']) 
                mapped_residues = map_alignment_to_residues(aligned_seq_pdb, aligned_seq_query, residues_M)
                    
                df_tuples = pd.DataFrame(mapped_residues, columns=['resid', 'mhc_pdb', mhc_allele])
                contacts_TCR_MHC = pd.merge(
                        contacts_TCR_MHC, 
                        df_tuples[['resid', mhc_allele]], 
                        left_on='resid_to', 
                        right_on='resid', 
                        how='left')
                    
                contacts_TCR_MHC = contacts_TCR_MHC.drop(columns=['resid'])

                # Add TCR-MHC potential
                result_string += "\nCalculating TCR-MHC-I potential\n"
                contacts_TCR_MHC['potential'] = contacts_TCR_MHC.apply(
                        lambda row: get_potential(row, tcr_mhc_potential, tcr_id, mhc_allele), axis=1)
                    
                contacts_TCR_MHC['potential'] = pd.to_numeric(contacts_TCR_MHC['potential'], errors='coerce')
                total_contacts_mhc=len(contacts_TCR_MHC)
                total_score_mhc = contacts_TCR_MHC['potential'].sum()/total_contacts_mhc
                result_string +=f"Total score for TCR-MHC {mhc_allele} normalized by the total number of contacts ({total_contacts_mhc}): {total_score_mhc}"

            except Exception as e:
                print(f"Error processing MHC sequence for allele {mhc_allele}: {e}")
                mhc_seq = None

        else:
            result_string += "\nNo MHC allele provided. Skipping MHC processing.\n"
            total_score_mhc = None

        # Append results to the results container
        results["tcr_id"]=tcr_id
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
    parser.add_argument("-i", "--input_file", type=str, required=True, help='Input CSV file path with TCRs to process.')
    parser.add_argument("-e", "--epitope_seq", type=str, required=True, help='Epitope sequence')
    parser.add_argument("-a", "--mhc_allele", type=str, required=False, help='Allele. Format: A*02:01')
    parser.add_argument("-tcrp", "--tcr_potential_folder", type=str, required=True, help="TCR-peptide potential.")
    parser.add_argument("-tcrm", "--mhc_potential", type=str, required=True, help='TCR-MHC potential.')
    parser.add_argument("-o", "--output_file", type=str, required=False, help='Output CSV file path for scores')
    parser.add_argument("-w", "--max_workers", type=int, required=False, default=8, help='Num of workers.')
    parser.add_argument("-v", "--verbose", action="store_false", default=True, help='Verbose mode')
    args = parser.parse_args()

    # Load data
    print("Loading data...")
    tcr_df = pd.read_csv(args.input_file)
    mhc_df = pd.read_csv(os.path.join (DATA_DIR,"all_mhc_seqs.csv"))
    epitope_seq = args.epitope_seq
    epitope_length = len(epitope_seq)
    mhc_allele = args.mhc_allele + ":01:01"
    chain_dict = parse_general_file("./structures_annotation/general.txt")

    # Load models acording to epitope_length
    print("Loading models...")
    print("    1) Loading TCR-p models")
    
    if epitope_length <= 10:
        epitope_length_folder = f"-L{epitope_length}"
    else:
        epitope_length_folder = "-L11"

    # Get the CSV files from the appropriate folder
    csv_files = [
        os.path.join(args.tcr_potential_folder, folder, file)
        for folder in os.listdir(args.tcr_potential_folder)
        if folder.startswith("TCR-p") and folder.endswith(epitope_length_folder)
        for file in os.listdir(os.path.join(args.tcr_potential_folder, folder))
        if file.endswith(".csv")]
    
    tcr_p_potential = {}
    if not csv_files:
        print(f"No CSV files found in the directory: {folder_path}")
    else:
        for file in csv_files:
            try:
                file_name = os.path.basename(file)
                potential_number = file_name.split('_P')[1].split('.csv')[0]
                potential_key = f'tcr_p_potential_P{potential_number}'
                tcr_p_potential[potential_key] = pd.read_csv(file)
            except Exception as e:
                print(f"Failed to load {file}: {e}")
    print("    2) Loading TCR-MHC model")         
    tcr_mhc_potential = pd.read_csv(args.mhc_potential)

    # Process each TCR in parallel
    results = []
    with ProcessPoolExecutor(max_workers=int(args.max_workers)) as executor:
        for tcr_id, alpha_seq, beta_seq in zip(tcr_df['tcr_id'], tcr_df['alpha_seq'], tcr_df['beta_seq']):
            try:
                result = executor.submit(process_tcr, tcr_id, alpha_seq, beta_seq, epitope_seq, tcr_p_potential, tcr_mhc_potential, mhc_df, chain_dict, args.verbose, mhc_allele)
                results.append(result)
            except Exception as e:
                print(f"Error processing TCR {tcr_id}: {e}")
                continue
                
    # Collect and save results
    results_list = [r.result() for r in results]
    results_df = pd.DataFrame(results_list)
    
    # Save the results as a CSV in the specified file path
    results_df.to_csv(args.scores_out, index=False)
    print(f"Results saved to {args.scores_out}")
    
if __name__ == "__main__":
    main()