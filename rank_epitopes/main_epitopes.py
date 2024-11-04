# Main function to rank TCRs given an input neoantigenic sequence.

# Import libraries
import argparse
import pandas as pd
import os
import warnings  
import ast

from find_contact_map import *
from mapping import *
from extract_contacts import *
from TCRen_calc import get_TCRen
from select_nr_set import parse_general_file, extract_specific_sequences

warnings.filterwarnings("ignore")

def folder(arg):
    """
    Custom type for argparse to handle folder paths.
    """
    if not os.path.isdir(arg):
        raise argparse.ArgumentTypeError(f"The folder {arg} does not exist.")
    return arg

def main():
    parser = argparse.ArgumentParser(description='Calculate TCRen potential from TCR-pMHC complexes database.')
    parser.add_argument("-t", "--tcr_file", type=str, required=True, help='CSV file with TCR to process. Format: tcr_id,alpha_seq,beta_seq')
    parser.add_argument("-e", "--epitope_file", type=str, required=True, help='CSV file with epitopes to process.')
    parser.add_argument("-pp", "--peptide_potential", type=str, required=True, help='CSV file with TCR-peptide TCRen potential. Format: residue_from,residue_to,TCRen')
    parser.add_argument("-g", "--general_file", type=str, required=True, default='./structures_annotation/general.txt', help='Annotation file with chain information for each pdb_file of pdb_folder')
    parser.add_argument("-p", "--pdb_folder", type=folder, required = False, default='./pdb_files', help='Database of structures, contact maps avaliable')
    parser.add_argument("-metric", "--similarity_metric", type=str, required= False, choices=['TCRdist', 'Levenshtein'], help='Metric to decide which complex to use: TCRdist or Levenshtein')
    parser.add_argument("-s", "--scores_file", type=str, default='scores.csv', help='Output CSV file for scores (default: scores.csv)')
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    output_dir = "output"
    os.makedirs(output_dir, exist_ok=True)
    
    # Load data
    epitopes_df = pd.read_csv(args.epitope_file)
    tcr_p_potential = pd.read_csv(args.peptide_potential)
    general_df = pd.read_csv(args.general_file, sep='\t')
    tcr_df = pd.read_csv(args.tcr_file)

    epitopes = []
    scores_tcr_p = []

    # Chain dictionaries
    chain_dict = {}
    for pdb_id, group in general_df.groupby('pdb.id'):
        chains = {
            'tcra_chain': None,
            'tcrb_chain': None,
            'peptide_chain': None,
            'mhc_chain': None
        }
        
        for _, row in group.iterrows():
            if row['chain.component'] == 'TCR' and row['chain.type'] == 'TRA':
                chains['tcra_chain'] = row['chain.id']
            elif row['chain.component'] == 'TCR' and row['chain.type'] == 'TRB':
                chains['tcrb_chain'] = row['chain.id']
            elif row['chain.component'] == 'PEPTIDE':
                chains['peptide_chain'] = row['chain.id']
            elif row['chain.component'] == 'MHC' and row['chain.supertype'] == 'MHCI' and row['chain.type'] == 'MHCa':
                chains['mhc_chain'] = row['chain.id']
        
        chain_dict[pdb_id] = chains  
        #Format chain_dict = {"pdb_id": {"tcra_chain": "chain_id", "tcrb_chain": "B", "peptide_chain": "C", "mhc_chain": "D"}, "pdb_id": {"tcra_chain": "E"...
    
    seq_dict = parse_general_file(args.general_file)
    
    ##### PROCESSING INPUT TCR #####
    print(f"\n #### Processing TCR:")
    
    for epitope in epitopes_df['peptide']:
        epitope=epitope

        alpha_seq = tcr_df['alpha_seq'].values[0]
        beta_seq = tcr_df['beta_seq'].values[0]
        tcr_id= tcr_df['tcr_id'].values[0]
        

        if args.similarity_metric == "TCRdist":
            similarity_df = pd.read_csv("./structures_annotation/closest_tcr_all.csv")
            similar_tcr = similarity_df[similarity_df['tcr_name'] == str(tcr_id)]['closest_tcr'].values[0]

            if isinstance(similar_tcr, str):
                if similar_tcr.startswith("[") and similar_tcr.endswith("]"):
                    pdb_id_similar = ast.literal_eval(similar_tcr)
                else:
                    pdb_id_similar = [similar_tcr]
            else:
                pdb_id_similar = [similar_tcr]

            if len(pdb_id_similar) > 1:
                distances = []
                for pdb in pdb_id_similar:
                    pdb_file_path = os.path.join("./pdb_files", f"{pdb}.pdb")
                    if os.path.isfile(pdb_file_path):
                        a_seq, b_seq, e_seq = extract_specific_sequences(pdb_file_path, seq_dict)
                        distance = calculate_sequence_distance(epitope, e_seq)
                        distances.append(distance)
                    else:
                        print(f"File not found: {pdb_file_path}")

                if distances:
                    closest_pdb_index = distances.index(min(distances))
                    pdb_id_similar = pdb_id_similar[closest_pdb_index]
                    print(f"Match found for {tcr_id}: {pdb_id_similar}.")

            elif len(pdb_id_similar) == 1:
                pdb_id_similar = pdb_id_similar[0]
                print(f"Match found for {tcr_id}: {pdb_id_similar}.")
            else:
                print(f"No similar TCRs found for {tcr_id}. Try with Levenshtein distance.")
        
        elif not args.similarity_metric:      
            pdb_id_similar = tcr_id
            print(f"No similarity metric provided. Using the structural data {pdb_id_similar} .") 
        
        ##### EXTRACT CONTACTS FROM PDB FILE #####
        
        if pdb_id_similar:
            pdb_file_path = os.path.join(args.pdb_folder, f"{pdb_id_similar}.pdb")
            if not os.path.isfile(pdb_file_path):
                print(f"Error: The file does not exist at {pdb_file_path}")
            else:
                contacts_df = pd.read_csv(f"./contact_maps/{pdb_id_similar}_contacts.csv")
                chains = chain_dict.get(pdb_id_similar, {})
                if all(chains.values()):  # Ensure all chain identifiers are present
                    contacts_TCR_p, contacts_TCR_MHC = filter_contacts(
                            contacts_df,
                            chains['tcra_chain'],
                            chains['tcrb_chain'],
                            chains['peptide_chain'],
                            chains['mhc_chain'],
                            threshold=1)
            
            # Mapping with IMGT numbering
            print("Renumbering sequences with IMGT convention")
            alpha_pdb, beta_pdb, epitope_pdb = extract_specific_sequences(pdb_file_path, seq_dict)

            anarci_D = run_anarci(alpha_pdb, 'alpha')
            anarci_E = run_anarci(beta_pdb, 'beta')
            parsed_anarci_D = parse_anarci_output(anarci_D)
            parsed_anarci_E = parse_anarci_output(anarci_E)

            residues_D = extract_residues_and_resids(pdb_file_path, chains['tcra_chain'])
            residues_E = extract_residues_and_resids(pdb_file_path, chains['tcrb_chain'])

            mapping_D = map_imgt_to_original(parsed_anarci_D, residues_D)
            mapping_E = map_imgt_to_original(parsed_anarci_E, residues_E)

            imgt_mappings = {pdb_id_similar: {chains['tcra_chain']: mapping_D, chains['tcrb_chain']: mapping_E}}
            
            contacts_imgt_M = add_imgt_mappings(contacts_TCR_MHC, imgt_mappings)
            contacts_imgt_P = add_imgt_mappings(contacts_TCR_p, imgt_mappings)
            
            print("IMGT mapping done")

        else:
            print("No matching pdb_file found")

        ##### PROCESS INPUT TCR #####
        
        print(f"Parsing TCR {tcr_id}")
        tcr_alpha = run_anarci(alpha_seq, "alpha")
        tcr_beta = run_anarci(beta_seq, "beta")
        parsed_tcr_D = parse_anarci_output(tcr_alpha)
        parsed_tcr_E = parse_anarci_output(tcr_beta)
        imgt_dict_D = dict(parsed_tcr_D)
        imgt_dict_E = dict(parsed_tcr_E)

        contacts_TCR_p[tcr_id] = contacts_TCR_p.apply(lambda row: imgt_dict_D.get(row['imgt_from'], None) if row['chain_from'] == chains['tcra_chain']
                        else imgt_dict_E.get(row['imgt_from'], None) if row['chain_from'] == chains['tcrb_chain']
                        else None, axis=1)
        
        ##### PROCESS INPUT EPITOPES ####
        print(f"Processing epitope: {epitope}")
        contacts_TCR_p['epitope'] = contacts_TCR_p.apply(lambda row: map_epitope_residue(row, epitope), axis=1)
        print("Residues mapped")
        print("\nContact map TCR peptide \n", contacts_TCR_p.head()) # TO DEBUG
            
        # Add TCR-P TCRen potential
        print("Calculating TCRen potential")
        contacts_TCR_p['TCRen'] = contacts_TCR_p.apply(lambda row: get_TCRen(row, tcr_p_potential, tcr_id ,"epitope"), axis=1)
        contacts_TCR_p['TCRen'] = pd.to_numeric(contacts_TCR_p['TCRen'], errors='coerce')
        total_TCRen_p = contacts_TCR_p['TCRen'].sum()
        print(f"     -> TCRen score for TCR:{tcr_id} - Epitope {epitope}: {total_TCRen_p}")
            #print("Saving to csv")
            #contacts_TCR_p.to_csv(os.path.join(output_dir, f"{pdb_id_similar}_{tcr_id}_TCRen_p.csv"), index=False)
            #print(f"Mapped contacts saved to {os.path.join(output_dir, f'{pdb_id_similar}_{tcr_id}_TCRen_p.csv')}") 
        epitopes.append(epitope)
        scores_tcr_p.append(total_TCRen_p)  
   
        # Scores
    results_df = pd.DataFrame({
        'epitope': epitopes,
        'score_tcr_p': scores_tcr_p})
    
    final_csv_path = os.path.join(output_dir, args.scores_file)
    results_df.to_csv(final_csv_path, index=False)
    print(f"Final results saved to {final_csv_path}")

if __name__ == "__main__":
    main()
