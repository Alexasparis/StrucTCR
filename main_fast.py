# Main function to rank TCRs given an input neoantigenic sequence.

# Example of usage
# python3 main.py -p ./pdb_files/ -t ./input/input_TCRs.csv -e HPVGEADYFE -m ./input/input_MHCs.csv -pp ./model/TCRen_TCR_p.csv -mhcp ./model/TCRen_TCR_MHC.csv -g ./structures_annotation/general.txt -metric TCRdist -s scores.csv

# Import libraries
import argparse
import pandas as pd
import os
import warnings

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
    parser.add_argument("-p", "--pdb_folder", type=folder, default='./pdb_files', help='Folder with the PDB files used to train the model')
    parser.add_argument("-t", "--tcr_file", type=str, required=True, help='CSV file with TCRs to process. Format: tcr_id,alpha_seq,beta_seq')
    parser.add_argument("-e", "--epitope_seq", type=str, required=True, help='Epitope sequence')
    parser.add_argument("-m", "--mhc_file", type=str, required=True, help='CSV file with MHCs to process. Format: mhc_allele,mhc_seq')
    parser.add_argument("-pp", "--peptide_potential", type=str, default='./model/TCRen_TCR_p.csv', help='CSV file with TCR-peptide TCRen potential. Format: residue_from,residue_to,TCRen')
    parser.add_argument("-mhcp", "--mhc_potential", type=str, default='./model/TCRen_TCR_MHC.csv', help='CSV file with TCR-MHC TCRen potential. Format: residue_from,residue_to,TCRen')
    parser.add_argument("-g", "--general_file", type=str, default='./structures_annotation/general.txt', help='Annotation file with chain information for each pdb_file of pdb_folder')
    parser.add_argument("-metric", "--similarity_metric", type=str, default='TCRdist', choices=['TCRdist', 'Levenshtein'], help='Metric to decide which complex to use: TCRdist or Levenshtein')
    parser.add_argument("-s", "--scores_file", type=str, default='scores.csv', help='Output CSV file for scores (default: scores.csv)')

    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    output_dir = "output"
    os.makedirs(output_dir, exist_ok=True)
    
    # Load data
    tcr_df = pd.read_csv(args.tcr_file)
    mhc_df = pd.read_csv(args.mhc_file)
    epitope_seq = args.epitope_seq
    tcr_p_potential = pd.read_csv(args.peptide_potential)
    tcr_mhc_potential = pd.read_csv(args.mhc_potential)
    general_df = pd.read_csv(args.general_file, sep='\t')
    
    #Results
    tcr_ids = []
    mhc_alleles = []
    scores_tcr_p = []
    scores_tcr_mhc = []
    
    # Chain dictionary
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
    
    ##### PROCESSING INPUT TCRs #####
    
    print("\n#### Processing TCRs ####")
    
    for index, row in tcr_df.iterrows():
        
        #Get TCR data from pdb files
        tcr_id = row['tcr_id']
        print(f"\nProcessing TCR: {tcr_id}")
        alpha_seq = row['alpha_seq']
        beta_seq = row['beta_seq']
        pdb_id_similar = None
        
        # Compute similarity with TCRs form PDB files
        if args.similarity_metric == "TCRdist":
            v_gene_a, j_gene_a = get_germlines(alpha_seq)
            v_gene_b, j_gene_b = get_germlines(beta_seq)
            values = [v_gene_a, j_gene_a, v_gene_b, j_gene_b] 
            db_tcrdist = pd.read_csv("./TCRdist/alphabeta_gammadelta_db.tsv", sep="\t")
            result = db_tcrdist[db_tcrdist['id'].isin(values)]
            all_present = all(result['id'].isin(values))
            
            if all_present:
                tcr_dist_df = pd.read_csv("./structures_annotation/TCRdist_df.csv")
                pdb_id_similar = find_closest_tcr(tcr_dist_df, alpha_seq, beta_seq) # Parsed df excluding alleles not present in the db of the method
            else:
                print("Not possible to calculate TCRdist because alleles of query TCR are not in the db")
                print("Calculating similarities based on Levenshtein distances")
                df_distances = create_distance_matrix(df_lev, alpha_seq, beta_seq, epitope_seq)
                pdb_id_similar = get_min_combined_distance(df_distances)
                print(f"pdb_id_similar: {pdb_id_similar}")
                
        elif args.similarity_metric == "Levenshtein":
            #df_lev = pd.DataFrame(cdrs_results)
            df_lev=pd.read_csv("./structures_annotation/Levenshtein_df.csv")
            df_distances = create_distance_matrix(df_lev, alpha_seq, beta_seq, epitope_seq)
            pdb_id_similar = get_min_combined_distance(df_distances)
            print(f"pdb_id_similar: {pdb_id_similar}")
            
    ##### EXTRACT CONTACTS FROM PDB FILE #####
    
        if pdb_id_similar:
            print(f"Matching contact map found in {pdb_id_similar}")
            pdb_file_path = os.path.join(args.pdb_folder, f"{pdb_id_similar}.pdb")

            if pdb_id_similar not in processed_pdb_ids:
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
                            threshold=2
                        )

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

            contacts_TCR_p.to_csv(os.path.join(output_dir, f"{pdb_id_similar}_{tcr_id}_contacts_peptide.csv"), index=False)
            contacts_TCR_MHC.to_csv(os.path.join(output_dir, f"{pdb_id_similar}_{tcr_id}_contacts_mhc.csv"), index=False)
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
        contacts_TCR_MHC[tcr_id] = contacts_TCR_MHC.apply(lambda row: imgt_dict_D.get(row['imgt_from'], None) if row['chain_from'] == chains['tcra_chain']
                        else imgt_dict_E.get(row['imgt_from'], None) if row['chain_from'] == chains['tcrb_chain']
                        else None, axis=1)

        ##### PROCESS INPUT EPITOPE ####
        
        print(f"Processing epitope: {epitope_seq}")
        contacts_TCR_p['epitope'] = contacts_TCR_p.apply(lambda row: map_epitope_residue(row, epitope_seq), axis=1)
        print("Residues mapped")
        print("\nContact map TCR peptide \n", contacts_TCR_p.head())
        
	# Add TCR-P TCRen potential
        print("Calculating TCRen potential")
        contacts_TCR_p['TCRen'] = contacts_TCR_p.apply(lambda row: get_TCRen(row, tcr_p_potential, tcr_id ,"epitope"), axis=1)
        contacts_TCR_p['TCRen'] = pd.to_numeric(contacts_TCR_p['TCRen'], errors='coerce')
        total_TCRen_p = contacts_TCR_p['TCRen'].sum()
        print(f"TCRen score for TCR-peptide: {total_TCRen_p}")
        print("Saving to csv")
        contacts_TCR_p.to_csv(os.path.join(output_dir, f"{pdb_id_similar}_{tcr_id}_TCRen_p.csv"), index=False)
        print(f"Mapped contacts saved to {os.path.join(output_dir, f'{pdb_id_similar}_{tcr_id}_TCRen_p.csv')}") 
        
        ##### PROCESS INPUT MHCs #####
        
        for index, row in mhc_df.iterrows():
            seq_pdb = extract_sequences(pdb_file_path)
            mhc_allele = row['mhc_allele']
            mhc_seq = row['mhc_seq']    
            print(f"\nProcessing MHC allele: {mhc_allele}")
            aligned_seq_pdb, aligned_seq_query, score = global_alignment(seq_pdb[chains['mhc_chain']], mhc_seq)
            residues_M = extract_residues_and_resids(pdb_file_path, chains['mhc_chain']) 
            mapped_residues = map_alignment_to_residues(aligned_seq_pdb, aligned_seq_query, residues_M)
            df_tuples = pd.DataFrame(mapped_residues, columns=['resid', 'mhc_pdb', mhc_allele])
            contacts_TCR_MHC_updated = pd.merge(contacts_TCR_MHC, df_tuples[['resid', mhc_allele]], left_on='resid_to', right_on='resid', how='left')
            contacts_TCR_MHC_updated = contacts_TCR_MHC_updated.drop(columns=['resid'])
            print("Residues mapped")
            print("\nContact map TCR MHC \n", contacts_TCR_MHC_updated.head())
            
            # Add TCR-MHC TCRen potential
            print("Calculating TCRen potential")
            contacts_TCR_MHC_updated['TCRen'] = contacts_TCR_MHC_updated.apply(lambda row: get_TCRen(row, tcr_mhc_potential, tcr_id, mhc_allele), axis=1)
            contacts_TCR_MHC_updated['TCRen'] = pd.to_numeric(contacts_TCR_MHC_updated['TCRen'], errors='coerce')
            total_TCRen_mhc = contacts_TCR_MHC_updated['TCRen'].sum()
            print(f"TCRen score for TCR-MHC: {total_TCRen_mhc}")
            
            tcr_ids.append(tcr_id)
            mhc_alleles.append(mhc_allele)
            scores_tcr_p.append(total_TCRen_p)  
            scores_tcr_mhc.append(total_TCRen_mhc) 
            
            print("Saving to csv")
            contacts_TCR_MHC_updated.to_csv(os.path.join(output_dir, f"{pdb_id}_{tcr_id}_{mhc_allele}_TCRen_mhc.csv"), index=False)
            print(f"Mapped contacts saved to {os.path.join(output_dir, f'{pdb_id}_{tcr_id}_{mhc_allele}_TCRen_mhc.csv')}") 
            print("\n")
   
    # Scores
    results_df = pd.DataFrame({
    'tcr_id': tcr_ids,
    'mhc_allele': mhc_alleles,
    'score_tcr_p': scores_tcr_p,
    'score_tcr_MHC': scores_tcr_mhc})
    
    final_csv_path = os.path.join(output_dir, args.scores_file)
    results_df.to_csv(final_csv_path, index=False)
    print(f"Final results saved to {final_csv_path}")

if __name__ == "__main__":
    main()
