# Train_model

# This is the main function to train the model and extract the potential matrices from a given pdb_file database composed by TCR-pMHC complexes.
#Example of usage
#python3 -p ./pdb_nr [- d 5 -out_p ./model/TCRen_TCR_p.csv -out_m ./model/TCRen_TCR_MHC.csv -g ./structures_annotation/general.txt]

# Import libraries
import argparse
import os
import pandas as pd

from extract_contacts import *
from TCRen_calc import *

def folder(arg):
    """
    Custom type for argparse to handle folder paths.
    """
    if not os.path.isdir(arg):
        raise argparse.ArgumentTypeError(f"The folder {arg} does not exist.")
    return arg

def main():
    parser = argparse.ArgumentParser(description='Calculate TCRen potential from TCR-pMHC complexes database.')
    parser.add_argument("-p", "--pdbs", type=folder, required=True, help='Folder with PDB files to process.')
    parser.add_argument("-d", "--distance", type=float, default=5.0, help='Distance threshold for contacts. Default 5 Å.')
    parser.add_argument("-out_p", "--output_p", type=str, default='./model/TCRen_TCR_p.csv', help='Output CSV file name for TCR peptide potential.')
    parser.add_argument("-out_m", "--output_m", type=str, default='./model/TCRen_TCR_MHC.csv', help='Output CSV file name for TCR MHC potential.')
    parser.add_argument("-g", "--general", type=str, default='./structures_annotation/general.txt', help='Path to general.txt file.')

    args = parser.parse_args()
    
    # Ensure the output directory exists
    output_dir = 'model'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Read chain information from general.txt
    df = pd.read_csv(args.general, sep='\t')
    chain_dict = {}
    
    for pdb_id, group in df.groupby('pdb.id'):
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
    
    # Process PDB files
    pdb_files = [os.path.join(args.pdbs, f) for f in os.listdir(args.pdbs) if f.endswith('.pdb')]
    
    # Extract contacts from PDB files
    print("Extracting contacts from PDB files")
    contacts = extract_contacts(pdb_files, chain_dict, distance=args.distance)
    contacts.to_csv(os.path.join(output_dir, 'contacts_training.csv'), index=False)
    
    # Filter contacts
    print("Filtering contacts")
    
    #Initialize DataFrames to store all contacts
    all_contacts_TCR_p = pd.DataFrame()
    all_contacts_TCR_MHC = pd.DataFrame()
    
    # Use the chain_dict to filter contacts for each pdb_id
    for pdb_id, chains in chain_dict.items():
        if all(chains.values()):  # Ensure all chain identifiers are present
            contacts_TCR_p, contacts_TCR_MHC = filter_contacts(
                contacts[contacts['pdb_id'] == pdb_id],
                chains['tcra_chain'],
                chains['tcrb_chain'],
                chains['peptide_chain'],
                chains['mhc_chain'],
                threshold=2)
        
        # Accumulate the filtered contacts for each PDB
            all_contacts_TCR_p = pd.concat([all_contacts_TCR_p, contacts_TCR_p], ignore_index=True)
            all_contacts_TCR_MHC = pd.concat([all_contacts_TCR_MHC, contacts_TCR_MHC], ignore_index=True)
        else:
            print(f"Missing chain information for PDB ID: {pdb_id}. Skipping...")
    
    # Calculate TCRen potentials for all accumulated contacts
    print("Calculating TCRen potentials for all PDB IDs")
    print("Calculating TCR-peptide potential")
    data_TCR_p = calculate_TCRen(all_contacts_TCR_p, peptide=True)
    print("Calculating TCR-MHCI potential")
    data_TCR_MHC = calculate_TCRen(all_contacts_TCR_MHC, peptide=False)
    
    # Save the results to the specified output files
    data_TCR_p.to_csv(args.output_p, index=False)
    data_TCR_MHC.to_csv(args.output_m, index=False)
    
    print(f"TCRen potential for TCR-peptide saved to {args.output_p}")
    print(f"TCRen potential for TCR-MHC saved to {args.output_m}")
    print("Processing complete.")

if __name__ == "__main__":
    main()
