#!/usr/bin/env python3

# This script extracts contact maps from .pdb files of an input folder.
# The output is a folder with extracted contacts between TCRa-P, TCRb-P, TCRa-MHC, TCRb-MHC.

# Example of execution: python contact_maps_pdb.py -pdb pdb_files -out contact_maps -workers 7


import os
import concurrent.futures
import argparse
import pandas as pd
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))
from extract_contacts import extract_contacts
from utils import parse_general_file
from config import STRUCTURES_ANNOTATION_DIR

def validate_chain_columns(csv_file, chain_dict):
    """
    Verifies that the 'chain_from' and 'chain_to' columns in the CSV file contain exactly 2 unique strings each,
    and that these chains match the expected types (TCRα, TCRβ, peptide, MHC).
    """
    df = pd.read_csv(csv_file)
    if 'chain_from' not in df.columns or 'chain_to' not in df.columns:
        raise ValueError(f"The file {csv_file} does not contain 'chain_from' or 'chain_to' columns.")
    
    pdb_id = os.path.basename(csv_file).split('_')[0]
    
    if pdb_id not in chain_dict:
        chain_dict = {f'{pdb_id}': {'tcra_chain': 'D',
                                    'tcrb_chain': 'E',
                                    'peptide_chain': 'C',
                                    'b2m_chain': 'B',
                                    'mhc_chain': 'A'}}
        print(f"Chain info not found for {pdb_id}. Using default chains...")
    
    expected_chains = chain_dict.get(pdb_id)
    expected_chain_from = [expected_chains['tcra_chain'], expected_chains['tcrb_chain']]
    expected_chain_to = [expected_chains['peptide_chain'], expected_chains['mhc_chain']]
    
    unique_chain_from = df['chain_from'].unique()
    unique_chain_to = df['chain_to'].unique()
    
    is_valid = (len(unique_chain_from) == 2 and len(unique_chain_to) == 2 and
                all(chain in expected_chain_from for chain in unique_chain_from) and
                all(chain in expected_chain_to for chain in unique_chain_to))
    
    return is_valid

def process_pdb_file(pdb_file, pdb_dir, output_dir, chain_dict):
    if not pdb_file.endswith('.pdb'):
        return  # Skip non-PDB files

    pdb_id = os.path.basename(pdb_file).split('_')[0]
    model_number = os.path.basename(pdb_file).split('_')[1].split('.')[0]
    pdb_path = os.path.join(pdb_dir, pdb_file)
    print(pdb_id, model_number)
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f'{pdb_id}_{model_number}_contacts.csv')

    if os.path.exists(output_file):
        print(f"File {output_file} exists, omitting...")
        return

    if pdb_id in chain_dict:
        contacts_df = extract_contacts([pdb_path], chain_dict)
    else:
        print(f"Chain info not found for {pdb_id}. Using default chains...")
        chain_dict_local = {f'{pdb_id}_{model_number}': {
            'tcra_chain': 'D', 'tcrb_chain': 'E',
            'peptide_chain': 'C', 'b2m_chain': 'B',
            'mhc_chain': 'A'
        }}
        contacts_df = extract_contacts([pdb_path], chain_dict_local)
    
    temp_file = f"{output_file}.temp"
    contacts_df.to_csv(temp_file, index=False)
    
    if validate_chain_columns(temp_file, chain_dict):
        os.rename(temp_file, output_file)
        print(f"Saved valid contacts in {output_file}.")
    else:
        os.remove(temp_file)
        print(f"File {output_file} is not valid and has been discarded.")

def main():
    parser = argparse.ArgumentParser(description="Extract valid contact maps from PDB files.")
    parser.add_argument("-pdb", "--pdb_folder", required=True, help="Path to the folder containing PDB files.")
    parser.add_argument("-out", "--output_folder", required=True, help="Path to the output folder for contact maps.")
    parser.add_argument("-workers", "--num_workers", type=int, default=os.cpu_count() - 1,
                        help="Number of worker processes (default: max cores - 1).")
    args = parser.parse_args()
    path_general = os.path.join(STRUCTURES_ANNOTATION_DIR, "general.txt")
    chain_dict = parse_general_file(path_general)
    pdb_files = [f for f in os.listdir(args.pdb_folder) if f.endswith('.pdb')]

    with concurrent.futures.ProcessPoolExecutor(max_workers=args.num_workers) as executor:
        executor.map(process_pdb_file, pdb_files, [args.pdb_folder] * len(pdb_files), 
                     [args.output_folder] * len(pdb_files), [chain_dict] * len(pdb_files))

if __name__ == "__main__":
    main()

