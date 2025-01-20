#!/usr/bin/env python3
# Example of usage ./run.py -m main.py -i control_fold_1 -ocontrol_fold_1 -ofn scores_fold1 -mfp TCR-p_all -mfm TCR_MHC_all.csv -w 8

import os
import subprocess
import pandas as pd
import argparse
import warnings
warnings.filterwarnings("ignore")

from mapping import *
from utils import *
from extract_contacts import *

def process_file(main, input_folder, input_file, output_folder, output_file_name, peptide_potential_dir, model_file_m, seq_dict, workers):    
    """Process a single input CSV file."""
    pdb_id = input_file.split('_')[0]
    pdb_file_path = f'./pdb_files/{pdb_id}.pdb'  # Path to the PDB file
    output_file = os.path.join(output_folder, f'{output_file_name}_{pdb_id}.csv')  # Output file 

    # Check if the output file already exists to skip processing
    if os.path.exists(os.path.join("./output", output_file)):
        print(f"Output file {output_file} already exists. Skipping {pdb_id}.")
        return

    # Print progress message
    print(f"Processing PDB file: {pdb_file_path}")

    # Extract specific sequences from the PDB file
    try:
        tcra_seq, tcrb_seq, pep_seq = extract_specific_sequences(pdb_file_path, seq_dict, extract_sequences)
        print('Sequences extracted for:', pdb_id)
    except Exception as e:
        print(f"Error extracting sequences from {pdb_file_path}: {e}")
        return  # Terminate if there's an error

    # Command to execute the script main_fast.py
    alleles_df = pd.read_csv('./structures_annotation/mhc_annotation.csv')
    
    try:
        reference_allele = alleles_df.loc[alleles_df['pdb_id'] == pdb_id]['mhci_allele'].values[0]
        reference_allele = f"{reference_allele}"
    except IndexError:
        reference_allele = None

    if reference_allele and isinstance(reference_allele, str):
        command = [
            'python3', f'{main}',
            '-t', f'./input/{input_folder}/{input_file}',  # Use the current input file
            '-e', f'{pep_seq}',  # Peptide sequence
            '-a', f'{reference_allele}',  # Correct formatting for allele
            '-tp', f"{peptide_potential_dir}",  
            '-mp', f'{model_file_m}',
            '-w', f'{workers}',
            '-s', output_file
        ]
    else:
        command = [
            'python3', f'{main}',
            '-t', f'./input/{input_folder}/{input_file}',  # Use the current input file
            '-e', f'{pep_seq}',  # Peptide sequence
            '-tp', f"./model/{peptide_potential_dir}",  
            '-mp', f'./model/{model_file_m}',
            '-w', f'{workers}',
            '-s', output_file
        ]

    # Print progress message before executing the command
    print(f"Running command: {' '.join(command)}")

    # Use subprocess.run 
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    print(result.stdout)
    if result.stderr:
        print(f"Error: {result.stderr}")

    if result.returncode == 0:
        print("El comando se ejecutó correctamente")
    else:
        print(f"Hubo un error en la ejecución del comando. Código de error: {result.returncode}")

    # Print the output and errors if any
    print(result.stdout)
    if result.stderr:
        print(f"Error: {result.stderr}")

    # Check if the process finished successfully
    if result.returncode == 0:
        print(f"Successfully processed {pdb_id} and generated {output_file}")
    else:
        print(f"Error occurred while processing {pdb_id}")


def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Process CSV files and extract information from PDB files.")
    parser.add_argument('-m', '--main', required=True, help="Input python script to execute.")
    parser.add_argument('-i', '--input_folder', required=True, help="Input folder containing the CSV files.")
    parser.add_argument('-o', '--output_folder', required=True, help="Output folder to save processed files.")
    parser.add_argument('-ofn', '--output_file_name', required=True, help="Base name for the output files.")
    parser.add_argument('-mfp', '--peptide_potential_dir', required=True, help="Path to the first model file.")
    parser.add_argument('-mfm', '--model_file_m', required=True, help="Path to the second model file.")
    parser.add_argument('-w', '--workers', required=True, help="Workers")

    args = parser.parse_args()
    os.makedirs(f'./output/{args.output_folder}', exist_ok=True)
    
    # Initialize seq_dict from the general.txt file
    seq_dict = parse_general_file('./structures_annotation/general.txt')
    print("Initialized sequence dictionary from general.txt.")

    # List to store the input files to be processed
    input_files = [f for f in os.listdir(f'./input/{args.input_folder}') if f.endswith('.csv')]
    print(f"Found {len(input_files)} input files to process.")

    # Process files sequentially using a for loop
    for input_file in input_files:
        print(f"Started processing {input_file}...")
        process_file(args.main, args.input_folder, input_file, args.output_folder, args.output_file_name, args.peptide_potential_dir, args.model_file_m, seq_dict, args.workers)
        print(f"Finished processing {input_file}.")

    print("All files processed.")

if __name__ == "__main__":
    main()
