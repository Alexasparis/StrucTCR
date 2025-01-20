#!/usr/bin/env python3
# Example of usage python run_seqs.py -m main.py -i sequences  -mfp TCR-p_all -mfm TCR_MHC_all.csv -w 8

import os
import subprocess
import pandas as pd
import argparse
import warnings
warnings.filterwarnings("ignore")

from mapping import *
from utils import *
from extract_contacts import *

def process_file(main, input_folder, input_file, output_folder, output_file_name, peptide_potential_dir, model_file_m, workers):    
    """Process a single input CSV file."""
    pdb_id = str(input_file.split('_')[0])
    output_file = os.path.join(output_folder, f'{output_file_name}_{pdb_id}.csv')  # Output file 

    # Check if the output file already exists to skip processing
    if os.path.exists(os.path.join("./output", output_file)):
        print(f"Output file {output_file} already exists. Skipping {pdb_id}.")
        return

    # Print progress message
    print(f"Processing PDB file: {pdb_id}")

    # Extract specific sequences from the PDB file
    try:
        seqs_info=pd.read_csv("./seqs_info/vdjdb_filtered.csv")
        pep_seq = seqs_info.loc[seqs_info['TCR_name'] == int(pdb_id)]['Epitope'].values[0]
        print(f"Peptide sequence: {pep_seq}")
    except Exception as e:
        print(f"Error extracting sequences from {pdb_file_path}: {e}")
        return  # Terminate if there's an error
    
    # Command to execute the script main_fast.py
    try:
        reference_allele = seqs_info.loc[seqs_info['TCR_name'] == int(pdb_id)]['MHC_allele'].values[0]
        reference_allele = f"{reference_allele}"
        reference_allele = reference_allele.replace("HLA-", "")
        print(f"Reference allele: {reference_allele}")

    except IndexError:
        reference_allele = None

    if reference_allele and isinstance(reference_allele, str):
        command = [
            'python3', f'{main}',
            '-t', f'./input/{input_folder}/{input_file}',  # Use the current input file
            '-e', f'{pep_seq}',  # Peptide sequence
            '-a', f'{reference_allele}',  # Correct formatting for allele
            '-tp', f"./model/{peptide_potential_dir}",  
            '-mp', f'./model/{model_file_m}',
            '-w', f'{workers}',
            '-s', os.path.join("./output/", output_file)
        ]
    else:
        command = [
            'python3', f'{main}',
            '-t', f'./input/{input_folder}/{input_file}',  # Use the current input file
            '-e', f'{pep_seq}',  # Peptide sequence
            '-tp', f"./model/{peptide_potential_dir}",  
            '-mp', f'./model/{model_file_m}',
            '-w', f'{workers}',
            '-s', os.path.join("./output/", output_file)
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

# Set up argument parsing
def main():
    parser = argparse.ArgumentParser(description="Process CSV files and extract information from PDB files.")
    parser.add_argument('-m', '--main', required=True, help="Input python script to execute.")
    parser.add_argument('-i', '--input_folder', required=True, help="Input folder containing the CSV files.")
    parser.add_argument('-mfp', '--peptide_potential_dir', required=True, help="Path to the first model file.")
    parser.add_argument('-mfm', '--model_file_m', required=True, help="Path to the second model file.")
    parser.add_argument('-w', '--workers', required=True, help="Workers")

    args = parser.parse_args()
    os.makedirs(f'./output/{args.input_folder}', exist_ok=True)

    # List to store the input files to be processed
    input_files = [f for f in os.listdir(f'./input/{args.input_folder}') if f.endswith('.csv')]
    print(f"Found {len(input_files)} input files to process.")

    # Process files sequentially using a for loop
    for input_file in input_files:
        print(f"Started processing {input_file}...")
        process_file(args.main, args.input_folder, input_file, args.input_folder, args.input_folder, args.peptide_potential_dir, args.model_file_m, args.workers)
        print(f"Finished processing {input_file}.")

    print("All files processed.")

if __name__ == "__main__":
    main()
