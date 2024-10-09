import os
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
import pandas as pd

from find_contact_map import *
from mapping import *
from select_nr_set import *
from extract_contacts import *

def process_file(input_file, seq_dict):
    """Process a single input CSV file."""
    pdb_id = input_file.split('_')[0]  # Extract the PDB ID from the filename
    pdb_file_path = f'./pdb_files/{pdb_id}.pdb'  # Path to the PDB file
    output_file = f'./training_set/scores_training_all_{pdb_id}.csv'  # Output file

    # Check if the output file already exists to skip processing
    if os.path.exists(os.path.join("./output", output_file)):
        print(f"Output file {output_file} already exists. Skipping {pdb_id}.")
        return

    # Print progress message
    print(f"Processing PDB file: {pdb_file_path}")

    # Extract specific sequences from the PDB file
    try:
        tcra_seq, tcrb_seq, pep_seq = extract_specific_sequences(pdb_file_path, seq_dict)
        print('Sequences extracted for:', pdb_id)
    except Exception as e:
        print(f"Error extracting sequences from {pdb_file_path}: {e}")
        return  # Terminate if there's an error

    # Command to execute the script main_fast.py
    alleles_df = pd.read_csv('./structures_annotation/mhc_alleles.csv')
    
    try:
        reference_allele = alleles_df.loc[alleles_df['pdb_id'] == pdb_id]['mhci_allele'].values[0] + ":01:01"
    except IndexError:
        reference_allele = None

    if reference_allele and isinstance(reference_allele, str):
        command = [
            'python3', 'main_fast.py',
            '-p', './pdb_files/',
            '-t', f'./input/input_training/{input_file}',  # Use the current input file
            '-e', f'{pep_seq}',  # Peptide sequence
            '-m', './input/input_MHCs.csv',
            '-a', f'{reference_allele}',  # Correct formatting for allele
            '-pp', f'./model/TCRen_TCR_p_all.csv',
            '-mhcp', f'./model/TCRen_TCR_MHC_all.csv',
            '-g', './structures_annotation/general.txt',
            '-metric', 'TCRdist',
            '-s', output_file  # Output file
        ]
    else:
        command = [
            'python3', 'main_fast.py',
            '-p', './pdb_files/',
            '-t', f'./input/input_training/{input_file}',  # Use the current input file
            '-e', f'{pep_seq}',  # Peptide sequence
            '-m', './input/input_MHCs.csv',
            '-pp', f'./model/TCRen_TCR_p_all.csv',
            '-mhcp', f'./model/TCRen_TCR_MHC_all.csv',
            '-g', './structures_annotation/general.txt',
            '-metric', 'TCRdist',
            '-s', output_file  # Output file
        ]

    # Print progress message before executing the command
    print(f"Running command: {' '.join(command)}")

    # Execute the command
    try:
        subprocess.run(command, check=True)
        print(f"Successfully processed {pdb_id} and generated {output_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while processing {pdb_id}: {e}")

def main():
    # Initialize seq_dict from the general.txt file
    seq_dict = parse_general_file('./structures_annotation/general.txt')
    print("Initialized sequence dictionary from general.txt.")

    # List to store the input files to be processed
    input_files = [f for f in os.listdir('./input/input_training') if f.endswith('.csv')]
    print(f"Found {len(input_files)} input files to process.")

    # Use ProcessPoolExecutor to run processes in parallel
    with ProcessPoolExecutor(max_workers=8) as executor:
        # Submit jobs to be processed
        futures = {executor.submit(process_file, input_file, seq_dict): input_file for input_file in input_files}
        print("Started processing input files in parallel.")

        # Wait for all processes to finish
        for future in as_completed(futures):
            input_file = futures[future]
            try:
                future.result()  # Capture any exceptions that occur
                print(f"Successfully processed {input_file}.")
            except Exception as e:
                print(f"Exception occurred while processing {input_file}: {e}")

    print("All files processed.")

if __name__ == "__main__":
    main()