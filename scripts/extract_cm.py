import os
import concurrent.futures
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))
from extract_contacts import extract_contacts
from utils import parse_general_file
chain_dict = parse_general_file ("../data/structures_annotation/general.txt")

# Function to process a single PDB file
def process_pdb_file(pdb_file):
    if not pdb_file.endswith('.pdb'):
        return  # Skip non-PDB files

    pdb_id = os.path.basename(pdb_file).split('_')[0]  # Extract PDB ID
    model_number = os.path.basename(pdb_file).split('_')[1].split('.')[0]  
    pdb_path = f'../../one_model/{pdb_file}'  # Full path

    # Create output directory if not exists
    output_dir = '../../one_model/contact_maps'
    os.makedirs(output_dir, exist_ok=True)

    # Define the output file path
    output_file = f'{output_dir}/{pdb_id}_{model_number}_contacts.csv'

    # Skip processing if the file already exists
    if os.path.exists(output_file):
        print(f"File {output_file} exists, omitting...")
        return

    # Use chain_dict if available, otherwise use default values
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

    # Save the contact map
    contacts_df.to_csv(output_file, index=False)
    print(f"Saved contacts in {output_file}.")

if __name__ == "__main__":
    # Get a list of all PDB files in the directory
    pdb_files = [f for f in os.listdir('../../one_model/') if f.endswith('.pdb')]

    # Set the number of workers
    max_workers = max(1, os.cpu_count() - 1)  # Use max cores - 1

    # Run in parallel using ProcessPoolExecutor
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        executor.map(process_pdb_file, pdb_files)

