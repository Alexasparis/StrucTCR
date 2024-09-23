# Train_model

# This is the main function to train the model and extract the potential matrices from a given pdb_file database composed by TCR-pMHC complexes.
#Example of usage
#python3 -p ./database [- d 5 -out_p TCRen_TCR_p.csv -out_m TCRen_TCR_MHC.csv]

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
    parser.add_argument("-d", "--distance", type=float, default=5.0, help='Distance threshold for contacts. Default 5 A')
    parser.add_argument("-out_p", "--output_p", type=str, default='TCRen_TCR_p.csv', help='Output CSV file name for TCR peptide potential.')
    parser.add_argument("-out_m", "--output_m", type=str, default='TCRen_TCR_MHC.csv', help='Output CSV file name for TCR MHC potential.')

    args = parser.parse_args()
    
    # Ensure the output directory exists
    output_dir = 'model'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    pdb_files = [os.path.join(args.pdbs, f) for f in os.listdir(args.pdbs) if f.endswith('.pdb')]
    
    # Extract contacts from PDB files
    print("Extracting contacts from PDB files")
    contacts = extract_contacts(pdb_files, distance=args.distance)
    contacts.to_csv(os.path.join(output_dir, 'contacts_training.csv'), index=False)
    
    # Filter contacts
    print("Filtering contacts")
    contacts_TCR_p, contacts_TCR_MHC = filter_contacts(contacts)
    
    # Calculate TCRen potentials
    print("Calculating potentials")
    print("Calculating TCR_peptide potential")
    data_TCR_p = calculate_TCRen(contacts_TCR_p, peptide=True)
    print("Calculating TCR_MHCI potential")
    data_TCR_MHC = calculate_TCRen(contacts_TCR_MHC, peptide=False)
    
    # Save to CSV
    print("Saving results")
    data_TCR_p.to_csv(os.path.join(output_dir, args.output_p), index=False)
    data_TCR_MHC.to_csv(os.path.join(output_dir, args.output_m), index=False)
    
    print(f"TCRen potential for TCR peptide saved to {os.path.join(output_dir, args.output_p)}")
    print(f"TCRen potential for TCR MHC saved to {os.path.join(output_dir, args.output_m)}")

if __name__ == "__main__":
    main()


