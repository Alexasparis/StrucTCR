# strucTCR

This project aims to predict TCR-peptide binding based on structural data, using position-specific matrices derived from TCR-peptide-MHC interactions. The methodology is structure-based and incorporates the use of MHC (Major Histocompatibility Complex) information, along with sequence clustering to reduce redundancy in the structural dataset.

Project Overview:

The core of the project involves predicting how T-cell receptors (TCRs) bind to peptide-major histocompatibility complex class I (pMHC) complexes. This is achieved by analyzing TCR-pMHC crystals retrieved from the Protein Data Bank (PDB), as well as synthetic structures generated through a data augmentation process using AlphaFold3 from VDJdatabase. The model is based on the creation of a statistical potential derived from position-specific TCR-peptide contacts and TCR-MHC contacts.

Requirements:

Python 3.6 or later.
Required Python packages:
- numpy, pandas, scipy, TCRdist...

You can install the required Python packages with the following command:
pip install -r requirements.txt

Directory Structure
The project follows the following directory structure:
StructTCR/
│
├── data/
│   ├── pdb_files/                # PDB structure files
│   ├── contact_maps/             # Contact map files
│   ├── structures_annotation/    # Annotations and structure-related files
│   ├── all_mhc_seqs.csv          # Dataframe of MHC alleles, sequences 
│
├── models/                       # Pre-trained models
│
├── scripts/                      # Python scripts for processing
│   ├── select_nr_set.py	  # Script to select nr structures.
│   ├── contact_maps_pdb.py	  # Script to extract contact_maps.
│   ├── training.py		  # Script to train the model
│   ├── main.py			  # Main script to execute predictions
│
├── src/                          # Source code (utilities and mapping)
│   ├── utils.py                  # Utility functions
│   ├── potential_calc.py         # Functions to generate and extract potential.
│   ├── find_contact_map.py	  # Functions to find the most similar TCR.
│   ├── extract_contacts.py	  # Functions to extract and filter contacts
│   ├── mapping.py	  	  # Functions to map TCR, MHC, Peptide
│   ├── config.py
│
├── output/                       # Directory for output results
├── input/                        # Directory for input files
├── README.md                     # Project documentation
└── requirements.txt              # Python dependencies

How to Run the Script






