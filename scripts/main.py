#!/usr/bin/env python3

# This script is used to train the model given a database of contact_maps.


# Example of execution: python contact_maps_pdb.py -pdb pdb_files -out contact_maps -workers 7

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))