import os

# Dir paths
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

DATA_DIR = os.path.join(BASE_DIR, "data")
PDB_DIR = os.path.join(DATA_DIR, "pdb_files")
CONTACT_MAPS_DIR = os.path.join(DATA_DIR, "contact_maps")
STRUCTURES_ANNOTATION_DIR = os.path.join(DATA_DIR, "structures_annotation")

# Modelos y resultados
MODELS_DIR = os.path.join(BASE_DIR, "models")

# Scripts
SCRIPTS_DIR = os.path.join(BASE_DIR, "scripts")

# Directorio de salida de los resultados
OUTPUT_DIR = os.path.join(BASE_DIR, "output")
INPUT_DIR = os.path.join(BASE_DIR, "input")





    
