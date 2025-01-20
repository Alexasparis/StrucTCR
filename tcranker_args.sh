#!/bin/bash
#SBATCH --job-name=tcranker      # Job name
#SBATCH --output=logs/tcranker_%j.out # stdout file
#SBATCH --error=logs/tcranker_%j.err  # stderr file
#SBATCH --ntasks=1                # Number of tasks
#SBATCH --cpus-per-task=50        # CPUs per task
#SBATCH --time=05:00:00           # Execution time
#SBATCH -D .
#SBATCH --qos=gp_bscls            # QoS (dbug, gp_bscls)

RUN_SCRIPT=$1                   
MAIN_SCRIPT=$2                   
INPUT_FOLDER=$3                                
MODEL_FILE_P=$4           
MODEL_FILE_M=$5                   
WORKERS=${6:-50}                  


# Load conda environment
module load miniforge
source activate anarci

# Run script
python $RUN_SCRIPT -m $MAIN_SCRIPT \
                    -i $INPUT_FOLDER \
                    -mfp $MODEL_FILE_P \
                    -mfm $MODEL_FILE_M \
                    -w $WORKERS

# Example validations: sbatch -A bsc72 -q gp_bscls run_validations.py main.py control_fold_1 ./model/TCR-p_all ./model/TCR_MHC_all.csv 50
# Example sequences: sbatch -A bsc72 -q gp_bscls run_seqs.py main.py sequences ./model/TCR-p_all ./model/TCR_MHC_all.csv 50
# Deactivate env
conda deactivate
