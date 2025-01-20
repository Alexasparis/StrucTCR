import subprocess
import argparse

# EXAMPLE python run_validations.py -m main.py -i control_fold -mfp ./model/TCR-p_all -mfm ./model/TCR_MHC_all.csv -w 8

# Set up argparse to handle command-line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description="Run validation for each fold.")
    parser.add_argument('-m', '--main', type=str, help='Main script to run', required=True)
    parser.add_argument('-i', '--input_folder', type=str, help='Input folder path', required=True)
    parser.add_argument('-mfp', '--mfp_name', type=str, help='MFP file name', required=True)
    parser.add_argument('-mfm', '--mfm_name', type=str, help='MFM file name', required=True)
    parser.add_argument('-w', '--workers', type=int, help='max workers', required=True)
    return parser.parse_args()

# Function to execute the command in bash
def run_command(command):
    print(f"Running: {command}")
    result = subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.stdout:
        print("Standard Output:")
        print(result.stdout)
    else:
        print("No standard output generated.")
    
    # Print errors if there are any
    if result.stderr:
        print("Standard Error:")
        print(result.stderr)

# Main function to handle running the command for each fold
def main():
    # Parse arguments passed through the command line
    args = parse_arguments()
    
    # Loop through each fold (1 to 5)
    for fold in range(1, 6):
        # Dynamically create the arguments for each fold
        input_folder = f"{args.input_folder}_{fold}"
        output_folder = f"{args.input_folder}_{fold}"
        output_filename = f"scores_{args.input_folder}_{fold}"
        mfp_name = f"{args.mfp_name}_f{fold}"
        mfm_name=args.mfm_name
        workers=args.workers
        
        # Construct the command with parsed arguments
        command = f"python run.py -m {args.main} -i {input_folder} -o {output_folder} -ofn {output_filename} -mfp {mfp_name} -mfm {mfm_name} -w {workers}"
        
        # Call the function to run the constructed command
        run_command(command)

if __name__ == "__main__":
    main()


