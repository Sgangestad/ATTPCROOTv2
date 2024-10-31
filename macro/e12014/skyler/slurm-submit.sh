#!/bin/bash
#SBATCH --job-name=run_digi_fiss_job   # Job name
#SBATCH --output=run_digi_fiss_%j.out  # Standard output and error log
#SBATCH --error=run_digi_fiss_%j.err   # Error log
#SBATCH --ntasks=1                     # Number of tasks (processes)
#SBATCH --time=01:00:00                # Time limit hrs:min:sec
#SBATCH --mem=4G                       # Memory limit
#SBATCH --partition=standard           # Partition name

# Set the enviroment variables needed to run the script

# Set the environment variable if needed
export VMCWORKDIR=/path/to/your/vmcworkdir

# Run the ROOT macro
root -l -b -q 'run_digi_fiss.C(0, false)'

# Note: Adjust the arguments (0, false) as needed