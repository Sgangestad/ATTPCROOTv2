#!/bin/bash
#SBATCH --job-name=sim_fission_fit
#SBATCH --output=log/sim_fission_%j.out
#SBATCH --error=log/sim_fission_%j.err
#SBATCH --partition=cs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=1G
#SBATCH --time=UNLIMITED
#SBATCH --chdir=/home/faculty/aanthony/fission/adam/ATTPCROOTv2/macro/e12014/skyler

# This macro will generate fisison events and then fit them using the MC fitting code
# The number of events is set by the NUM_EVENTS variable
# The number of threads to use when running each job is set by the SLURM_CPUS_ON_NODE variable
# The run number is set to be the same as the job number (the %j variable)
# The output will be saved to the log directory of the current directory

NUM_EVENTS=1
SAVE_RAW_EVENT=true

echo "For run number  ${SLURM_JOB_ID}"
echo "Starting to simulate fission events at  $(date)"

# The following command will generate the fission events
root -l -q -b "run_sim_fiss_slurm.C(${SLURM_JOB_ID}, ${NUM_EVENTS})"

echo "Finished to simulating fission events at  $(date)"
echo "Starting to fit fission events at  $(date) with ${SLURM_CPUS_ON_NODE} threads"

# The following command will fit the fission events
root -l -q -b "run_sim_fit_slurm.C(${SLURM_JOB_ID}, ${SLURM_CPUS_ON_NODE}, ${SAVE_RAW_EVENT})"

echo "Finished fitting fission events at  $(date)"