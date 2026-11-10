#!/bin/bash
#SBATCH --job-name=test_job       # Job name
#SBATCH --output=test_job.out     # Standard output file
#SBATCH --error=test_job.err      # Standard error file
#SBATCH --time=00:05:00           # Time limit (HH:MM:SS)
#SBATCH --partition=general       # Partition (change as needed)
#SBATCH --nodes=1                 # Number of nodes
#SBATCH --ntasks=1                # Number of tasks
#SBATCH --cpus-per-task=1         # Number of CPU cores per task


# Run a test command
echo "Running on $(hostname)"
sleep 30
echo "Job completed successfully."
python sample_code.py
