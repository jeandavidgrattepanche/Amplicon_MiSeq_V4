#!/bin/bash
#SBATCH --job-name=miseq_p1
#SBATCH --output=miseq_p1_%j.out
#SBATCH --error=miseq_p1_%j.err

#SBATCH --partition=RM
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64GB
#SBATCH --time=24:00:00

: "${PROJECT:?Set PROJECT to your pipeline directory}"
SAMPLE_FILE="${1:?Provide a sample-list file}"
RAW_DATA_DIR="${2:?Provide the raw-data directory}"


# Go to submission directory
cd $SLURM_SUBMIT_DIR

set -euo pipefail

echo "====================================="
echo "Job started on $(date)"
echo "Running on $(hostname)"
echo "====================================="

# Clean module environment
module purge

# Load modules
module load anaconda3
module load bbmap

# Activate environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate Amplicon311

# Define project path
export PROJECT=ADD_PATH_TO_PROJECT


# Check software
echo "Python:"
which python

echo "vsearch:"
which vsearch

echo "swarm:"
which swarm

echo "bbduk:"
which bbduk.sh

echo "bbmerge:"
which bbmerge-auto.sh

python --version
vsearch --version
swarm --version
bbduk.sh --version

echo "SLURM_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK"
echo "SLURM_MEM_PER_NODE=$SLURM_MEM_PER_NODE"

# Run pipeline
python $PROJECT/Amplicon/MiSeq_pipeline_part1_PSC.py \
    $SAMPLE_FILE \
    $RAW_DATA_DIR \
    --minlen 400 \
    --sample_threads 8 \
    --tool_threads 2

echo "====================================="
echo "Job finished on $(date)"
echo "====================================="