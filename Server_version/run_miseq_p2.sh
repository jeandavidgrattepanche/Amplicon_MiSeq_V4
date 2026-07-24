#!/bin/bash
#SBATCH --job-name=miseq_p2
#SBATCH --output=miseq_p2_%j.out
#SBATCH --error=miseq_p2_%j.err

#SBATCH --partition=RM
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64GB
#SBATCH --time=48:00:00


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

echo "water:"
which water

echo "needle:"
which needle

python --version
vsearch --version
swarm --version
command -v water >/dev/null 2>&1 || { echo "water not found"; exit 1; }
command -v needle >/dev/null 2>&1 || { echo "needle not found"; exit 1; }
water --version
needle --version

echo "SLURM_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK"
echo "SLURM_MEM_PER_NODE=$SLURM_MEM_PER_NODE"
echo "Working directory: $(pwd)"

echo "====================================="
echo "Pipeline started at $(date)"
echo "====================================="
# Run pipeline
python $PROJECT/Amplicon/MiSeq_pipeline_part2_PSC.py \
    $SAMPLE_FILE \
    $RAW_DATA_DIR \
    --threads 16 \
    --swarm_d 1 \
    --read_cutoff 100 \
    --tool water \
    --assign_taxonomy \
    --idmin 97 \
    --qcov 80 \
    --taxa P \
    --filter_taxa \
    --force_pool \
    --force_swarm \
    --force_chimera \
    --force_contaminant \
    --force_taxonomy
    
echo "====================================="
echo "Job finished on $(date)"
echo "====================================="