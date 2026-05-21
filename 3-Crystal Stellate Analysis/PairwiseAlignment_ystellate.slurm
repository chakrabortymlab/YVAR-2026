#!/bin/bash
#SBATCH --job-name=mafft_linsi
#SBATCH --output=mafft_combined_%j.out
#SBATCH --error=mafft_combined_%j.err
#SBATCH --time=24:00:00          
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48       
#SBATCH --mem=360G               
set -euo pipefail

module purge
module load GCC/13.2.0 MAFFT/7.526-with-extensions

echo "Starting MAFFT L-INS-i alignment at $(date)"

# 3. Execute MAFFT
# --localpair and --maxiterate 1000 triggers the L-INS-i algorithm
mafft --localpair --maxiterate 1000 --thread $SLURM_CPUS_PER_TASK PCKR_SuSte_Combined.fasta > PCKR_SuSte_Combined_Aligned.fasta

echo "Alignment complete! Saved to PCKR_SuSte_Combined_Aligned.fasta at $(date)"
