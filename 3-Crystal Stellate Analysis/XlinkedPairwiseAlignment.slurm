#!/bin/bash
#SBATCH --job-name=mafft_ste_x
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=360G
#SBATCH --time=72:00:00 
#SBATCH --output=mafft_ste_x_%j.out
#SBATCH --error=mafft_ste_x_%j.err

set -euo pipefail

module purge
module load GCC/13.2.0 MAFFT/7.526-with-extensions

# Step 1: Pool all the final alignment outputs into one master FASTA
MASTER_FASTA="All_Strains_X_Stellate_Units_Master.fasta"
echo "Concatenating individual strain FASTAs..."
cat Final_*_X_Stellate_BLAST_Units.fasta > $MASTER_FASTA

TOTAL_SEQ=$(grep -c "^>" $MASTER_FASTA)

# Step 2: Run MAFFT L-INS-i
# L-INS-i is invoked using --localpair and --maxiterate 1000
ALIGNED_FASTA="All_Strains_X_Stellate_Aligned.fasta"

echo "Running MAFFT L-INS-i alignment on $SLURM_CPUS_PER_TASK threads..."
mafft --localpair \
      --maxiterate 1000 \
      --thread $SLURM_CPUS_PER_TASK \
      $MASTER_FASTA > $ALIGNED_FASTA

echo "Alignment complete!"
date
