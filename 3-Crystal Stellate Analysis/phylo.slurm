#!/bin/bash
#SBATCH --job-name=ISO1_Tree
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16      
#SBATCH --mem=168G              
#SBATCH --time=04:00:00       

set -euo pipefail

module purge
module load GCC/13.2.0  OpenMPI/4.1.6 IQ-TREE/2.3.6      

# Run IQ-TREE. Assumes you have isolated repeat units for a single strain or individual
iqtree2 -s ISO1_onlyALIGNED.fasta -m MFP -B 1000 -bnni -T 16 --prefix ISO1_ML_Tree
