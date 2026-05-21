#!/bin/bash
#SBATCH --job-name=juicer_CantonS
#SBATCH --output=juicer_%j.out
#SBATCH --error=juicer_%j.err
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=168G
#SBATCH --time=67:00:00     #Juicer takes a long time!

set -euo pipefail

module purge
module load GCC/12.2.0 BWA/0.7.17 SAMtools/1.17 Juicer/1.6-Java-11-CUDA-11.7.0-CPU
juicer=/sw/eb/sw/Juicer/1.6-GCCcore-12.2.0-Java-11-CUDA-11.7.0-CPU/CPU/scripts/juicer.sh
juicer_dir=/sw/eb/sw/Juicer/1.6-GCCcore-12.2.0-Java-11-CUDA-11.7.0-CPU/CPU

# Run Juicer with absolute paths
$juicer \
  -g CantonS \
  -d $(pwd) \
  -p $(pwd)/Delgado.chrom.sizes \
  -y $(pwd)/CANTONJUICE_DpnII.txt \
  -z $(pwd)/wholegenomeDelgado.fasta \
  -t 32 \
  -D $juicer_dir
