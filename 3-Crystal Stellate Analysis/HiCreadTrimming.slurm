#!/bin/bash
#SBATCH --job-name=fastp_ISO1
#SBATCH --output=fastp_%j.out
#SBATCH --error=fastp_%j.err
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00

module purge
module load GCC/13.2.0 fastp/0.23.4

# Run fastp
fastp \
  -i SRR33554829_1.fastq \
  -I SRR33554829_2.fastq \
  -o SRR33554829_1.trimmed.fastq \
  -O SRR33554829_2.trimmed.fastq \
  --detect_adapter_for_pe \
  --length_required 100 \
  -w 16 \
  -h SRR33554829_report.html \
  -j SRR33554829_report.json

echo "fastp processing complete for ISO-1"
