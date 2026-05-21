# Project Component Name

Brief description of this part of the project.

Example:

This directory contains scripts and workflows used for repeat annotation of the Drosophila melanogaster Y chromosome assembly.

---

# Overview

Short paragraph describing:
- what this analysis step does
- what inputs it expects
- what outputs it generates
- where it fits in the larger workflow

Example:

These scripts identify and summarize repetitive elements in the assembled Y chromosome. The workflow begins with RepeatMasker annotation and ends with generation of summary statistics and manuscript figures.

---

# Directory Structure

```text
repeat_annotation/
├── README.md
├── config/
├── scripts/
├── intermediate/
├── results/
└── figures/
```

---

# Workflow Summary

1. Run RepeatMasker
2. Parse RepeatMasker output
3. Generate repeat summaries
4. Create manuscript figures

---

# Scripts

## 1. RepeatMasker Processing

### run_repeatmasker.sh

Runs RepeatMasker on the assembled genome sequence.

**Input**
- assembly.fa

**Output**
- masked genome FASTA
- RepeatMasker annotation files

**Example usage**

```bash
sbatch run_repeatmasker.sh
```

---

### parse_repeatmasker.py

Parses RepeatMasker output and summarizes repeat family composition.

**Requirements**
- Python 3
- pandas

**Input**
- repeats.out

**Output**
- repeat_summary.tsv

**Example usage**

```bash
python parse_repeatmasker.py repeats.out
```

---

# 2. Figure Generation

### plot_repeat_content.R

Generates repeat composition plots used in the manuscript.

**Input**
- repeat_summary.tsv

**Output**
- PDF and PNG figures

**Example usage**

```bash
Rscript plot_repeat_content.R
```

---

# Software Requirements

| Software | Version |
|---|---|
| RepeatMasker | 4.1 |
| Python | 3.10 |
| R | 4.3 |

---

# Notes

Additional notes, caveats, or warnings.

Example:
- Large intermediate files are not tracked in GitHub.
- Raw sequencing data are stored separately.
- SLURM scripts assume TAMU HPRC environment.

---

# Contact

Name and email for questions related to this workflow.