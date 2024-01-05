#!/bin/bash
# Split phage into 37-mers
python3 splitPhage.py -i phage.fasta -o phage.split.fasta

# Extract a specific 37-mer
index=34
python3 extractSpacer.py -i phage.split.fasta -s $index -o output.fasta
