
# Find CRISPR Spacers from Phage Genome Tool

## Purpose
This tool automates identifying CRISPR spacers from a given phage genome through BLAST searches against a genome database.

## Installation
Run `make install` to install necessary Python dependencies and to set up the scripts. Ensure BLAST is installed on your system.

## How to run spacePhage

- Prepare your genome data in FASTA format in a folder (e.g., `test_db`).

```bash
spacePhage test.fasta test_db
```

### Output
- The output is a CSV file (`blast_results.csv`) with the BLAST search results, including CRISPR spacers identified from the phage genome.

## Running the scripts manually:

### Prerequisites
- Python 3
  - `brew install python3`
- NCBI's BLAST+ software
  - `brew install blast`
- BioPython library
  - `pip3 install biopython`

### Steps

#### 1. Create a BLAST Database
   - Prepare your genome data in FASTA format in a folder (e.g., `test_db`).
   - Use `createBlastDB.py` to create a BLAST database from these FASTA files.
   - This script requires the path to your FASTA files folder and a name for the new BLAST database.
   ```bash
   DB=blastdb/test_blastdb
   python3 createBlastDB.py -i test_db -o $DB
   ```

#### 2. Analyze Phage Genome
   - Prepare your phage genome sequence in a FASTA file (e.g., `test.fasta`).
   - Use `spacePhage.py` to analyze the phage genome. This script splits the phage genome into 37-mers and performs a BLAST search against the database.
   - Provide the path to your phage FASTA file, the BLAST database, and an output CSV file name.
   ```bash
   python3 spacePhage.py -p test.fasta -d $DB -o blast_results.csv
   ```




