
#!/bin/bash

# Check for correct number of arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <phage_fasta_file> <genome_fasta_folder>"
    exit 1
fi

PHAGE_FASTA=$1
GENOME_FOLDER=$2
DB_NAME="blastdb/blast_database"
OUTPUT_FILE="blast_results.csv"

# Create BLAST database from genome FASTA files
python3 createBlastDB.py -i "$GENOME_FOLDER" -o "$DB_NAME"

# Run analysis on phage FASTA file
python3 spacePhage.py -p "$PHAGE_FASTA" -d "$DB_NAME" -o "$OUTPUT_FILE"

echo "BLAST analysis completed. Results saved in $OUTPUT_FILE"
echo "Cleaning temporary files"
rm -rf blastdb
rm ${GENOME_FOLDER}/combined.fasta