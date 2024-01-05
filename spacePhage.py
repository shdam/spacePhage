import subprocess
import argparse
import csv
import os
from Bio import SeqIO

def split_sequence(sequence, length):
    return [sequence[i:i+length] for i in range(len(sequence) - length + 1)]

def extract_spacer(sequence, start, end):
    return sequence[start-1:end]

def run_blast(phage_file, db_name, output_file, min_score=37, spacer_length=37):
    # Splitter phage sekvensen i 37-mers
    phage_seq = SeqIO.read(phage_file, 'fasta')
    spacers = split_sequence(str(phage_seq.seq), spacer_length)
    query_sequences = {f'{phage_seq.id}_pos{i+1}': spacer for i, spacer in enumerate(spacers)}

    # Midlertidig fil til at gemme de splittede sekvenser
    temp_query_file = 'temp_query.fasta'
    with open(temp_query_file, 'w') as file:
        for id, seq in query_sequences.items():
            file.write(f'>{id}\n{seq}\n')

    # Tilpasset output format: query id, subject id, alignment length, start, end, score, e-value
    outfmt = "6 qseqid sseqid length qlen sstart send score evalue"

    # Kør BLAST-kommando og fang output
    cmd = f'blastn -query {temp_query_file} -db {db_name} -outfmt "{outfmt}"'
    result = subprocess.run(cmd, capture_output=True, text=True, shell=True)

    # Behandle og gem output som CSV
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerow(['qseqid', 'sseqid', 'length', 'qlen', 'sstart', 'send', 'score', 'evalue', 'spacer_sequence']) # Header
        for line in result.stdout.splitlines():
            row = line.split('\t')
            if int(row[6]) >= min_score:  # Tjekker om scoren er 37 eller højere
                spacer_sequence = query_sequences[row[0]]
                writer.writerow(row + [spacer_sequence])

    # Fjern midlertidig fil
    os.remove(temp_query_file)

def main():
    parser = argparse.ArgumentParser(description='Run BLAST search with phage sequence, split into spacers, and output results as a CSV file.')
    parser.add_argument('-p', '--phage', required=True, help='FASTA file containing the phage sequence')
    parser.add_argument('-d', '--database', required=True, help='Name of the BLAST database to search against')
    parser.add_argument('-o', '--output', required=True, help='CSV file to write BLAST search results')
    parser.add_argument('-s', '--min_score', type=int, default=37, help='Minimum score to include in the results')

    args = parser.parse_args()
    run_blast(args.phage, args.database, args.output, args.min_score)

if __name__ == "__main__":
    main()

