import subprocess
import argparse
import csv
import os
from Bio import SeqIO

def split_sequence(sequence, length):
    return [sequence[i:i+length] for i in range(len(sequence) - length + 1)]

def extract_sequence(sequence, start, end, flank_length=3):
    # Henter flankerende nukleotider
    pre_flank = sequence[max(0, start-flank_length-1):max(0, start-1)]
    post_flank = sequence[end:min(end+flank_length, len(sequence))]
    print(pre_flank)
    print(post_flank)
    return sequence[start-1:end], pre_flank, post_flank

def run_blast(phage_file, db_name, genome_folder, output_file, min_score=36, spacer_length=36):
    # Load and process the phage sequence
    phage_seq = SeqIO.read(phage_file, 'fasta')
    spacers = split_sequence(str(phage_seq.seq), spacer_length)
    query_sequences = {f'{phage_seq.id}_pos{i+1}': spacer for i, spacer in enumerate(spacers)}

    # Create a temporary file for the query sequences
    temp_query_file = 'temp_query.fasta'
    with open(temp_query_file, 'w') as file:
        for id, seq in query_sequences.items():
            file.write(f'>{id}\n{seq}\n')

    # Run BLAST and capture the output
    outfmt = "6 qseqid sseqid length qlen sstart send score evalue"
    cmd = f'blastn -query {temp_query_file} -db {db_name} -outfmt "{outfmt}"'
    result = subprocess.run(cmd, capture_output=True, text=True, shell=True)

    # Read combined FASTA file
    combined_fasta_path = os.path.join(genome_folder, 'combined.fasta')
    combined_sequences = SeqIO.to_dict(SeqIO.parse(combined_fasta_path, 'fasta'))

    # Process and save the output as CSV
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerow(['qseqid', 'sseqid', 'length', 'qlen', 'sstart', 'send', 'score', 'evalue', 'spacer_sequence', '3nt_pre', '3nt_post', 'bacterial_pre', 'bacterial_post'])
        for line in result.stdout.splitlines():
            row = line.split('\t')
            if int(row[6]) >= min_score:
                spacer_sequence = query_sequences[row[0]]
                spacer_id = row[0]
                spacer_start = int(spacer_id.split('_pos')[-1]) - 1
                pre_flank = phage_seq.seq[max(0, spacer_start - 3):spacer_start]
                post_flank = phage_seq.seq[spacer_start + spacer_length:spacer_start + spacer_length + 3]
                #writer.writerow(row + [phage_seq.seq[spacer_start:spacer_start + spacer_length], pre_flank, post_flank])

                if row[1] in combined_sequences:
                    if int(row[5]) < int(row[4]):
                        start=int(row[5])-1
                        end=int(row[4])
                    else:
                        start=int(row[4])-1
                        end=int(row[5])
                    bacterial_pre = str(combined_sequences[row[1]].seq[start-29:start])
                    bacterial_post = str(combined_sequences[row[1]].seq[end:end+29])
                else:
                    bacterial_sequence = 'Sequence not found'
                writer.writerow(row + [phage_seq.seq[spacer_start:spacer_start + spacer_length], pre_flank, post_flank] + [bacterial_pre, bacterial_post])

    # Remove the temporary file
    os.remove(temp_query_file)


def main():
    parser = argparse.ArgumentParser(description='Run BLAST search with phage sequence, split into spacers, and output results as a CSV file.')
    parser.add_argument('-p', '--phage', required=True, help='FASTA file containing the phage sequence')
    parser.add_argument('-d', '--database', required=True, help='Name of the BLAST database to search against')
    parser.add_argument('-o', '--output', required=True, help='CSV file to write BLAST search results')
    parser.add_argument('-s', '--min_score', type=int, default=36, help='Minimum score to include in the results')
    parser.add_argument('-g', '--genome_folder', required=True, help='Folder containing the combined genome FASTA file')

    args = parser.parse_args()
    run_blast(args.phage, args.database, args.genome_folder, args.output, args.min_score)

if __name__ == "__main__":
    main()
