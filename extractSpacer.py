import argparse
from Bio import SeqIO

def extract_sequence_by_index(filename, index):
    for i, record in enumerate(SeqIO.parse(filename, 'fasta')):
        if i == index - 1:  # Justerer til 0-baseret indeksering
            return record.id, record.seq
    return None, None

def main():
    parser = argparse.ArgumentParser(description='Extract a specific sequence from a FASTA file based on its index.')
    parser.add_argument('-i', '--input', required=True, help='Input fasta file')
    parser.add_argument('-s', '--index', type=int, required=True, help='Index (1-based) of the sequence to extract')
    parser.add_argument('-o', '--output', required=True, help='Output file for the extracted sequence')

    args = parser.parse_args()

    sequence_id, sequence = extract_sequence_by_index(args.input, args.index)
    if sequence is None:
        print("No sequence found at the specified index.")
        return

    with open(args.output, 'w') as file:
        file.write(f'>{sequence_id}\n')
        file.write(str(sequence) + '\n')

if __name__ == "__main__":
    main()

