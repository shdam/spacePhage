import argparse
from Bio import SeqIO

def split_sequence(sequence, length):
    return [sequence[i:i+length] for i in range(0, len(sequence), length)]

def main():
    parser = argparse.ArgumentParser(description='Split a fasta sequence into smaller sequences.')
    parser.add_argument('-i', '--input', required=True, help='Input fasta file')
    parser.add_argument('-o', '--output', required=True, help='Output file for the split sequences')

    args = parser.parse_args()

    fag_seq = SeqIO.read(args.input, 'fasta')
    splitted_sequences = split_sequence(str(fag_seq.seq), 37)

    with open(args.output, 'w') as file:
        for seq in splitted_sequences:
            file.write(seq + '\n')

if __name__ == "__main__":
    main()

