import argparse
from Bio import SeqIO

def split_sequence(sequence, length):
    return [(i, sequence[i:i+length]) for i in range(len(sequence) - length + 1)]

def main():
    parser = argparse.ArgumentParser(description='Split a fasta sequence into 37 bp spacers with position information.')
    parser.add_argument('-i', '--input', required=True, help='Input fasta file')
    parser.add_argument('-o', '--output', required=True, help='Output file for the spacers')

    args = parser.parse_args()

    fag_seq = SeqIO.read(args.input, 'fasta')
    input_id = fag_seq.id
    spacers = split_sequence(str(fag_seq.seq), 37)

    with open(args.output, 'w') as file:
        for start, spacer in spacers:
            spacer_id = f'>{input_id}_pos{start+1}\n'
            file.write(spacer_id + spacer + '\n')

if __name__ == "__main__":
    main()

