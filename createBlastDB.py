import os
import subprocess
import argparse
import glob

def combine_fasta_files(input_folder, output_file):
    # Find all FASTA files (.fasta and .fa) in the input folder
    fasta_files = glob.glob(os.path.join(input_folder, '*.fasta')) + glob.glob(os.path.join(input_folder, '*.fa'))

    # Concatenate the content of all found FASTA files into one output file
    with open(output_file, 'w') as outfile:
        for file in fasta_files:
            with open(file, 'r') as infile:
                outfile.write(infile.read())
                outfile.write('\n')  # Ensure separation between files


def create_blast_db(combined_fasta, db_name, db_type='nucl'):
    cmd = f'makeblastdb -in {combined_fasta} -dbtype {db_type} -out {db_name}'
    subprocess.run(cmd, shell=True)

def main():
    parser = argparse.ArgumentParser(description='Create a BLAST database from a folder of FASTA files.')
    parser.add_argument('-i', '--input_folder', required=True, help='Folder containing FASTA files')
    parser.add_argument('-o', '--db_name', required=True, help='Name for the BLAST database')
    parser.add_argument('-t', '--db_type', choices=['nucl', 'prot'], default='nucl', help='Type of the database (nucl or prot)')

    args = parser.parse_args()

    # Kombinerer FASTA-filer til én enkelt fil
    combined_fasta = os.path.join(args.input_folder, 'combined.fasta')
    combine_fasta_files(args.input_folder, combined_fasta)

    # Opretter BLAST-database
    create_blast_db(combined_fasta, args.db_name, args.db_type)

if __name__ == "__main__":
    main()

