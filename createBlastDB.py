import os
import subprocess
import argparse

def combine_fasta_files(input_folder, output_file):
    with open(output_file, 'w') as outfile:
        for filename in os.listdir(input_folder):
            if filename.endswith('.fasta'):
                filepath = os.path.join(input_folder, filename)
                with open(filepath, 'r') as infile:
                    for line in infile:
                        outfile.write(line)

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

