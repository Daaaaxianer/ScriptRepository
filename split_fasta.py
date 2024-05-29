import argparse
import os
from Bio import SeqIO

def split_sequences(input_file, output_dir):
    for record in SeqIO.parse(input_file, "fasta"):
        sequence_id = record.id
        output_filename = os.path.join(output_dir, f"{sequence_id}.fasta")
        
        with open(output_filename, "w") as out_handle:
            SeqIO.write(record, out_handle, "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split FASTA sequences into individual files.")
    parser.add_argument("-fa", "--fastaSeq", type = str, help = "Input FASTA file")
    parser.add_argument("-odir","--output_dir", default="split_out", help="Output directory (default: split_out)")
    args = parser.parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    split_sequences(args.fastaSeq, args.output_dir)