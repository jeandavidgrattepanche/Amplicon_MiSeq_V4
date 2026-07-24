#!/usr/bin/env python3

from Bio import SeqIO
from pathlib import Path
import sys

def load_sample_mapping(list_file):
    sample_map = {}
    with open(list_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                sample_map[parts[0]] = parts[1]
    return sample_map

def pool_sequences(fasta_file, sample_map, output_file):
    sample_code = Path(fasta_file).stem.split('_')[0]
    sample_name = sample_map.get(sample_code, sample_code)

    records = []
    with open(fasta_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            record.id = f"{sample_name}_{record.id}"
            record.description = ""
            records.append(record)

    with open(output_file, "a") as out_handle:
        SeqIO.write(records, out_handle, "fasta")

def process_pool_rename(fasta_file, list_file):
    sample_map = load_sample_mapping(list_file)
    output_file = str(Path(fasta_file).parent / "readpooled.fas")
    pool_sequences(fasta_file, sample_map, output_file)
    return output_file

def main():
    if len(sys.argv) != 3:
        print("Usage: python 1_pool_rename_vHPC.py <seqfile.fasta> <List_samples.txt>")
        sys.exit(1)
    process_pool_rename(sys.argv[1], sys.argv[2])
#     print(f"[INFO] Pooled {fasta_file} into {output_file}")

if __name__ == "__main__":
    main()
