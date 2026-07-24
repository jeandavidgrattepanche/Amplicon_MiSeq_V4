#!/usr/bin/env python3

import sys
import time
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO
from tqdm import tqdm

def parse_otu_table(otu_file, read_cutoff):
    abddict = defaultdict(list)
    IDdict = {}
    abunlist = set()

    with open(otu_file) as f:
        for line in f:
            if line.startswith("OTU"):
                parts = line.strip().split("\t")
                otu_id = parts[0]
                count = int(parts[1])
                if count >= read_cutoff:
                    IDdict[otu_id] = count
                    abddict[count].append(otu_id)
                    abunlist.add(count)
    return IDdict, abddict, sorted(abunlist, reverse=True)

def extract_sequences(seqfile, valid_ids):
    fastadict = {}
    for record in SeqIO.parse(seqfile, "fasta"):
        if record.id in valid_ids:
            fastadict[record.id] = str(record.seq)
    return fastadict

def write_output(chimera_dir, IDdict, abddict, abun_sorted, fastadict):
    out_path = Path(chimera_dir) / ("Seq_reads_test.fas")
    with open(out_path, "w") as out:
        for size in tqdm(abun_sorted, desc="Writing filtered sequences"):
            for otu_id in abddict[size]:
                if otu_id in fastadict:
                    seq = fastadict[otu_id]
                    out.write(f">{otu_id};size={IDdict[otu_id]};\n{seq}\n")
def pre_uchime( chimera_dir, seqfile, OTUfile, readcutoff):
  
    print("[INFO] Parsing OTU table...")
    IDdict, abddict, abun_sorted = parse_otu_table(OTUfile, readcutoff)

    print(f"[INFO] Extracting sequences (total OTUs: {len(IDdict)})...")
    fastadict = extract_sequences(seqfile, IDdict)

    print("[INFO] Writing filtered FASTA file...")
    write_output(chimera_dir, IDdict, abddict, abun_sorted, fastadict)

def post_uchime(seqfile):
    outpath = seqfile.parent.parent / "chimeras" / "Seq_reads_nochimera_nosingleton_renamed.fas"
    outpath.parent.mkdir(parents=True, exist_ok=True)

    with open(outpath, "w") as out_handle:
        for record in SeqIO.parse(seqfile, "fasta"):
            header = record.description.replace(';size=', '_').replace(';', 'r')
            out_handle.write(f">{header}\n{record.seq}\n")

    print(f"[INFO] Output written to: {outpath}")


def main():
#     pre_uchime(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]))
#     post_uchime(sys.argv[1])
    print("Use pre_uchime() and post_uchime() through pipeline wrapper.")

if __name__ == "__main__":
    main()


