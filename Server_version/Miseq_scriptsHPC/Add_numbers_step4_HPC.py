#!/usr/bin/env python3

"""
Optimized Add Numbers Script for PSC
- Uses multiprocessing for large OTU maps
- Utilizes 16 cores and 48GB RAM
- Shows progress bars via tqdm
"""

import sys
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO
from tqdm import tqdm
import multiprocessing as mp

def load_sample_list(list_file):
    with open(list_file) as f:
        return [line.strip().split('\t')[1] for line in f if line.strip() and not line.startswith("MiSeqcode")]

def process_otu_line(line):
    # Worker for one OTU line
    parts = line.strip().split('\t')
    if len(parts) < 2:
        return None
    otu_id = parts[0]
    reads = parts[1:]
    total_reads = sum(int(r.split(";size=")[1]) if ";size=" in r else 1 for r in reads)
    sample_counts = defaultdict(int)
    for read in reads:
        sample = '_'.join(read.replace(" ", "").split('_')[:-1])
        sizes = int(read.split(";size=")[1]) if ";size=" in read else 1
        sample_counts[sample] += sizes
    return (otu_id, total_reads, dict(sample_counts))

def process_numbers(swarm_fasta, swarm_txt, samplelist, threads):
    sample_list = load_sample_list(samplelist)
    output_txt = swarm_txt.parent / ("OTUtable_temp.txt")
    output_fasta = swarm_fasta.parent / (swarm_fasta.stem + "_nosingleton.fas")
    # Read all OTU lines first (I/O bound)
    n_lines = sum(1 for line in open(swarm_txt) if line.strip())
    with open(swarm_txt) as f:
        keep_set = set()
        with open(output_txt, 'w') as out:
            header = ['SWARM', 'TotalReads'] + sample_list
            out.write('\t'.join(header) + '\n')
    
            with mp.Pool(threads) as pool:
                for res in tqdm(pool.imap_unordered(process_otu_line, f, chunksize = 100), total=n_lines, desc="Processing OTUs"):
                    if not res:
                        continue
                    otu_id, total_reads, sample_counts = res
                    if total_reads <= 1:
                        continue
                    keep_set.add(otu_id)
                    row = [otu_id, str(total_reads)]
                    row.extend(str(sample_counts.get(s, 0)) for s in sample_list)
                    out.write('\t'.join(row) + '\n')

    print(f"[INFO] Output written to {output_txt}")
    # Write filtered FASTA
    kept = 0
    with open(output_fasta, 'w') as outfa:
        for record in tqdm(SeqIO.parse(swarm_fasta, 'fasta'), desc="Writing FASTA"):
            otu_id = record.id.split()[0]
            if otu_id in keep_set:
                SeqIO.write(record, outfa, 'fasta')
                kept += 1
    print(f"[INFO] Wrote {kept} sequences to {output_fasta}")


def main():
    process_numbers(Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]), int(sys.argv[4]))

if __name__ == "__main__":
    main()
