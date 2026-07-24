#!/usr/bin/env python3

"""
3_postSwarm_vHPC.py - Optimized for Mac Studio M4 Max

This script post-processes SWARM clustering results by reconstructing the representative sequence and member sequence lists
for each OTU cluster. It uses parallel processing and progress tracking to fully leverage Mac Studio M4 Max performance.
"""

from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from Bio import SeqIO
from tqdm import tqdm
import sys

GLOBAL_READDICT = None
GLOBAL_UNIDICT = None


def parse_dereplicated_sequences(fasta_path):
    return {record.description: str(record.seq) for record in SeqIO.parse(fasta_path, 'fasta')}


def parse_uc_map(map_path):
    unidict = {}
    with open(map_path) as infile:
        for line in infile:
            parts = line.strip().split('\t')
            if parts[0] == 'S':
                unidict.setdefault(parts[8].split(';')[0], [parts[8].strip()])
            elif parts[0] == 'H':
                unidict.setdefault(parts[9].split(';')[0], []).append(parts[8].strip())
    return unidict

def init_worker(readdict, unidict):

    global GLOBAL_READDICT
    global GLOBAL_UNIDICT

    GLOBAL_READDICT = readdict
    GLOBAL_UNIDICT = unidict
    

def process_swarm_line(i, line):
    representative = line.strip().split(" ")[0]
#     with outseq_path.open('a') as outseq, outlist_path.open('a') as outlist:
    seq_text = (f">OTU{i}\t{representative}\n{GLOBAL_READDICT.get(representative, '')}\n")
    amplicons = line.strip().split(" ")
    seqs = ['\t'.join( GLOBAL_UNIDICT.get(amp.split(';')[0], [])) for amp in amplicons]
    list_text = (f"OTU{i}\t" + "\t".join(seqs) + "\n")
    return seq_text, list_text

def process_postswarm(swarm_output_path, derep_map_path,  derep_seq_path, threads) :
    output_dir = derep_map_path.parent
    outseq_path = output_dir / "SWARM_postout.fas"
    outlist_path = output_dir / "SWARM_postout.txt"

    outseq_path.unlink(missing_ok=True)
    outlist_path.unlink(missing_ok=True)

    print("[INFO] Parsing dereplicated sequences and map file...")
    readdict = parse_dereplicated_sequences(derep_seq_path)
    unidict = parse_uc_map(derep_map_path)
    with open(swarm_output_path) as f:
        swarm_lines = f.readlines()
#     for i, line in tqpm(enumerate(swarm_lines, start=1),total=len(swarm_lines), desc="Processing OTUs")
#         process_swarm_line(1, line, readdict, unidict, outseq_path, outlist_path)

    print(f"[INFO] Processing {len(swarm_lines)} OTUs in parallel...")
    with open(outseq_path, "w") as outseq, open(outlist_path,"w") as outlist:
        with ProcessPoolExecutor(max_workers=threads, initializer=init_worker, initargs=(readdict, unidict)) as executor:
            futures = [
                executor.submit(process_swarm_line, i + 1, line)
                for i, line in enumerate(swarm_lines)
            ]
            for future in tqdm(as_completed(futures), total=len(futures), desc="Processing OTUs"):
                seq_text, list_text = future.result()
                outseq.write(seq_text)
                outlist.write(list_text)

    print("[DONE] OTU reconstruction complete.")
    return  outseq_path, outlist_path

def main():
    if len(sys.argv) != 5:
        print("Usage: python3 postSwarm_step3_HPC.py <swarm_output> <derep_map> <derep_fasta> <threads>")
        sys.exit(1)

    process_postswarm(Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]), int(sys.argv[4]))
#     swarm_output_path = Path(sys.argv[1])
#     derep_map_path = Path(sys.argv[2])
#     derep_seq_path = Path(sys.argv[3])
#     thread = int(sys.argv[4])


if __name__ == "__main__":
    main()
