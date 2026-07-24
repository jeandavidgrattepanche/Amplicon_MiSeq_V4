#!/usr/bin/env python3

from Bio import SeqIO
import time
from pathlib import Path
from tqdm import tqdm
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
import argparse


def run_alignment(tool, ref_seq, target_id, target_seq):
    cmd = [
        tool,
        f"asis::{ref_seq}",
        f"asis::{target_seq}",
        "-gapopen=10", "-gapextend=0.5",
        "-auto", "-stdout"
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"[ERROR] {tool} alignment failed for {target_id}:\n{result.stderr}")
    return result.stdout


def parse_alignment_output(output):
    midline = ""
    flag = False
    for line in output.splitlines():
        if "#=======================================" in line:
            flag = True
            continue
        elif "#---------------------------------------" in line:
            flag = False
        elif flag and line.startswith(" "):
            midline += line.strip()

    gap = ident = sim = diff = 0
    midline += "*"
    midline2, id_count = "", 0
    for char in midline:
        if id_count < 4:
            if char == '|':
                id_count += 1
            else:
                id_count = 0
        else:
            midline2 += char

    for char in midline2:
        if char == ' ':
            gap += 1
        elif char == '|':
            ident += 1
            sim += 1
        elif char == ':':
            sim += 1
        elif char == '.':
            diff += 1

    aln_len = gap + sim + diff
    sim_ratio = sim / aln_len * 100 if aln_len > 0 else 0
    return aln_len, ident, sim, sim_ratio


def process_sequence(tool, ref_id, ref_seq, target_id, target_seq):
    output = run_alignment(tool, ref_seq, target_id, target_seq)
    aln_len, ident, sim, sim_ratio = parse_alignment_output(output)
    return target_id, target_seq, aln_len, ident, sim, sim_ratio


def process_contaminant_filter(seqfile, tool, threads ):
    t0 = time.time()

    outpath = seqfile.parent.parent
    chimera_dir = outpath / "chimeras"
    chimera_dir.mkdir(parents=True, exist_ok=True)

    seqdict = {rec.id: str(rec.seq) for rec in SeqIO.parse(seqfile, "fasta")}
    all_ids = list(seqdict.keys())
    max_id = max(all_ids, key=lambda x: int(x.split('_')[1].split('r')[0]))
    ref_seq = seqdict[max_id]

    with open(seqfile.with_name(seqfile.stem + "_nocont.fasta"), "w") as nocont_f, \
        open(seqfile.with_name(seqfile.stem + "_cont.fasta"), "w") as cont_f, \
        open(seqfile.with_name(seqfile.stem + "_pairwise_out_scores.csv"), "w") as score_f:
        nocont_f.write(f">{max_id}\n{seqdict[max_id]}\n")
        with ThreadPoolExecutor(max_workers=threads) as executor:
            futures = {
                executor.submit(process_sequence, tool, max_id, ref_seq, sid, seqdict[sid]): sid
                for sid in all_ids if sid != max_id
            }
            for future in tqdm(as_completed(futures), total=len(futures), desc=f"Aligning with {tool}"):
                sid, seq, aln_len, ident, sim, sim_ratio = future.result()
                score_f.write(f"{max_id},{sid},len: {aln_len},identity: {ident},similarity: {sim},similarity%: {sim_ratio:.2f}\n")

                if sim_ratio < 50:
                    cont_f.write(f">{sid}\n{seq}\n")
                else:
                    nocont_f.write(f">{sid}\n{seq}\n")

    print(f"[INFO] Done. Cleaned FASTA written to: {seqfile.stem}_nocont.fasta")
    elapsed = time.time() - t0
    print(f"[DONE] Finished in {elapsed:.2f} seconds.")

def main():
#     parser = argparse.ArgumentParser(description="Contaminant removal using pairwise alignment (water or needle)")
#     parser.add_argument("input_fasta", help="Input fasta file")
#     parser.add_argument("--tool", choices=["water", "needle"], default="water", help="Alignment tool to use")
#     parser.add_argument("--threads", type=int, default=8, help="Number of parallel threads")
#     args = parser.parse_args()
    print("Use process_contaminant_filter() through pipeline wrapper.")
    
if __name__ == "__main__":
    main()
