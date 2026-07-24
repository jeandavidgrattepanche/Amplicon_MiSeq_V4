#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MiSeq Pipeline (Part 1 of 2)
Optimized for PSC
Author: Jean-David Grattepanche
Refactor & Optimization: ChatGPT
"need:
    - module load anaconda3
    - conda activate Amplicon" on PSC bridges2
    - module load bbmap
    
This script performs:
- Reads the sample list
- Pools and renames raw FASTA reads from multiple samples
- Prepares data for downstream clustering and chimera detection

===========================
🔧 Requirements to run:
===========================
Python 3.11+ with the following modules installed:
    biopython
    tqdm
    pandas


===========================
🛠️ External tools in PATH:
===========================
    BBTools: bbmap / bbduk / bbmerge
    vsearch 
    swarm 

===========================
💡 Usage example:
===========================
python MiSeq_pipeline_part1_PSC.py sample_list.txt RawDataFolder

This will create pooled and renamed FASTA files inside a folder named 'output/temp' (or similar),
ready for processing in Part 2 of the pipeline.
"""



import os
import sys
import subprocess
import argparse
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
from Bio import SeqIO
from Bio.Seq import Seq
import multiprocessing
print(f"SLURM CPUs: {os.environ.get('SLURM_CPUS_PER_TASK', 'NA')}")


def run_command(cmd, desc=None):
    if desc:
        print(f"[STEP] {desc}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"[ERROR] Command failed: {cmd}\n{result.stderr}")
    return result

def clean_sequence(seq):
    seq = str(seq)
    # Trim Ns from both ends
    start = 0
    end = len(seq)
    while start < end and seq[start] == 'N':
        start += 1
    while end > start and seq[end - 1] == 'N':
        end -= 1
    return seq[start:end]

def process_fasta(input_fasta):
    input_path = Path(input_fasta)
    base_name = input_path.stem.replace('.fasta', '')
    out_noN = input_path.with_name(f"{base_name}_noN.fasta")
#     out_withN = input_path.with_name(f"{base_name}_withNv2.fas")

    noN_seqs = []
    withN_seqs = []
    with open(input_path, "r") as handle:
        records = SeqIO.parse(handle, "fasta")
        for record in records:
            new_seq = Seq(clean_sequence(record.seq))
            if len(new_seq) > new_seq.count('N'):
                record.seq = new_seq
                if 'N' in new_seq:
                    withN_seqs.append(record)
                else:
                    noN_seqs.append(record)
            else:
                print(f"[SKIP] Only N: {record.id}")

    if noN_seqs:
        with open(out_noN, "w") as out1:
            SeqIO.write(noN_seqs, out1, "fasta")
    return out_noN
#     if withN_seqs:
#         with open(out_withN, "w") as out2:
#             SeqIO.write(withN_seqs, out2, "fasta")


def process_sample(sample, pathA, pathB, pathC, outputpath, minlen, tool_threads):
    rawfiles = list(Path(pathA).glob(f"{sample}*"))
    r1 = next((f for f in rawfiles if "R1" in f.name), None)
    r2 = next((f for f in rawfiles if "R2" in f.name), None)

    if not r1 or not r2:
        print(f"[WARN] Missing R1 or R2 for sample {sample}")
        return sample, False

    merged_fq = Path(pathB) / f"{sample}_merged.fastq"
    filtered_fa = Path(pathB) / f"{sample}_clean1.fasta"
    merged_noN = Path(pathB) / f"{sample}_clean1_noN.fasta"
    cleaned_fa = Path(pathB) / f"{sample}_clean2.fasta"
    trimmed_fa = Path(pathC) / f"{sample}_trimmed.fasta"
    swarm_out = Path(pathB) / f"{sample}_swarm.txt"
    stats_out = Path(pathB) / f"{sample}_swarm_stats.txt"
    if swarm_out.exists():
        print(f"[SKIP] {sample} already completed")
        return sample, True

    # Merge reads using BBMerge (assumes bbmerge.sh is in PATH)
    cmd_merge = f"bbmerge-auto.sh -Xmx4g threads={tool_threads} ecct extend2=150 loose iterations=5 in1={r1} in2={r2} out={merged_fq}"
    run_command(cmd_merge, desc=f"Merging reads with BBMerge for {sample}")

    # filter and dereplicate 1 for older vsearch --derep_fulllength  and --output for newer vsearch (>v2.20) --fastx_uniques and --fastaout
    cmd_filter = f"vsearch --fastx_uniques {merged_fq} --sizeout --fasta_width 0 --fastaout {filtered_fa} --threads {tool_threads}"
    run_command(cmd_filter, desc=f"Filtering and dereplicating pass1 {sample}")

	# remove N
#     cmd = f"python $PROJECT/Amplicon/Miseq_scripts/ext_remove_N_in_seqfile_M4max.py {filtered_fa}"
#     run_command(cmd)
    merged_noN = process_fasta(filtered_fa) 

    # filter and dereplicate 2
    cmd_filter2 = f"vsearch --derep_fulllength {merged_noN} --sizein --sizeout --fasta_width 0 --output {cleaned_fa} --threads {tool_threads}"
    run_command(cmd_filter2, desc=f"Filtering and dereplicating pass2 {sample}")
	

    # Clean with BBDuk (assumes bbduk.sh is in PATH)
    cmd_clean = f"bbduk.sh -Xmx4g threads={tool_threads} in={cleaned_fa} out={trimmed_fa} maq=20 minlen={minlen}"
    run_command(cmd_clean, desc=f"Cleaning reads with BBDuk for {sample}")


    # SWARM
    cmd_swarm = f"swarm -f -t {tool_threads} -s {stats_out} -d 1 -z {trimmed_fa} > {swarm_out}"
    run_command(cmd_swarm, desc=f"Clustering OTUs with SWARM for {sample}")

    return sample, True


def main():
    parser = argparse.ArgumentParser(description="MiSeq Pipeline (Part 1 of 2) Optimized for PSC using BBTools")
    parser.add_argument("samplefile", help="Tab-delimited sample list file")
    parser.add_argument("rawdata_folder", help="Raw data folder containing R1 and R2 files")
    parser.add_argument("--minlen", type=int, default=400, help="Minimum length for reads (default: 400)")
    parser.add_argument("--sample_threads", type=int, default=4)
    parser.add_argument("--tool_threads", type=int, default=2)
    args = parser.parse_args()

    pathA = Path(args.rawdata_folder)
    base = Path(args.samplefile).stem
    outputpath = pathA.parent / f"{base}_output"
    pathB = outputpath / "temp"
    pathC = outputpath / "RDP2"
    outputpath.mkdir(parents=True, exist_ok=True)
    pathB.mkdir(parents=True, exist_ok=True)
    pathC.mkdir(parents=True, exist_ok=True)

    with open(args.samplefile) as f:
        listsamp = [line.strip().split("\t")[0] for line in f if line.strip() and not line.startswith("MiSeqcode")]

    print(f"[INFO] Running on {len(listsamp)} samples with {args.sample_threads} samples at a time using {args.tool_threads} threads each")

    results = []
    with ProcessPoolExecutor(max_workers=args.sample_threads) as executor:
        futures = {executor.submit(process_sample, s, pathA, pathB, pathC, outputpath, args.minlen, args.tool_threads): s for s in listsamp}
        for future in tqdm(as_completed(futures), total=len(futures), desc="Pipeline Progress"):
            sample, success = future.result()
            results.append((sample, success))

    # Summary
    completed = sum(1 for _, success in results if success)
    failed = len(results) - completed
    print(f"\n[SUMMARY] Completed: {completed}, Failed: {failed} out of {len(results)} samples")


if __name__ == "__main__":
    main()
