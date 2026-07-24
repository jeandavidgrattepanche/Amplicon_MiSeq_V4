#!/usr/bin/env python3

"""
contaminant_filter_consensus_HPC.py

Consensus-based contaminant filtering for barcode amplicons.

Strategy:
---------
1. Select top-N most abundant OTUs
2. Compute pairwise similarities among top OTUs
3. Identify dominant coherent cluster
4. Choose medoid OTU as reference
5. Compare all OTUs to reference
6. Remove divergent sequences

Designed for:
- barcode amplicons
- off-target amplification detection
- mixed-domain contamination
- HPC execution

Dependencies:
-------------
vsearch
biopython
tqdm

Example:
--------
python contaminant_filter_consensus_HPC.py \
    input.fasta \
    --top_n 10 \
    --sim_cutoff 60 \
    --cluster_cutoff 85 \
    --threads 2
"""

import argparse
import time
import subprocess
from pathlib import Path
from Bio import SeqIO
from tqdm import tqdm
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed


# ============================================================
# Utilities
# ============================================================

def parse_abundance(seq_id):

    try:
        return int(seq_id.split("_")[1].split("r")[0])
    except:
        return 1


def run_vsearch_pairwise(seq1, seq2):

    cmd = [
        "vsearch",
        "--allpairs_global",
        "-",
        "--acceptall",
        "--threads", "1"
    ]

    fasta_input = f">A\n{seq1}\n>B\n{seq2}\n"

    result = subprocess.run(
        cmd,
        input=fasta_input,
        capture_output=True,
        text=True
    )

    return result.stdout


def compute_identity(seq1, seq2):

    cmd = [
        "vsearch",
        "--usearch_global",
        "-",
        "--db", "-",
        "--id", "0.0",
        "--blast6out", "/dev/stdout",
        "--acceptall",
        "--threads", "1"
    ]

    fasta_query = f">query\n{seq1}\n"
    fasta_db = f">target\n{seq2}\n"

    result = subprocess.run(
        cmd,
        input=fasta_query + fasta_db,
        capture_output=True,
        text=True
    )

    if result.stdout.strip():

        try:
            identity = float(result.stdout.split("\t")[2])
            return identity

        except:
            return 0.0

    return 0.0


# ============================================================
# Top OTU selection
# ============================================================

def load_sequences(fasta_file):

    seqdict = {}

    for record in SeqIO.parse(fasta_file, "fasta"):

        seqdict[record.id] = str(record.seq)

    return seqdict


def select_top_otus(seqdict, top_n):

    sorted_ids = sorted(
        seqdict.keys(),
        key=parse_abundance,
        reverse=True
    )

    return sorted_ids[:top_n]


# ============================================================
# Pairwise matrix
# ============================================================

def pairwise_worker(args):

    id1, seq1, id2, seq2 = args

    sim = compute_identity(seq1, seq2)

    return id1, id2, sim


def build_similarity_matrix(seqdict, top_ids, threads):

    tasks = []

    for i in range(len(top_ids)):

        for j in range(i + 1, len(top_ids)):

            id1 = top_ids[i]
            id2 = top_ids[j]

            tasks.append((
                id1,
                seqdict[id1],
                id2,
                seqdict[id2]
            ))

    simdict = defaultdict(dict)

    with ThreadPoolExecutor(max_workers=threads) as executor:

        futures = {
            executor.submit(pairwise_worker, t): t
            for t in tasks
        }

        for future in tqdm(
            as_completed(futures),
            total=len(futures),
            desc="Top OTU pairwise similarities"
        ):

            id1, id2, sim = future.result()

            simdict[id1][id2] = sim
            simdict[id2][id1] = sim

    return simdict


# ============================================================
# Dominant cluster
# ============================================================

def identify_consensus_cluster(simdict, top_ids, cluster_cutoff):

    cluster_sizes = {}

    for otu in top_ids:

        neighbors = [
            other
            for other, sim in simdict[otu].items()
            if sim >= cluster_cutoff
        ]

        cluster_sizes[otu] = len(neighbors)

    best_otu = max(cluster_sizes, key=cluster_sizes.get)

    cluster = [
        other
        for other, sim in simdict[best_otu].items()
        if sim >= cluster_cutoff
    ]

    cluster.append(best_otu)

    return list(set(cluster))


# ============================================================
# Medoid selection
# ============================================================

def choose_medoid(simdict, cluster):

    avg_similarity = {}

    for otu in cluster:

        sims = []

        for other in cluster:

            if otu == other:
                continue

            sims.append(simdict[otu].get(other, 0))

        avg_similarity[otu] = (sum(sims) / len(sims) if sims else 0)

    medoid = max(avg_similarity, key=avg_similarity.get)

    return medoid


# ============================================================
# Global filtering
# ============================================================

def filter_worker(args):

    otu_id, seq, ref_seq, sim_cutoff = args

    sim = compute_identity(seq, ref_seq)

    keep = sim >= sim_cutoff

    return otu_id, seq, sim, keep


def filter_sequences(seqdict, ref_id, sim_cutoff, threads):

    ref_seq = seqdict[ref_id]

    tasks = [
        (otu_id, seq, ref_seq, sim_cutoff)
        for otu_id, seq in seqdict.items()
    ]

    kept = []
    removed = []

    with ThreadPoolExecutor(max_workers=threads) as executor:

        futures = {
            executor.submit(filter_worker, t): t
            for t in tasks
        }

        for future in tqdm(
            as_completed(futures),
            total=len(futures),
            desc="Filtering sequences"
        ):

            otu_id, seq, sim, keep = future.result()

            if keep:
                kept.append((otu_id, seq, sim))
            else:
                removed.append((otu_id, seq, sim))

    return kept, removed


# ============================================================
# Output
# ============================================================

def write_output(output_file, records):

    with open(output_file, "w") as out:

        for otu_id, seq, sim in records:

            out.write(f">{otu_id};sim={sim:.2f}\n")
            out.write(f"{seq}\n")


# ============================================================
# Main
# ============================================================

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("input_fasta")

    parser.add_argument(
        "--top_n",
        type=int,
        default=10
    )

    parser.add_argument(
        "--sim_cutoff",
        type=float,
        default=60
    )

    parser.add_argument(
        "--cluster_cutoff",
        type=float,
        default=70
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=8
    )

    args = parser.parse_args()

    fasta_file = Path(args.input_fasta)

    t0 = time.time()

    output_keep = fasta_file.with_name(
        fasta_file.stem + "_nocont.fasta"
    )

    output_remove = fasta_file.with_name(
        fasta_file.stem + "_cont.fasta"
    )

    print("[INFO] Loading sequences...")
    seqdict = load_sequences(fasta_file)

    print("[INFO] Selecting top OTUs...")
    top_ids = select_top_otus(
        seqdict,
        args.top_n
    )

    print("[INFO] Building similarity matrix...")
    simdict = build_similarity_matrix(
        seqdict,
        top_ids,
        args.threads
    )
#     print(simdict)
    
    print("[INFO] Identifying dominant cluster...")
    cluster = identify_consensus_cluster(
        simdict,
        top_ids,
        args.cluster_cutoff
    )

    print(f"[INFO] Dominant cluster size: {len(cluster)}")

    print("[INFO] Selecting medoid OTU...")
    medoid = choose_medoid(simdict, cluster)

    print(f"[INFO] Medoid OTU: {medoid}")

    print("[INFO] Filtering sequences...")
    kept, removed = filter_sequences(
        seqdict,
        medoid,
        args.sim_cutoff,
        args.threads
    )

    print("[INFO] Writing outputs...")

    write_output(output_keep, kept)
    write_output(output_remove, removed)

    print(f"[INFO] Kept: {len(kept)}")
    print(f"[INFO] Removed: {len(removed)}")
    elapsed = time.time() - t0
    print(f"[DONE] Finished in {elapsed:.2f} seconds.")


if __name__ == "__main__":
    main()