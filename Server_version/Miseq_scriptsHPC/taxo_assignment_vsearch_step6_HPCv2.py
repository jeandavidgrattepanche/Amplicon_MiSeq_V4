#!/usr/bin/env python3

from Bio import SeqIO
from pathlib import Path
from tqdm import tqdm
import subprocess
import sys
import os
blastdict = {}
taxdict = {}

def run_vsearch(query_file, output_path, SSU_DB, idmin, qcov, maxacc, threads):
    cmd = [
        "vsearch",
        "--usearch_global", str(query_file),
        "--db", str(SSU_DB),
        "--strand", "both", #plus in v1
        "--id", str(idmin / 100),
        "--threads", str(threads),
        "--query_cov", str(qcov / 100),
        "--maxaccepts", str(maxacc),
        "--userout", str(output_path / "VsearchBLAST.tsv"),
        "--userfields query+target+id+alnlen+mism+gaps+qcov+tcov+bits",
    ]
    print("[INFO] Running VSEARCH:", " ".join(map(str, cmd)))
    subprocess.run(cmd, check=True)

blastdict.clear()
taxdict.clear()

def load_db(SSU_DB):
    for ref in SeqIO.parse(SSU_DB, "fasta"):
        taxdict[ref.id] = ref.description.split("\t")[0].replace(" ", "_")

def parse_vsearch_output(blast_file):
    with open(blast_file, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue  # skip malformed lines
            try:
                query, subject, pid = parts[0], parts[1], float(parts[2])
            except ValueError:
                print(f"[WARNING] Skipping malformed line: {line.strip()}")
                continue  # skip lines with invalid float conversion

            tax = taxdict.get(subject, "unknown")
            new_entry = f"{tax}\t{line.strip()}"
            if query not in blastdict:
                blastdict[query] = new_entry
            if query in blastdict:
#                 print(new_entry,'\n\n',pid, '\t', blastdict[query])
                if pid > float(blastdict[query].split('\t')[3]):
                    blastdict[query] = new_entry

def write_results(output_path, ngs_reduced_file, taxa):
    tokeep = []

    all_query_ids = set()
    for rec in SeqIO.parse(ngs_reduced_file, "fasta"):
        all_query_ids.add(rec.id)

    clean_path = output_path / "VsearchBLAST_clean.tsv"
    group_path = output_path / "VsearchBLAST_inGroup.tsv"
    fasta_out = output_path / "Seq_reads_inGroup.fasta"

    with open(clean_path, "w") as out_clean, open(group_path, "w") as out_group:
        for query in tqdm(sorted(all_query_ids), desc="Writing BLAST TSV results"):
            if query in blastdict:
                result = blastdict[query]
                out_clean.write(result + "\n")
                if taxa is None or taxa in result:
                    out_group.write(result + "\n")
                    tokeep.append(query.split('_')[0])
            else:
                out_clean.write(f"NO_HIT\t{query}\tNO_HITinSILVA\n")

    with open(fasta_out, "w") as out_seq:
        for rec in tqdm(SeqIO.parse(ngs_reduced_file, "fasta"), desc="Writing ingroup FASTA"):
            if rec.id in blastdict and (taxa is None or taxa in blastdict[rec.id]):
                out_seq.write(f">{rec.description}\n{rec.seq}\n")

    return tokeep
    

def write_map_subset(otu_map_file, tokeep, output_path):
    tokeep_set = set(tokeep)  # Ensure fast lookup
    out_file = output_path / "Seq_map_inGroup.txt"

    with open(otu_map_file) as infile, open(out_file, "w") as outmap:
        kept = 0

        # Count total lines first for accurate progress bar (optional but nice)
        with open(otu_map_file) as f:
            total_lines = sum(1 for _ in f)

        for line in tqdm(infile, total=total_lines, desc="Filtering OTU map"):
            parts = line.split('\t', 1)
            if parts[0] in tokeep_set:
                outmap.write(line)
                kept += 1

        print(f"[INFO] {kept}/{total_lines} lines kept ({(kept / total_lines) * 100:.1f}%)")


def filter_input_fasta(ngs_file, read_cutoff):
    reduced_file = Path(ngs_file).with_name(Path(ngs_file).stem + "_reduced.fas")
    with open(reduced_file, "w") as out:
        for seq in SeqIO.parse(ngs_file, "fasta"):
            try:
                if int(seq.description.split("_")[1].replace('r', '')) > int(read_cutoff):
                    out.write(f">{seq.description}\n{seq.seq}\n")
            except (IndexError, ValueError):
                continue
    return reduced_file

def process_taxonomic_assignment(ngs_file,otu_map,db_path,taxa,idmin,qcov,maxacc, readcutoff,threads):
    blastdict.clear()
    taxdict.clear()

    output_path = Path(ngs_file).parent.parent / "taxonomic_assignment"
    output_path.mkdir(parents=True, exist_ok=True)

    load_db(db_path)

    reduced = filter_input_fasta(ngs_file, readcutoff)
    run_vsearch(reduced, output_path,db_path, idmin, qcov, maxacc, threads)
    parse_vsearch_output(output_path / "VsearchBLAST.tsv")
    tokeep = write_results(output_path, reduced, taxa)
    write_map_subset(otu_map, tokeep, output_path)


if __name__ == "__main__":
    print("Use process_taxonomic_assignment() through pipeline wrapper.")