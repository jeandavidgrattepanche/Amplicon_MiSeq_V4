#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MiSeq Prokaryotic Pipeline (Part 2 of 3)
Optimized for PSC bridges2
Author: Jean-David Grattepanche
Refactor & Optimization: ChatGPT

need:
    - module load anaconda3
    - conda activate Amplicon 
    - sbatch run_miseq_p2.sh 

This script performs:
- Renaming and pooling FASTA reads
- OTU clustering with SWARM
- Chimera detection with VSEARCH
- Optional taxonomy assignment with VSEARCH BLAST
- Outputs OTU tables and filtered FASTA files

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
    vsearch
    swarm
    EMBOSS water or needle (for contaminant removal)

===========================
💡 Usage example:
===========================
python MiSeq_pipeline_part2_PSC.py sample_list.txt RawDataFolder \\
    --swarm_d 1 --read_cutoff 100 --taxa [E or P] --assign_taxonomy --idmin 97 --qcov 80 --tool water

Output files will be written in a folder named '<samplelist>_output'.
"""


import os
import subprocess
import argparse
from pathlib import Path
from tqdm import tqdm
import multiprocessing


def get_parser():
    parser = argparse.ArgumentParser(
        description=(
            "MiSeq Pipeline (Part 2 of 3)\n"
            "Optimized for PSC "
            "This script processes dereplicated sequences, clusters OTUs, "
            "removes chimeras, assigns taxonomy, and creates OTU tables."
        ),
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=(
            "Example:\n"
            "  python MiSeq_pipeline_part2_PSC.py sample_list.txt RawDataFolder \\\n"
            "         --swarm_d 1 --read_cutoff 100 --assign_taxonomy --idmin 97 --qcov 80 --tool water\n\n"
            "Output files will be written in a folder called '<samplelist>_output' created next to your input."
        )
    )

    parser.add_argument("samplefile", help="Path to the sample list text file (one sample per line, tab-delimited)")
    parser.add_argument("rawdata_folder", help="Path to the folder containing raw FASTA files")
    parser.add_argument("--swarm_d", type=int, default=1, help="SWARM clustering distance (default: 1)")
    parser.add_argument("--read_cutoff", type=int, default=1, help="Minimum number of reads per OTU to keep (default: 1)")
    parser.add_argument("--assign_taxonomy", action="store_true", help="Enable taxonomy assignment using VSEARCH BLAST (default: off)")
    parser.add_argument("--idmin", type=float, default=0.0, help="Minimum BLAST identity percentage (default: 0.0)")
    parser.add_argument("--qcov", type=float, default=70.0, help="Minimum BLAST query coverage percentage (default: 70.0)")
    parser.add_argument("--maxacc", type=int, default=1, help="Maximum accepted BLAST result (default: 1)")
    parser.add_argument("--tool", choices=["water", "needle"], default="water", help="Tool to use for decontamination (default: water)")
    parser.add_argument("--threads", type=int, default=multiprocessing.cpu_count(), help="Number of threads to use (default: 10)" )
    parser.add_argument("--taxa", choices=["P", "E"], required=True, help="Indicate if your target taxa are [E]ukaryotes or [P]rokaryotes")
    parser.add_argument("--filter_taxa", action="store_true", help="Enable filtering by taxonomy using --taxa: if E then P are removed, if P then E are removed (default: off)" )
    parser.add_argument("--force_pool", action="store_true")
    parser.add_argument("--force_swarm", action="store_true")
    parser.add_argument("--force_chimera", action="store_true")
    parser.add_argument("--force_contaminant", action="store_true")
    parser.add_argument("--force_taxonomy", action="store_true")
    return parser


def run_command(cmd, desc=None):
    if desc:
        print(f"[INFO] Running: {desc}", flush=True)
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.stdout:
        print(result.stdout)
    if result.returncode != 0:
        if result.stderr:
            print(result.stderr)
        raise RuntimeError(f"[ERROR] Command failed:\n{cmd}")
    return result

def makesinglefastafile(file, pathC, listsample):
    from Miseq_scriptsHPC.pool_rename_step1_HPC import process_pool_rename
    process_pool_rename(Path(pathC)/file, listsample)

def parallel_makesinglefastafile(pathC, listsample):
    files = [f for f in os.listdir(pathC) if f.endswith(".fasta")]
    for f in tqdm(files, desc="Renaming FASTA files", unit="file"):
        makesinglefastafile(f, pathC, listsample)

def pick_otu_swarm(dSWARM, pathC, outputpath, listsample, readcutoff, threads):
    otu_dir = Path(outputpath) / "OTUs"
    otu_dir.mkdir(parents=True, exist_ok=True)

    print(f"[INFO] Pick OTUs with SWARM using {threads} threads")
    run_command(
        f"vsearch --derep_fulllength {pathC}/readpooled.fas --sizein --sizeout --strand both "
        f"--fasta_width 0 --output {otu_dir}/dereplicated_seqfile.fas --uc {otu_dir}/dereplicated_seqfile.map.txt --threads {threads}",
        desc="Dereplicating sequences with VSEARCH"
    )

    run_command(f"swarm -f -t {threads} -s {otu_dir}/statSWARM -d {dSWARM} -z {otu_dir}/dereplicated_seqfile.fas > {otu_dir}/derepseqfile_output.swarm", desc="Clustering OTUs with SWARM"  )
    from Miseq_scriptsHPC.postSwarm_step3_HPC import process_postswarm
    process_postswarm(Path(otu_dir)/"derepseqfile_output.swarm", Path(otu_dir)/"dereplicated_seqfile.map.txt", Path(otu_dir)/"dereplicated_seqfile.fas", threads)
    from Miseq_scriptsHPC.Add_numbers_step4_HPC import process_numbers
    process_numbers(Path(otu_dir)/"SWARM_postout.fas", Path(otu_dir)/"SWARM_postout.txt", listsample, threads)

def clean_otus(outputpath, readcutoff, threads):
    print("[INFO] Preparing files for Chimera check")
    otu_dir = Path(outputpath) / "OTUs"
    chimera_dir = Path(outputpath) / "chimeras"
    os.makedirs(chimera_dir, exist_ok=True)
 
    from Miseq_scriptsHPC.cleaning_step5_HPC import pre_uchime
    pre_uchime(Path(chimera_dir), Path(otu_dir)/"SWARM_postout_nosingleton.fas", Path(otu_dir)/"OTUtable_temp.txt", readcutoff)

    chimera_cmd = f"vsearch --uchime3_denovo {chimera_dir}/Seq_reads_test.fas --nonchimera {chimera_dir}/Seq_reads_nochimera_nosingleton.fas --uchimeout {chimera_dir}/chimeratable.txt --threads {threads}"
    run_command(chimera_cmd, desc="Removing chimeras with VSEARCH")
 
    from Miseq_scriptsHPC.cleaning_step5_HPC import post_uchime
    post_uchime(Path(chimera_dir)/"Seq_reads_test.fas")

def remove_contaminants(outputpath, tool, threads):
    chimera_dir = Path(outputpath) / "chimeras"
    from Miseq_scriptsHPC.contaminant_removal_step5b_HPC import process_contaminant_filter
    process_contaminant_filter(Path(chimera_dir)/"Seq_reads_nochimera_nosingleton_renamed.fas", tool, threads)


def run_blast(AssTaxo, FilTaxa, taxa, outputpath, idmin, qcov, maxacc, readcutoff, threads):
    if not AssTaxo:
        print("[INFO] Skipping taxonomic assignment")
        return

    input_fasta = (Path(outputpath)/ "chimeras"/ "Seq_reads_nochimera_nosingleton_renamed_nocont.fasta")
    otu_map = (Path(outputpath)/ "OTUs"/ "SWARM_postout.txt")

    print("[INFO] Running BLAST for taxonomic assignment")
    if taxa == "P":
        db_path = (Path("db_v4")/ "SILVA_138.2_backbone_clean.fasta") 
#         db_path = (Path("db_v4")/ "SILVA_138.2_SSURef_NR99_tax_silva_filtered.fasta") 
        taxa_filter = "Bacteria" if FilTaxa else None
    elif taxa == "E":
        db_path = (Path("db_v4")/ "pr2_version_5.0.0_SSU_UTAX.fasta") 
        taxa_filter = "Eukaryota" if FilTaxa else None
    else:
        raise ValueError(f"Unknown taxa mode: {taxa}")
    print(f"[INFO] Starting taxonomy assignment")
    from Miseq_scriptsHPC.taxo_assignment_vsearch_step6_HPCv2 import process_taxonomic_assignment
    process_taxonomic_assignment(input_fasta,otu_map,db_path,taxa_filter,float(idmin),float(qcov),int(readcutoff),threads)    
    print(f"[INFO] Completed taxonomy assignment")
    
def checkpoint_ok(path):
    return path.exists() and path.stat().st_size > 0

def main():
    parser = get_parser()
    args = parser.parse_args()

    pathA = Path(args.rawdata_folder)
    base = Path(args.samplefile).stem
    outputpath = pathA.parent / f"{base}_output"
    pathC = outputpath / "RDP2"

    pooled = pathC / "readpooled.fas"
    if checkpoint_ok(pooled) and not args.force_pool:
        print("[INFO] Found pooled FASTA, skipping pooling step")
    else:
        parallel_makesinglefastafile(pathC, args.samplefile)

    swarm_out = Path(outputpath) / "OTUs" / "SWARM_postout.fas"
    if checkpoint_ok(swarm_out) and not args.force_swarm:
        print("[INFO] Existing SWARM output detected, skipping clustering")
    else:
        pick_otu_swarm(args.swarm_d, pathC, outputpath, args.samplefile, args.read_cutoff, args.threads)

    nochim = Path(outputpath)/ "chimeras"/ "Seq_reads_nochimera_nosingleton.fas"
    if checkpoint_ok(nochim) and not args.force_chimera:
        print("[INFO] Chimera filtering already completed")
    else:
        clean_otus(outputpath, args.read_cutoff, args.threads)

    nocont = Path(outputpath)/ "chimeras"/ "Seq_reads_nochimera_nosingleton_renamed_nocont.fasta"
    if checkpoint_ok(nocont) and not args.force_contaminant:
        print("[INFO] contaminant filtering already completed")
    else:
        remove_contaminants(outputpath, args.tool, args.threads)

    tax_out = Path(outputpath)/ "taxonomic_assignment"/ "VsearchBLAST_clean.tsv"
    
    if checkpoint_ok(tax_out) and not args.force_taxonomy:
        print("[INFO] Taxonomy assignment already completed")
    else:
        run_blast(args.assign_taxonomy, args.filter_taxa, args.taxa, outputpath, args.idmin, args.qcov, args.maxacc, args.read_cutoff, args.threads)


if __name__ == "__main__":
    main()
