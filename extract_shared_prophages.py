#!/usr/bin/env python3
"""
extract_shared_prophages.py

Purpose:
--------
Extract integrated prophage regions from bacterial genomes based on
overlapping PhiSpy and Phigaro predictions. For each overlapping prophage,
the script extracts the union region from the corresponding genome FASTA file.

Inputs:
-------
1. metadata.tsv
   - Tab-separated file mapping biosample_id to genome_id
   - Required columns: genome_id, biosample_id

2. combined_overlap.tsv
   - Output file containing overlapping PhiSpy–Phigaro predictions
   - Must include coordinates for both tools

3. genomes_dir/
   - Directory containing genome FASTA files
   - Genome files must be named as: <genome_id>.fna

4. output_dir/
   - Directory where extracted prophage FASTA files will be written

Outputs:
--------
- One FASTA (.fna) file per extracted prophage region
- Duplicate regions are retained and suffixed with _dupN if needed

Usage:
------
python extract_prophages.py <metadata.tsv> <combined_overlap.tsv> <genomes_dir> <output_dir>

Example:
--------
python extract_prophages.py metadata.tsv combined_overlap.tsv genomes/ extracted_prophages/
Author: Aqib Javaid
"""

import os
import sys
import csv
from collections import defaultdict
from Bio import SeqIO

def main():
    if len(sys.argv) != 5:
        print("Usage: python extract_prophages.py <metadata.tsv> <combined_overlap.tsv> <genomes_dir> <output_dir>")
        sys.exit(1)

    metadata_file = sys.argv[1]
    combined_file = sys.argv[2]
    genomes_dir = sys.argv[3]
    output_dir = sys.argv[4]

    os.makedirs(output_dir, exist_ok=True)

    FLANK = 0  # 5 kb flanking sequence on each side

    # --- 1. Load genome_id ↔ biosample_id mapping ---
    biosample_to_genome = {}
    with open(metadata_file) as meta_fh:
        reader = csv.DictReader(meta_fh, delimiter="\t")
        for row in reader:
            genome_id = row["genome_id"].strip()
            biosample_id = row["biosample_id"].strip()
            biosample_to_genome[biosample_id] = genome_id

    # --- 2. Load prophage predictions and compute union ---
    prophages = []
    with open(combined_file) as comb_fh:
        reader = csv.DictReader(comb_fh, delimiter="\t")
        for row in reader:
            phigaro_contig = row["phigaro_contig"].strip()
            phigaro_id = row["phigaro_id"].strip()
            biosample_id = phigaro_id.split("_")[0]
            genome_id = biosample_to_genome.get(biosample_id)

            if not genome_id:
                print(f"[WARNING] No genome_id found for biosample {biosample_id}")
                continue

            try:
                start = min(int(row["phispy_start"]), int(row["phigaro_start"]))
                end = max(int(row["phispy_end"]), int(row["phigaro_end"]))
            except ValueError:
                print(f"[WARNING] Invalid coordinates for {phigaro_id}, skipping row")
                continue

            # Treat start=0 as 1
            start = max(1, start)

            prophages.append({
                "contig": phigaro_contig,
                "genome_id": genome_id,
                "start": start,
                "end": end,
                "phigaro_id": phigaro_id
            })

    # --- 3. Group prophages by genome ---
    prophages_by_genome = defaultdict(list)
    for p in prophages:
        prophages_by_genome[p["genome_id"]].append(p)

    # --- 4. Extraction with duplicates and flanks ---
    duplicate_count = defaultdict(int)

    for genome_id, ph_list in prophages_by_genome.items():
        genome_path = os.path.join(genomes_dir, f"{genome_id}.fna")
        if not os.path.exists(genome_path):
            print(f"[WARNING] Genome file not found: {genome_path}")
            continue

        contigs = SeqIO.to_dict(SeqIO.parse(genome_path, "fasta"))

        for p in ph_list:
            if p["contig"] not in contigs:
                print(f"[WARNING] Contig {p['contig']} not found in {genome_id}")
                continue

            contig_len = len(contigs[p["contig"]].seq)

            # Clamp coordinates and include ±5kb flanking regions
            flank_start = max(1, p["start"] - FLANK)
            flank_end = min(contig_len, p["end"] + FLANK)

            if flank_start > flank_end:
                print(f"[WARNING] Adjusted start > end for {p['phigaro_id']}, skipping")
                continue

            key = (p["contig"], flank_start, flank_end)
            duplicate_count[key] += 1
            dup_suffix = "" if duplicate_count[key] == 1 else f"_dup{duplicate_count[key]-1}"

            out_path = os.path.join(output_dir, f"{p['phigaro_id']}{dup_suffix}.fna")
            subseq = contigs[p["contig"]].seq[flank_start-1:flank_end]

            with open(out_path, "w") as out_fh:
                out_fh.write(f">{p['phigaro_id']}|{p['contig']}|{flank_start}-{flank_end}{dup_suffix}\n")
                for i in range(0, len(subseq), 80):
                    out_fh.write(str(subseq[i:i+80]) + "\n")

            print(f"[INFO] Extracted {p['phigaro_id']}{dup_suffix} → {out_path}")


if __name__ == "__main__":
    main()
