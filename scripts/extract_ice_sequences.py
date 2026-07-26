#!/usr/bin/env python3
# Extract standalone ICE nucleotide sequences for BLAST -- ICEfinder2 outputs
# per-ICE HTML and gene JSON files but not standalone FASTA, so these are
# sliced out of the merged chromosome using the coordinate ranges recorded
# in each isolate's _ICEsum.json. Requires biopython (icefinder2_env).

import json, glob, os
from Bio import SeqIO

ice_result_dir = os.path.expanduser("~/shared-team/ronnie.dir/ICEfinder2_linux/result")
recon_dir      = os.path.expanduser("~/shared-team/ronnie.dir/07_mge_sharing/01_mob_recon")
out_fasta      = os.path.expanduser("~/shared-team/ronnie.dir/07_mge_sharing/03_ICE_detection/all_ICE_sequences.fasta")

with open(out_fasta, "w") as out:
    for summary in glob.glob(f"{ice_result_dir}/**/*_ICEsum.json", recursive=True):
        isolate_base = os.path.splitext(os.path.basename(
            summary.replace("_ICEsum.json", "")))[0]
        # match isolate name to mob_recon directory (strip suffix cruft
        # added by the merge/ICEfinder2 steps, e.g. "_chromosome_merged")
        chr_fasta = None
        isolate_label = None
        for d in os.listdir(recon_dir):
            if isolate_base in d:
                candidate = os.path.join(recon_dir, d, "chromosome.fasta")
                if os.path.exists(candidate):
                    chr_fasta = candidate
                    isolate_label = d
                    break
        if not chr_fasta:
            print(f"WARNING: chromosome not found for {isolate_base}")
            continue
        chrom = list(SeqIO.parse(chr_fasta, "fasta"))
        chrom_seq = "".join(str(r.seq) for r in chrom)
        with open(summary) as f:
            ices = json.load(f)
        if not isinstance(ices, list):
            ices = [ices]
        for ice in ices:
            loc = ice.get("location", "")
            if not loc:
                continue
            start, end = (int(x) - 1 for x in loc.split(".."))
            ice_seq = chrom_seq[start:end]
            ice_id = ice.get("detail", ice.get("region", "ICE"))
            out.write(f">{isolate_label}__{ice_id}\n{ice_seq}\n")

print("Done. Sequences written to", out_fasta)
