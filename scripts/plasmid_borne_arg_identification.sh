#!/bin/bash
# Confirm plasmid-borne ARGs directly: run AMRFinderPlus on each isolate's
# concatenated plasmid FASTAs alone (from mob_recon), so any AMR hit found
# is unambiguously plasmid-located rather than inferred from contig-level
# molecule_type alone.

MOB_DIR="07_mge_sharing/01_mob_recon"
OUT_DIR="07_mge_sharing/08_plasmid_amr"
mkdir -p "$OUT_DIR"

for isolate_dir in "$MOB_DIR"/*/; do
    isolate=$(basename "$isolate_dir")
    plasmid_combined="${OUT_DIR}/${isolate}_plasmids_combined.fasta"
    cat "$isolate_dir"plasmid_*.fasta > "$plasmid_combined" 2>/dev/null
    if [ ! -s "$plasmid_combined" ]; then
        echo "No plasmid FASTAs for $isolate -- skipping"
        rm -f "$plasmid_combined"
        continue
    fi
    echo "AMRFinderPlus: $isolate"
    amrfinder \
        --nucleotide "$plasmid_combined" \
        --plus \
        --output "${OUT_DIR}/${isolate}_plasmid_amr.tsv"
done
echo "All done."
