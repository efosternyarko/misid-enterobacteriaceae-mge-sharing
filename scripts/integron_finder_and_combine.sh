#!/bin/bash
# Run IntegronFinder on chromosome and plasmid sequences, then combine all
# per-isolate summaries into one TSV.
#
# FIX APPLIED: IntegronFinder nests its own results subdirectory inside
# whatever --outdir you give it (e.g. Result_Integron_Finder_<name>/), so
# the true isolate name is TWO directory levels above the .integrons file,
# not one. The original single-dirname approach
# (`isolate=$(basename "$(dirname "$summary")")`) returned IntegronFinder's
# internal subdirectory name instead of the isolate name. Using `find`
# with double-dirname fixes this.

# Run IntegronFinder on chromosome sequences
for chr_fasta in 03_MGE_sharing/01_mob_recon/*/chromosome.fasta; do
    isolate=$(basename "$(dirname "$chr_fasta")")
    integron_finder \
        --cpu 4 \
        --pdf \
        --outdir "03_MGE_sharing/05_integrons/$isolate/" \
        "$chr_fasta"
done

# Also run on plasmid sequences
for plasmid_fasta in 03_MGE_sharing/01_mob_recon/*/plasmid_*.fasta; do
    isolate=$(basename "$(dirname "$plasmid_fasta")")
    plasmid=$(basename "$plasmid_fasta" .fasta)
    integron_finder \
        --cpu 4 \
        --outdir "03_MGE_sharing/05_integrons/${isolate}_${plasmid}/" \
        "$plasmid_fasta"
done

# Combine all integron summaries -- double-dirname to skip past
# IntegronFinder's own results subdirectory and reach the isolate name
find 03_MGE_sharing/05_integrons -name "*.integrons" | while read -r summary; do
    isolate=$(basename "$(dirname "$(dirname "$summary")")")
    awk -v iso="$isolate" 'NR>1 {print iso "\t" $0}' "$summary"
done > 03_MGE_sharing/05_integrons/all_integrons_combined_clean.tsv

wc -l 03_MGE_sharing/05_integrons/all_integrons_combined_clean.tsv
