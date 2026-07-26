#!/bin/bash
# AMR gene detection (organism-specific AMRFinderPlus) and plasmid typing
# (MOB-suite mob_recon). Run from the project root, misid_ecoli environment.

# --- AMR gene detection ---
mkdir -p metadata
# metadata/species_organism_map.tsv is a manually curated two-column TSV:
# <isolate_id>\t<organism>  (organism must be one AMRFinderPlus recognises,
# e.g. Klebsiella_pneumoniae for all KPSC isolates, Escherichia for E. coli)

amrfinder -u   # download/update the AMRFinderPlus database

while IFS=$'\t' read -r isolate organism; do
    amrfinder -n "00_assemblies/${isolate}.fasta" \
        --organism "$organism" \
        -o "05_amr/${isolate}_amr.tsv"
done < metadata/species_organism_map.tsv

# --- Plasmid typing with MOB-suite ---
for f in 00_assemblies/*.fasta; do
    name=$(basename "${f%.fasta}")
    mob_recon --infile "$f" --outdir "06_plasmids/${name}"
done
