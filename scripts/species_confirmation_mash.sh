#!/bin/bash
# Species confirmation via Mash distance against the RefSeq sketch database.
# Run from the project root with the misid_ecoli conda environment active.

conda activate misid_ecoli

# Sketch all assemblies
mash sketch -o 02_species/sketch 00_assemblies/*fasta
cd 02_species/

# Download RefSeq genome sketch (once)
wget https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh

# Use mash dist to identify top 5 species hit for each assembly
mash dist refseq.genomes.k21s1000.msh sketch.msh \
    | sort -k3 -n | awk 'seen[$2]++ < 5' | sort -k2,2 -k3,3n > mash_results.txt

# Extract unique accession numbers and query species using ncbi datasets
awk '{print $1}' mash_results.txt | grep -oE 'GCF_[0-9]+\.[0-9]+' | sort -u > unique_accessions.txt
datasets summary genome accession --inputfile unique_accessions.txt --as-json-lines \
    | jq -r '[.accession, .organism.organism_name] | @tsv' > species_summary.txt

# Left join to remediate de-duplication of accession numbers performed by
# ncbi datasets, and merge back with the mash distance output
awk -F'\t' '
NR==FNR { species[$1]=$2; next }
{
    match($1, /GCF_[0-9]+\.[0-9]+/);
    acc = substr($1, RSTART, RLENGTH);
    print $2 "\t" acc "\t" species[acc] "\t" $3
}' species_summary.txt mash_results.txt > final_merged_results.txt

echo "Done. Review final_merged_results.txt, then upload KPSC assemblies to"
echo "Pathogenwatch and run Speciator for disambiguation, flagging any discrepancies."
