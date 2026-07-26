#!/bin/bash
# MLST typing (species-specific schemes) and Bakta genome annotation.
# Run from the project root with the misid_ecoli conda environment active.

conda activate misid_ecoli

# Rename assembly files with organism prefix (once only)
cd 00_assemblies
for file in {R18,R8,R5,K59,K34,E66,E55,E47,E17,C56,C55,C43,E18,C52,C11}.fasta; do
    mv "$file" "EC_$file"
done
for file in {K31,K28,K4,A02,C6,C22}.fasta; do
    mv "$file" "KQ_$file"
done
for file in {E59,E60}.fasta; do
    mv "$file" "KV_$file"
done
cd ..

# E. coli -- Achtman scheme
mlst --scheme ecoli_achtman_4 00_assemblies/EC_*.fasta > 04_mlst/ecoli_mlst_achtman.tsv

# Klebsiella complex
mlst --scheme klebsiella 00_assemblies/KQ_*.fasta 00_assemblies/KV_*.fasta > 04_mlst/kleb_mlst.tsv

# Genome annotation with Bakta
for f in 00_assemblies/*.fasta; do
    name=$(basename "${f%.fasta}")
    bakta --db /path/to/DB --prefix "$name" \
        --output "03_annotation/${name}" "$f"
done
