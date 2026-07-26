#!/bin/bash
# Clermont phylogroup typing (ezclermont) and virulence factor screening
# (ABRicate, Ecoli_VF database) for the E. coli isolates.

mkdir -p 08_virulence/clermont
for file in 00_assemblies/05_ecoli/*.fasta; do
    base_name=$(basename "$file" .fasta)
    echo "Processing $base_name..."
    ezclermont "$file" > "08_virulence/clermont/${base_name}_clermont.txt"
done

echo -e "Isolate\tPhylotype" > 08_virulence/clermont/summary_phylotypes.txt
for file in 08_virulence/clermont/EC_*_clermont.txt; do
    base_name=$(basename "$file" _clermont.txt)
    phylotype=$(cat "$file")
    echo -e "${base_name}\t${phylotype}" >> 08_virulence/clermont/summary_phylotypes.txt
done

# ABRicate screen using Ecoli_VF
abricate --db ecoli_vf 00_assemblies/*.fasta > 08_virulence/vfdb_all.tsv
