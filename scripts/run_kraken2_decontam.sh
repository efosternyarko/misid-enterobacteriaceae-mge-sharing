#!/bin/bash

# Decontaminate assemblies using Kraken2 contig-level classification
# Usage: bash run_kraken2_decontam.sh
# Place this script in the same directory as Ronnie_genomes/
#
# Requires:
#   - kraken2 (conda install -c bioconda kraken2)
#   - krakentools (pip install krakentools)
#   - checkm2 (for post-filtering QC)
#   - A Kraken2 database (standard or k2_standard); update KRAKEN2_DB below

GENOME_DIR="Ronnie_genomes"
KRAKEN2_DB="/path/to/kraken2_db"    # <- UPDATE THIS
KRAKEN2_OUT_DIR="02_kraken2"
FILTERED_DIR="03_filtered_assemblies"
CHECKM2_OUT_DIR="04_qc_filtered"
THREADS=8

# Contaminated samples to attempt rescue (contamination 5-50%)
# Heavily contaminated samples excluded (C8 104%, R9 101%, R12 99%, K10 88%, K14 84%, K27 56%)
declare -A TARGET_TAXIDS=(
    ["C49"]="562"       # Escherichia coli
    ["R16"]="562"
    ["R2"]="562"
    ["R17"]="1463165"   # Klebsiella quasipneumoniae
    ["C25"]="548"       # Klebsiella aerogenes
)

mkdir -p "$KRAKEN2_OUT_DIR" "$FILTERED_DIR"

for sample in "${!TARGET_TAXIDS[@]}"; do
    # Find the assembly across species subdirectories
    assembly=$(find "$GENOME_DIR" -name "${sample}.fasta" | head -1)
    if [[ -z "$assembly" ]]; then
        echo "WARNING: no assembly found for $sample, skipping"
        continue
    fi

    target_taxid="${TARGET_TAXIDS[$sample]}"

    echo "Running Kraken2 on $sample (target taxid: $target_taxid)..."
    kraken2 --db "$KRAKEN2_DB" \
        --threads "$THREADS" \
        --output "${KRAKEN2_OUT_DIR}/${sample}.out" \
        --report "${KRAKEN2_OUT_DIR}/${sample}_report.txt" \
        "$assembly"

    echo "Extracting reads matching target taxid $target_taxid for $sample..."
    extract_kraken_reads.py \
        -k "${KRAKEN2_OUT_DIR}/${sample}.out" \
        -s "$assembly" \
        -o "${FILTERED_DIR}/${sample}.fasta" \
        -t "$target_taxid" \
        --include-children \
        --report "${KRAKEN2_OUT_DIR}/${sample}_report.txt"

    echo "Done: $sample"
done

echo "Running CheckM2 on decontaminated assemblies..."
checkm2 predict --threads "$THREADS" \
    --input "$FILTERED_DIR/" \
    --output-directory "$CHECKM2_OUT_DIR/" \
    --extension fasta

echo "All done. Review $CHECKM2_OUT_DIR/quality_report.tsv before proceeding."
