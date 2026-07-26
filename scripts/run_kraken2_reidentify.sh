#!/bin/bash

# Re-identify and re-filter assemblies after discovering the initial Kraken2
# decontamination pass targeted the WRONG dominant organism for several
# samples (they were misidentified in the original species label, not just
# contaminated). Usage: bash run_kraken2_reidentify.sh
# Run only after run_kraken2_decontam.sh and after inspecting its Kraken2
# reports to confirm the true dominant taxon per sample.

GENOME_DIR="Ronnie_genomes"
KRAKEN2_DB="/path/to/kraken2_db"    # <- UPDATE THIS
KRAKEN2_OUT_DIR="02_kraken2"
FILTERED_DIR="03_filtered_assemblies"
CHECKM2_OUT_DIR="04_qc_filtered"
THREADS=8

# Kraken2 reports show these samples are MISIDENTIFIED, not just contaminated.
# Target taxids updated to the DOMINANT organism in each assembly.
# Previous labels -> True dominant organism:
#   C49  E. coli          -> Klebsiella (85% K. pneumoniae)
#   R16  E. coli          -> Klebsiella (88% K. pneumoniae)
#   R2   E. coli          -> Klebsiella (94% K. variicola)
#   R17  K. quasipneumoniae -> Proteus mirabilis (88%)
#   C25  K. aerogenes     -> Klebsiella (80% K. pneumoniae dominant)
declare -A TARGET_TAXIDS=(
    ["C49"]="570"    # Klebsiella genus (--include-children captures K. pneumoniae + all strains)
    ["R16"]="570"    # Klebsiella genus
    ["R2"]="570"     # Klebsiella genus (K. variicola dominant)
    ["R17"]="583"    # Proteus genus (--include-children captures P. mirabilis)
    ["C25"]="570"    # Klebsiella genus
)

mkdir -p "$KRAKEN2_OUT_DIR" "$FILTERED_DIR"

for sample in "${!TARGET_TAXIDS[@]}"; do
    assembly=$(find "$GENOME_DIR" -name "${sample}.fasta" | head -1)
    if [[ -z "$assembly" ]]; then
        echo "WARNING: no assembly found for $sample, skipping"
        continue
    fi

    target_taxid="${TARGET_TAXIDS[$sample]}"

    echo "Re-extracting $sample against corrected target taxid $target_taxid..."
    extract_kraken_reads.py \
        -k "${KRAKEN2_OUT_DIR}/${sample}.out" \
        -s "$assembly" \
        -o "${FILTERED_DIR}/${sample}.fasta" \
        -t "$target_taxid" \
        --include-children \
        --report "${KRAKEN2_OUT_DIR}/${sample}_report.txt"

    echo "Done: $sample"
done

echo "Running CheckM2 on re-identified assemblies..."
checkm2 predict --threads "$THREADS" \
    --input "$FILTERED_DIR/" \
    --output-directory "$CHECKM2_OUT_DIR/" \
    --extension fasta

echo "All done. Four assemblies (C49, R16, R17, R2) are expected to still show"
echo "completeness < 90%; C25 needs a further species-specific re-filter (see protocol.qmd)."
