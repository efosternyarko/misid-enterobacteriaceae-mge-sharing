#!/bin/bash
# Quality-filter AMRFinderPlus results: drop PARTIALX (frameshifted, likely
# non-functional) hits, and require >=70% coverage and >=90% identity.
# Column positions ($13 Method, $16 % Coverage, $17 % Identity) match this
# AMRFinderPlus version's TSV layout -- re-check column numbers with
# `head -1 05_amr/<isolate>_amr.tsv | tr '\t' '\n' | cat -n` if you upgrade
# AMRFinderPlus and this stops matching.

mkdir -p 05_amr_filtered

for file in 05_amr/*.tsv; do
    [ -e "$file" ] || continue
    filename=$(basename "$file")
    awk 'BEGIN {FS="\t"; OFS="\t"}
        NR==1 {print; next}
        # Keep row IF: method is not PARTIALX AND coverage >= 70 AND identity >= 90
        $13 != "PARTIALX" && $16 >= 70 && $17 >= 90 {print}
    ' "$file" > "05_amr_filtered/$filename"
    echo "Filtered: $filename"
done

echo "All clinical isolates processed. Quality-filtered results are saved in 05_amr_filtered/"
