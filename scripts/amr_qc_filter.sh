#!/bin/bash
# Quality-filter AMRFinderPlus results: drop PARTIALX (frameshifted, likely
# non-functional) hits, and require >=70% coverage and >=90% identity.
#
# Column positions are looked up by header name, not hardcoded — AMRFinderPlus's
# column count shifts depending on invocation (e.g. whether --protein/--gff were
# given alongside --nucleotide, which adds/omits the Contig id/Start/Stop/Strand
# columns), so a fixed column number for one run can silently be wrong for another.

mkdir -p 05_amr_filtered

for file in 05_amr/*.tsv; do
    [ -e "$file" ] || continue
    filename=$(basename "$file")
    awk 'BEGIN {FS="\t"; OFS="\t"}
        NR==1 {
            for (i=1; i<=NF; i++) {
                if ($i == "Method")                       method_col = i
                if ($i == "% Coverage of reference")      cov_col    = i
                if ($i == "% Identity to reference")      id_col     = i
            }
            if (!method_col || !cov_col || !id_col) {
                print "ERROR: expected column(s) not found in header of " FILENAME > "/dev/stderr"
                exit 1
            }
            print; next
        }
        $method_col != "PARTIALX" && $cov_col >= 70 && $id_col >= 90 {print}
    ' "$file" > "05_amr_filtered/$filename"
    echo "Filtered: $filename"
done

echo "All clinical isolates processed. Quality-filtered results are saved in 05_amr_filtered/"
