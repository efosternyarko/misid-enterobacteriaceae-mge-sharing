#!/bin/bash
# Aggregate the per-isolate quality-filtered AMRFinderPlus results
# (05_amr_filtered/*.tsv, produced by amr_qc_filter.sh) into one combined
# TSV, needed as input for step7_summary_table.R. Same pattern as
# step2b_aggregate_mobtyper.sh.

OUT_DIR="07_mge_sharing/07_amr_summary"
OUT="$OUT_DIR/all_amr_filtered_combined.tsv"
mkdir -p "$OUT_DIR"

header_file=$(find 05_amr_filtered -name "*.tsv" | head -1)
head -1 "$header_file" | awk '{print "isolate\t" $0}' > "$OUT"

for f in 05_amr_filtered/*.tsv; do
    isolate=$(basename "$f" _amr.tsv)
    tail -n +2 "$f" | awk -v iso="$isolate" 'NF > 1 {print iso "\t" $0}' >> "$OUT"
done

wc -l "$OUT"
