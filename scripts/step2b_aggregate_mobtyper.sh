#!/bin/bash
# Aggregate mob_recon's per-isolate mobtyper_results.txt (already generated
# with --run_typer during mob_recon, so primary/secondary cluster IDs are
# already assigned -- mob_cluster is not needed, it would just re-query the
# same reference database) across all 126 isolates into one combined TSV.

RECON_DIR="$HOME/shared-team/ronnie.dir/07_mge_sharing/01_mob_recon"
OUT="$HOME/shared-team/ronnie.dir/07_mge_sharing/02_mob_cluster/all_mobtyper_results.tsv"
mkdir -p "$(dirname "$OUT")"

# Write header from first file found
header_file=$(find "$RECON_DIR" -name "mobtyper_results.txt" | head -1)
head -1 "$header_file" | awk '{print "isolate\t" $0}' > "$OUT"

# Append all isolates with isolate ID prepended
for isolate_dir in "$RECON_DIR"/*/; do
    isolate=$(basename "$isolate_dir")
    f="${isolate_dir}mobtyper_results.txt"
    [ -f "$f" ] || continue
    tail -n +2 "$f" | awk -v iso="$isolate" 'NF > 1 {print iso "\t" $0}' >> "$OUT"
done

# Verify: expect one row per plasmid across all 126 isolates
wc -l "$OUT"
