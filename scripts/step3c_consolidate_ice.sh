#!/bin/bash
# Consolidate ICEfinder2 per-isolate result JSONs into one place, and
# report which isolates had ICEs detected. Run after run_icefinder2_batch.sh.

ICE_RESULT_DIR="$HOME/shared-team/ronnie.dir/ICEfinder2_linux/result"
OUT_DIR="$HOME/shared-team/ronnie.dir/07_mge_sharing/03_ICE_detection"
mkdir -p "$OUT_DIR"

for isolate_dir in "$ICE_RESULT_DIR"/*/; do
    isolate=$(basename "$isolate_dir")
    summary="${isolate_dir}${isolate}_ICEsum.json"
    [ -f "$summary" ] && cp "$summary" "$OUT_DIR/${isolate}_ICEsum.json"
done

echo "Isolates with predicted ICEs:"
ls "$OUT_DIR/"*_ICEsum.json 2>/dev/null | while read -r f; do
    count=$(python -c "import json; d=json.load(open('$f')); print(len(d))" 2>/dev/null)
    [ "$count" -gt 0 ] && echo "  $(basename "$f" _ICEsum.json): $count ICE(s)"
done
