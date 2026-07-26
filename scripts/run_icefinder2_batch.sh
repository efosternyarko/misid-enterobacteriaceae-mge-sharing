#!/bin/bash
# Run ICEfinder2 on all merged misidentified-isolate chromosomes.
# ICEfinder2 takes 1-2 minutes per genome -- run this inside `screen` so
# the job survives a disconnected session:
#   screen -S icefinder2
#   (run the commands below)
#   Ctrl+A then D to detach; `screen -r icefinder2` to reattach.
#
# Only -i (input FASTA) and -t Single are valid flags for ICEfinder2 -- do
# NOT add -o or --config, they do not exist. Output goes to result/ inside
# ICEfinder2_linux/; config is read from config.ini in the same directory.

conda activate icefinder2_env
cd "$HOME/shared-team/ronnie.dir/ICEfinder2_linux" || exit 1

MERGED_DIR="$HOME/shared-team/ronnie.dir/07_mge_sharing/03_ICE_detection/merged_chromosomes"

for fasta in "$MERGED_DIR"/*.fasta; do
    sample=$(basename "$fasta" .fasta)
    if [ -d "result/${sample}" ]; then
        echo "SKIP: $sample (already done)"
        continue
    fi
    echo "Running: $sample"
    python ICEfinder2.py -i "$fasta" -t Single
    echo "Done: $sample"
done

echo "All samples complete."
