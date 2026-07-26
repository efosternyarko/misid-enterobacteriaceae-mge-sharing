#!/bin/bash
# Symlink mob_recon outputs for the 23 misidentified isolates into the
# combined 01_mob_recon/ directory used by the MGE-sharing analysis,
# appending the species suffix that drives species parsing downstream.
# Run once only -- re-running fails on existing symlinks with a harmless
# "file exists" error (use ln -sf if you need to re-run).

PLASMIDS_DIR="$HOME/shared-team/ronnie.dir/06_plasmids"
RECON_DIR="$HOME/shared-team/ronnie.dir/07_mge_sharing/01_mob_recon"
mkdir -p "$RECON_DIR"

for isolate_dir in "$PLASMIDS_DIR"/*/; do
    isolate=$(basename "$isolate_dir")
    case "${isolate:0:2}" in
        EC) species="Ecoli" ;;
        KQ) species="Kquasipneumoniae" ;;
        KV) species="Kvariicola" ;;
        KA) species="Kaerogenes" ;;
        *) echo "Unknown prefix for $isolate -- skipping"; continue ;;
    esac
    ln -s "$(realpath "$isolate_dir")" "$RECON_DIR/${isolate}_${species}"
done

# Verify (expect 23 misidentified isolate directories)
ls "$RECON_DIR" | grep -v Kpneumoniae | wc -l
