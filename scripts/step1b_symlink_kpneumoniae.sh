#!/bin/bash
# Symlink the pre-processed Mills et al. 2024 K. pneumoniae mob_recon
# outputs into the same combined 01_mob_recon/ directory as the
# misidentified isolates (step1a_symlink_misidentified.sh). Run once only.

KPN_DIR="$HOME/shared-team/ronnie.dir/07_mge_sharing/kpneumoniae_mills2024"
RECON_DIR="$HOME/shared-team/ronnie.dir/07_mge_sharing/01_mob_recon"

for strain_dir in "$KPN_DIR"/*/; do
    strain=$(basename "$strain_dir")
    isolate="${strain}_Kpneumoniae"
    outdir="$RECON_DIR/$isolate"
    mkdir -p "$outdir"

    # Symlink chromosome (rename to mob_recon convention)
    ln -s "${strain_dir}${strain}_chromosome.fasta" "$outdir/chromosome.fasta"

    # Symlink plasmids (strip strain prefix to match mob_recon naming)
    for pls in "${strain_dir}${strain}_plasmid_"*.fasta; do
        [ -f "$pls" ] || continue
        pname=$(basename "$pls" | sed "s/^${strain}_//")
        ln -s "$pls" "$outdir/$pname"
    done

    # Symlink support files
    for f in contig_report.txt mobtyper_results.txt; do
        [ -f "${strain_dir}$f" ] && ln -s "${strain_dir}$f" "$outdir/$f"
    done
done

# Verify: expect 103 subdirectories
ls "$RECON_DIR" | wc -l
