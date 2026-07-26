#!/bin/bash
# Merge multi-contig chromosome FASTAs (misidentified isolates only -- the
# 103 K. pneumoniae isolates are not run through ICEfinder2) into a single
# contig per isolate, separated by 100 N spacers, since ICEfinder2 requires
# one merged contig per input FASTA.
#
# Do NOT use `seqkit concat` for this -- it concatenates across files by
# matching sequence IDs, it does not merge contigs within a single file.
# Empty output resulted from that approach; use the biopython script below
# instead. Requires the icefinder2_env conda environment (biopython 1.79).

RECON_DIR="$HOME/shared-team/ronnie.dir/07_mge_sharing/01_mob_recon"
MERGED_DIR="$HOME/shared-team/ronnie.dir/07_mge_sharing/03_ICE_detection/merged_chromosomes"
mkdir -p "$MERGED_DIR"

for isolate_dir in "$RECON_DIR"/*/; do
    isolate=$(basename "$isolate_dir")
    [[ "$isolate" == *"Kpneumoniae"* ]] && continue   # misidentified isolates only
    chr="${isolate_dir}chromosome.fasta"
    [ -f "$chr" ] || continue

    python3 -c "
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
seqs = list(SeqIO.parse(sys.argv[1], 'fasta'))
merged = SeqRecord(Seq(('N'*100).join(str(s.seq) for s in seqs)),
                   id=sys.argv[2], description='')
SeqIO.write(merged, sys.argv[3], 'fasta')
print(f'{sys.argv[2]}: {len(seqs)} contig(s), {len(merged.seq):,} bp')
" "$chr" "${isolate}_merged" "$MERGED_DIR/${isolate}_chromosome_merged.fasta"
done

# Verify -- all files must be non-zero
ls -lh "$MERGED_DIR"
