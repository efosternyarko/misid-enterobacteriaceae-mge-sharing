#!/bin/bash
# Pairwise all-vs-all BLASTn of extracted ICE sequences to identify shared
# ICEs across species at >=95% nucleotide identity and >=80% query coverage.

makeblastdb \
    -in 03_MGE_sharing/03_ICE_detection/all_ICE_sequences.fasta \
    -dbtype nucl \
    -out 03_MGE_sharing/04_ICE_blast/ICE_db

blastn \
    -query 03_MGE_sharing/03_ICE_detection/all_ICE_sequences.fasta \
    -db 03_MGE_sharing/04_ICE_blast/ICE_db \
    -perc_identity 95 \
    -qcov_hsp_perc 80 \
    -outfmt "6 qseqid sseqid pident qcovs length mismatch gapopen qstart qend sstart send evalue bitscore" \
    -out 03_MGE_sharing/04_ICE_blast/ICE_pairwise_blast.tsv \
    -num_threads 8

# Remove self-hits
awk '$1 != $2' 03_MGE_sharing/04_ICE_blast/ICE_pairwise_blast.tsv \
    > 03_MGE_sharing/04_ICE_blast/ICE_blast_no_selfhits.tsv
