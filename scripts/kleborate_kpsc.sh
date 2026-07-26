#!/bin/bash
# KpSC (K. quasipneumoniae + K. variicola) virulence characterisation with
# Kleborate.

kleborate \
  -a 00_assemblies/KQ_*.fasta 00_assemblies/KV_*.fasta \
  -p kpsc \
  -o 08_virulence \
  --trim_headers
mv 08_virulence/klebsiella_pneumo_complex_output.txt 08_virulence/kleborate_kpsc.tsv
