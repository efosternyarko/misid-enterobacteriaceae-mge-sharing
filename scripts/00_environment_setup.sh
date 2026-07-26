#!/bin/bash
# Conda environment setup and project directory structure.
# Run once when setting up a fresh working environment.

# --- misid_ecoli: main environment for QC, typing, annotation, AMR/plasmid tools ---
conda create -n misid_ecoli python=3.11
conda activate misid_ecoli
conda install -c bioconda -c conda-forge \
  mash kraken2 mlst prokka bakta amrfinderplus \
  kaptive kleborate virulencefinder \
  chewbbaca blast fastp
conda install -c bioconda -c conda-forge mob_suite "pandas<2.0"
conda install -c bioconda -c conda-forge ncbi-datasets-cli jq

# --- checkm2_py312: CheckM2 needs a newer Python than the main environment ---
conda create -n checkm2_py312 -c bioconda -c conda-forge python=3.12 checkm2 kraken2 krakentools
conda activate checkm2_py312
checkm2 database --setdblocation /home/jovyan/shared-team/ronnie.dir/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd
checkm2 testrun

# --- ectyper_env: ECTyper has its own dependency constraints ---
conda create -n ectyper_env -c bioconda -c conda-forge python=3.10 ectyper

# --- icefinder2_env: ICEfinder2 and its dependencies ---
conda create -n icefinder2_env -c bioconda -c conda-forge python=3.8 \
    hmmer blast kraken2 seqkit prodigal prokka \
    macsyfinder defense-finder biopython ete3
conda activate icefinder2_env
# download and extract ICEfinder2.1 from https://github.com/EBI-Metagenomics/icefinder2/
tar -xvzf ICEfinder2.1_linux.tar.gz
# find the paths for config.ini
which hmmsearch
which blastn
which blastp
which kraken2
which seqkit
which prodigal
which prokka
which defense-finder
which macsyfinder
# edit config.ini by completing the absolute paths to tools identified above
nano config.ini
# BioPython >= 1.80 removed Bio.SeqUtils.GC, which ICEfinder2 uses internally --
# downgrade before running ICEfinder2 (see merge_misidentified_chromosomes.sh)
conda install -c conda-forge biopython=1.79 -y

# --- Project directory structure ---
mkdir -p 00_assemblies 01_qc 02_species 03_annotation 04_mlst 05_amr \
    06_plasmids 07_mge_sharing 08_virulence 09_phylogeny 10_figures \
    11_excluded_assemblies
