#!/bin/bash
# E. coli O:H serotyping (ECTyper) and K-locus typing (Kaptive, two-pass:
# G2/G3 groups first, then G1/G4 on isolates left untypeable by the first pass).

# --- O:H serotyping ---
ectyper -i 00_assemblies/EC_*.fasta -o 08_virulence/ectyper/ -opid 70

# --- K-locus typing: Step 1, G2/G3 ---
cd 08_virulence/ || exit 1
unzip EC-K-typing-G1G4-main.zip
    # download from https://github.com/efosternyarko/EC-K-typing-G1G4
mkdir -p kaptive_G23
kaptive assembly \
    ./EC-K-typing-G1G4-main/DB/EC-K-typing_group2and3_v3.0.0.gbk \
    /shared/team/ronnie.dir/00_assemblies/EC_*.fasta \
    -o kaptive_G23/kaptive_results.tsv

# --- K-locus typing: Step 2, G1/G4 on isolates "Untypeable" in Step 1 ---
grep -w "Untypeable" kaptive_G23/kaptive_results.tsv | cut -f1 > untypeable_ids.txt
mkdir -p untypeable/
while IFS= read -r id; do
    cp "/shared/team/ronnie.dir/00_assemblies/${id}.fasta" untypeable/
done < untypeable_ids.txt

mkdir -p kaptive_G14
kaptive assembly \
    EC-K-typing-G1G4-main/DB/EC-K-typing_group1and4_v1.3.gbk \
    untypeable/*.fasta \
    --scores kaptive_G14/kaptive_scores.tsv \
    -t 8

python3 EC-K-typing-G1G4-main/scripts/normalise_kaptive_scores.py \
    --db EC-K-typing-G1G4-main/DB/EC-K-typing_group1and4_v1.3.gbk \
    --in kaptive_G14/kaptive_scores.tsv \
    --out kaptive_G14/kaptive_results_norm.tsv
