#!/bin/bash
set +u
source ~/miniconda3/bin/activate CIRIquant_env

prep_CIRIquant -i data/ciri_de_samples.txt \
    --lib results/3_ciri/4_ciri_de/library_info.csv \
    --circ results/3_ciri/4_ciri_de/circRNA_info.csv \
    --bsj results/3_ciri/4_ciri_de/circRNA_bsj.csv \
    --ratio results/3_ciri/4_ciri_de/circRNA_ratio.csv > results/3_ciri/4_ciri_de/prep_CIRIquant.log 2>&1
 
prepDE.py -i data/ciri_de_stringtie_gene.txt \
    -g results/3_ciri/4_ciri_de/gene_count_matrix.csv \
    -t results/3_ciri/4_ciri_de/transcript_count_matrix.csv > results/3_ciri/4_ciri_de/prepDE.log 2>&1

CIRI_DE_replicate \
    --lib  results/3_ciri/4_ciri_de/library_info.csv \
    --bsj  results/3_ciri/4_ciri_de/circRNA_bsj.csv \
    --gene results/3_ciri/4_ciri_de/gene_count_matrix.csv \
    --out  results/3_ciri/4_ciri_de/circRNA_de.csv \
    --out2 results/3_ciri/4_ciri_de/gene_de.csv > results/3_ciri/4_ciri_de/CIRI_DE_replicate.log 2>&1
