#!/bin/bash
mkdir -p /Users/JK/Desktop/genomad/genomad_output/
genomad end-to-end \
    /Users/JK/Desktop/spades/COBRA_M1_S1/COBRA_M1_S1_final.fasta \
    /Users/JK/Desktop/genomad/genomad_output/ \
    /Users/JK/Desktop/genomad/genomad_db/genomad_db \
    --cleanup \
    --splits 8

cat /Users/JK/Desktop/genomad/genomad_output/COBRA_M1_S1_final_marker_classification/COBRA_M1_S1_final_marker_classification.tsv # check_important output file
cat /Users/JK/Desktop/genomad/genomad_output/COBRA_M1_S1_final_annotate/COBRA_M1_S1_final_taxonomy.tsv # check_annotate output file