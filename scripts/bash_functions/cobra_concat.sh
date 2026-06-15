#!/bin/bash
SECONDS=0
eval "$(conda shell.bash hook)"
conda activate cobra # Inicializa conda # Activar el entorno con cobra-meta instalado

cobra-meta \
    -q /Users/JK/Desktop/spades/viral_query_contigs.fasta \ #
    -f /Users/JK/Desktop/spades/Assembly/trimm_cat_fastq_unmapped/M1_S1/contigs.fasta \
    -a metaspades \
    -mink 21 \
    -maxk 55 \
    -m /Users/JK/Desktop/bowtie2/spades_maps/m1_sorted_mapped.bam \
    -c /Users/JK/Desktop/bowtie2/spades_maps/coverage_for_cobra.tsv \
    -t 8 \
    -o /Users/JK/Desktop/spades/COBRA_M1_S1