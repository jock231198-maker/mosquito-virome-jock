#!/bin/bash

diamond blastx \
    -q /Users/JK/Desktop/spades/contigs_1000bp.fasta \
    -d /Users/JK/Desktop/diamond/diamond_db/viral_refseq_proteins.dmnd \
    -o /Users/JK/Desktop/diamond/diamond_results.tsv \
    -f 6 qseqid sseqid pident length evalue bitscore stitle \
    -e 1e-5 \
    --more-sensitive \
    -p 8

# Correr Diamond con formato extendido
diamond blastx \
    -q /Users/JK/Desktop/spades/contigs_1000bp.fasta \
    -d /Users/JK/Desktop/diamond/diamond_db/viral_refseq_proteins.dmnd \
    -o /Users/JK/Desktop/diamond/diamond_results_extended.tsv \
    -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle qcovhsp scovhsp \
    -e 1e-5 \
    --more-sensitive \
    -p 8

# Agregar header al archivo
echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstitle\tqcovhsp\tscovhsp" | \
cat - /Users/JK/Desktop/diamond/diamond_results_extended.tsv > /tmp/temp_diamond.tsv && \
mv /tmp/temp_diamond.tsv /Users/JK/Desktop/diamond/diamond_results_extended.tsv