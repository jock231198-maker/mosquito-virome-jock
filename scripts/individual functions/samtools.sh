#!/bin/bash
# Ejemplo de comando para extraer solo las lecturas no mapeadas de un archivo BAM usando samtools
samtools view -b -f 4 mapped.bam > unmapped.bam

# Ejemplo de comando para convertir un archivo BAM a FASTQ usando samtools
samtools fastq unmapped.bam > unmapped.fastq


rnaviralspades.py --12 /Users/JK/Desktop/bowtie2/MapBow/trimm_cat_fastq/M1_S1.unmapped.fastq -o /Users/JK/Desktop/spades