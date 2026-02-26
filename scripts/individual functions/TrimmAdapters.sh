#!/bin/bash
conda activate trimmomatic_env
INPUT_DIR="/Users/JK/Desktop/JKViromedata/Aa" #fastq base
OUTPUT_DIR="/Users/JK/Desktop/SecRef4TrimmAdapters " # carpeta salida

for R1 in "$INPUT_DIR"/*_R1.001.fastq
do
    BASENAME=$(basename "$R1" R1.001.fastq)
    R2="$INPUT_DIR/${BASENAME}_R2.001.fastq"
    
   trimmomatic PE -threads 10 -phred33 \
    "$R1" "$R2" \
    "$OUTPUT_DIR/${BASENAME}_R1_paired.fastq.gz" \
    "$OUTPUT_DIR/${BASENAME}_R1_unpaired.fastq.gz" \
    "$OUTPUT_DIR/${BASENAME}_R2_paired.fastq.gz" \
    "$OUTPUT_DIR/${BASENAME}_R2_unpaired.fastq.gz" \
  ILLUMINACLIP:/Users/Parsimony/miniconda3/envs/trimmomatic/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa:2:30:10:2:keepBothReads \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:25 MINLEN:35
done
conda deactivate
