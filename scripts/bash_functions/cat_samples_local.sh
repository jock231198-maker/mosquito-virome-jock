#!/bin/bash
#Concatenate samples in the samples path
for sample in $(ls *_L001_R1_001.fastq.gz | cut -d'_' -f1-2 | sort -u); do
outdir="/Users/JK/Desktop/JKViromedata/cat_fastq/${sample}/"
mkdir -p "$outdir"
echo "cat $sample"
    r1_l001="${sample}_L001_R1_001.fastq.gz"
    r1_l002="${sample}_L002_R1_001.fastq.gz"
    r2_l001="${sample}_L001_R2_001.fastq.gz"
    r2_l002="${sample}_L002_R2_001.fastq.gz"
    r1_out="${sample}_R1_001.fastq.gz"
    r2_out="${sample}_R2_001.fastq.gz"
    
    # Concatenar R1 si ambos archivos existen
    if [[ -f "$r1_l001" && -f "$r1_l002" ]]; then
        cat "$r1_l001" "$r1_l002" > ${outdir}/"$r1_out"
        echo "Created: $r1_out"
    else
        echo "WARNING: Missing files for $sample R1"
    fi

    # Concatenar R2 si ambos archivos existen
    if [[ -f "$r2_l001" && -f "$r2_l002" ]]; then
        cat "$r2_l001" "$r2_l002" > ${outdir}/"$r2_out"
        echo "Created: $r2_out"
    else
        echo "WARNING: Missing files for $sample R2"
    fi
done