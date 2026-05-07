#!/bin/bash
find "$input_dir" -type f -name "*fastq.gz" | while read -r file; do

sample=$(basename "$file" | cut -d'_' -f1-2)
endname=$(basename "$file" | cut -d'_' -f3-$(echo "$file" | tr '_' '\n' | wc -l))
outdir="/Users/JK/Desktop/bowtie2/MapBow/${result_from}/"
mkdir -p "$outdir"
echo "mapping $sample"

read1="$input_dir/${sample}_${endname}"
read2="$input_dir/${sample}_R2_001_paired.fastq.gz"



find "/Users/JK/Desktop/jkviromedata/cat_fastq" -type f -name "*fastq.gz" | while read -r file; do

sample=$(basename "$file" | cut -d'_' -f1-2)
endname=$(basename "$file" | cut -d'_' -f4-)
read=$(basename "$file" | cut -d'_' -f3)

echo "mapping $sample"
echo "read: $read"
echo "endname: $endname"

if [[ $read == "R1" ]]; then
read1="$input_dir/${sample}_R1_${endname}"
fi

if [[ $read == "R2" ]]; then
read2="$input_dir/${sample}_R2_${endname}"
fi

echo "read1: $read1"
echo "read2: $read2"
done

if [[ -f "$read1" && -f "$read2" ]]; then
echo "exist $read1 and $read2"
else
echo "Missing files for $sample R1 and $sample R2"
fi
done