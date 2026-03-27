#!/bin/bash
#bowtiemap
input_dir="$1"
find "$input_dir" -type f -name "*_R1_001.fastq.gz" | while read -r file; do
  sample=$(basename "$file" | cut -d'_' -f1-2)
  outdir="/Users/JK/Desktop/Bow/MapBow/${sample}/"
  mkdir -p "$outdir"
  echo "mapping $sample"

read1="${sample}_R1_001.fastq.gz"
read2="${sample}_R2_001.fastq.gz"
if [[ -f "$read1" && -f "$read2" ]]; then
bowtie2 -x /Users/JK/Desktop/Bow/IndexBow/Aae_index/Aae_index -1 "$input_dir"/"$read1" -2 "$input_dir"/"$read2" -S /${outdir}/${sample}.sam
else
echo "Missing files for $sample R1 and $sample R2"
fi
done