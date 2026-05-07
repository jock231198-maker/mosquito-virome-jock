#!/bin/bash
# STAR align
SECONDS=0
input_dir="$1"
find "$input_dir" -type f -name "*.fastq" | while read -r file; do
sample=$(basename "$file" .fastq) 
outdir="/Users/JK/Desktop/STARresults/STARalign/${sample}/"
mkdir -p "$outdir"
echo "mapping $sample"

read1="$input_dir/${sample}.fastq"
read2="$input_dir/${sample}.fastq"
if [[ -f "$read1" && -f "$read2" ]]; then
  STAR \
        --runThreadN 8 \
        --runMode alignReads \
        --genomeDir /Users/JK/Desktop/STARresults/Star_index_albo \
        --readFilesIn "$read1" "$read2" \
        --outFileNamePrefix "${outdir}${sample}_"
else
echo "Missing files for $sample R1 and $sample R2"
fi
done      

duration=$SECONDS
echo "Process finished in $(($duration / 60)) minutes and $(($duration %  60)) seconds."      

find "/Users/JK/Desktop/Prueba" -type f -name "*.fastq" | while read -r file; do