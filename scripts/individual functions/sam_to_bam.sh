#!/bin/bash
#samtools sam to bam
#runlike this ./sam_to_bam.sh "/directory/to/my/viromedata/"
SECONDS=0

input_dir="$1"
find "$input_dir" -type f -name "*.sam" | while read -r file; do
  sample=$(basename "$file" | cut -d'_' -f1-2)
  outdir="/Users/JK/Desktop/Bow/MapBow/${sample}/"
  echo "bamming $sample"

read1="$input_dir/${sample}.sam"
read2="$input_dir/${sample}.sam"
if [[ -f "$read1" && -f "$read2" ]]; then
samtools view -bS "${outdir}/${sample}.sam" > "${outdir}/${sample}.bam"
else
echo "Missing files for $sample R1 and $sample R2"
fi
done
duration=$SECONDS
echo "Process finished in $(($duration / 60)) minutes and $(($duration % 60)) seconds."