#!/bin/bash
#samtools sam to bam
#runlike this ./sam_to_bam.sh "/directory/to/my/viromedata/"
SECONDS=0
input_dir="$1"
find "$input_dir" -type f -name "*.sam" | while read -r file; do
  sample=$(basename "$file" .sam)
  outdir="/Users/JK/Desktop/Bow/MapBow/${sample}/"
  echo "bamming $sample"

read1="$input_dir/${sample}.sam"
if [[ -f "$read1" ]]; then
samtools view -bS "${outdir}/${sample}.sam" > "${outdir}/${sample}.bam"
else
echo "Missing files for $sample"

fi
done
duration=$SECONDS
echo "Process finished in $(($duration / 60)) minutes and $(($duration % 60)) seconds."
