#!/bin/bash
# extract_unmapped.sh
SECONDS=0
eval "$(conda shell.bash hook)"
conda activate bow

# run like this: ./extract_unmapped.sh "/directory/to/my/bam/" "experiment_name"

input_dir="$1"
result_from="$2"
outdir="/Users/JK/Desktop/bowtie2/UnmappedFastq/${result_from}/"
mkdir -p "$outdir"

find "$input_dir" -type f -name "*.bam" | while read -r bam; do

  sample=$(basename "$bam" .bam)
  echo ""
  echo "=== Processing $sample ==="

  OUT_R1="${outdir}/${sample}_unmapped_R1.fastq.gz"
  OUT_R2="${outdir}/${sample}_unmapped_R2.fastq.gz"

  samtools view -@ 4 -b -f 12 -F 256 "$bam" \
    | samtools sort -@ 4 -n - \
    | samtools fastq \
        -@ 4 \
        -1 "$OUT_R1" \
        -2 "$OUT_R2" \
        -0 /dev/null \
        -s /dev/null \
        -n

  echo "  Done: $sample"
  echo "  $OUT_R1"
  echo "  $OUT_R2"

done

duration=$SECONDS
echo ""
echo "Process finished in $(($duration / 3600)) hour(s), $(($duration % 3600 / 60)) minute(s) and $(($duration % 60)) second(s)."
conda deactivate
```