#!/bin/bash
# ---------------------------------------------------------------------------
# extract_unmapped_worker.sh  — extraer pares NO mapeados a FASTQ (un BAM por elemento)
#
# Prep + envío:
#   find mapped -type f -name "*.bam" | sort > bam_list.txt
#   N=$(wc -l < bam_list.txt)
#   bsub -J "unmap[1-$N]%10" -o "logs/unmap.%J.%I.log" -e "logs/unmap.%J.%I.err" \
#        -q normal -n 4 -M 8000 \
#        -R "select[mem>8000] rusage[mem=8000] span[hosts=1]" \
#        "./extract_unmapped_worker.sh bam_list.txt unmapped_fastq"
# ---------------------------------------------------------------------------
set -euo pipefail

bam_list="${1:?Uso: extract_unmapped_worker.sh <bam_list.txt> <output_dir>}"
outdir="${2:?Falta output_dir}"

THREADS="${LSB_DJOB_NUMPROC:-4}"
CONDA_ENV="${CONDA_ENV:-bowtie2}"

eval "$(conda shell.bash hook)"
conda activate "$CONDA_ENV"

idx="${LSB_JOBINDEX:?Enviar como job array}"
bam=$(sed -n "${idx}p" "$bam_list")
sample=$(basename "$bam" .bam)
mkdir -p "$outdir"

OUT_R1="$outdir/${sample}_unmapped_R1.fastq.gz"
OUT_R2="$outdir/${sample}_unmapped_R2.fastq.gz"

echo "=== Unmapped $sample (host: $(hostname)) ==="
# -f 12 : read y mate NO mapeados (0x4 + 0x8) ; -F 256 : excluye alineamientos secundarios
samtools view -@ "$THREADS" -b -f 12 -F 256 "$bam" \
  | samtools sort -@ "$THREADS" -n - \
  | samtools fastq -@ "$THREADS" \
      -1 "$OUT_R1" -2 "$OUT_R2" \
      -0 /dev/null -s /dev/null -n

echo "Done: $sample"
echo "  $OUT_R1"
echo "  $OUT_R2"
