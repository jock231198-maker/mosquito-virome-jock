#!/bin/bash
# ---------------------------------------------------------------------------
# extract_unmapped_worker.sh — extraer pares NO mapeados a FASTQ (un BAM por elemento)
#
# Prep + envío:
#   source .../config.sh
#   find "$SCRATCH/mapped" -type f -name "*.bam" | sort > "$SCRATCH/bam_list.txt"
#   N=$(wc -l < "$SCRATCH/bam_list.txt")
#   bsub -J "unmap[1-$N]%10" \
#        -o "$LOGS_DIR/unmap.%J.%I.log" -e "$LOGS_DIR/unmap.%J.%I.err" \
#        -q normal -n 4 -M 8000 \
#        -R "select[mem>8000] rusage[mem=8000] span[hosts=1]" \
#        "./extract_unmapped_worker.sh $SCRATCH/bam_list.txt"
#
# SALIDA: $SCRATCH/unmapped_fastq/  (entrada del ensamblaje; borrable tras ensamblar)
# Tras este paso, los BAM de $SCRATCH/mapped/ ya no se necesitan.
# ---------------------------------------------------------------------------
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

bam_list="${1:?Uso: extract_unmapped_worker.sh <bam_list.txt> [output_dir]}"
outdir="${2:-$SCRATCH/unmapped_fastq}"

THREADS="${LSB_DJOB_NUMPROC:-4}"
activate_env "$ENV_SAMTOOLS"

idx="${LSB_JOBINDEX:?Enviar como job array}"
bam=$(sed -n "${idx}p" "$bam_list")
[[ -n "$bam" ]] || { echo "ERROR: linea $idx vacia en $bam_list"; exit 1; }
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

# Conteo de lecturas no mapeadas -> resultado permanente (tabla de rendimiento viral)
countdir="$RESULTS_DIR/unmapped_counts"
mkdir -p "$countdir"
n=$(( $(zcat "$OUT_R1" | wc -l) / 4 ))
printf '%s\t%s\n' "$sample" "$n" > "$countdir/${sample}_unmapped_pairs.tsv"

echo "Done: $sample  ($n pares no mapeados)"
echo "  $OUT_R1"
echo "  $OUT_R2"
