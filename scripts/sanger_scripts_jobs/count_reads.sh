#!/bin/bash
# ---------------------------------------------------------------------------
# count_reads.sh  — contar lecturas en los FASTQ no-mapeados y generar reporte
#
# Escanea una carpeta (no necesita array). Para muchas muestras puedes enviarlo:
#   bsub -o "logs/count.%J.log" -e "logs/count.%J.err" \
#        -q normal -n 6 -M 4000 \
#        -R "select[mem>4000] rusage[mem=4000] span[hosts=1]" \
#        "./count_reads.sh unmapped_fastq M57"
# ---------------------------------------------------------------------------
set -euo pipefail
SECONDS=0

input_dir="${1:?Uso: count_reads.sh <input_dir> <result_from>}"
result_from="${2:?Falta result_from}"
THREADS="${LSB_DJOB_NUMPROC:-6}"
CONDA_ENV="${CONDA_ENV:-bow}"

eval "$(conda shell.bash hook)"
conda activate "$CONDA_ENV"

outdir="$input_dir"
OUTPUT_TXT="${outdir}/read_counts.txt"
: > "$OUTPUT_TXT"

count_reads() {
  local f="$1"
  if [[ "$f" == *.gz ]]; then
    pigz -dc -p "$THREADS" "$f" | awk 'NR%4==1' | wc -l
  else
    awk 'NR%4==1' "$f" | wc -l
  fi
}

{
  echo "=== Read Count Report: ${result_from} ==="
  echo "Generated: $(date)"
  echo ""
} | tee -a "$OUTPUT_TXT"

find "$input_dir" -type f -name "*R1*.fastq.gz" | while read -r fastq_file; do
  sample=$(basename "$fastq_file" _unmapped_R1.fastq.gz)
  OUT_R1="${outdir}/${sample}_unmapped_R1.fastq.gz"
  OUT_R2="${outdir}/${sample}_unmapped_R2.fastq.gz"

  if [[ ! -f "$OUT_R1" || ! -f "$OUT_R2" ]]; then
    echo "  [SKIP] $sample — archivos FASTQ no encontrados" | tee -a "$OUTPUT_TXT"
    continue
  fi

  echo "=== Contando: $sample ===" | tee -a "$OUTPUT_TXT"

  # Temporales ÚNICOS con mktemp (antes /tmp/_r1_count fijo se pisaba entre
  # trabajos paralelos en el mismo nodo)
  tmp_r1=$(mktemp); tmp_r2=$(mktemp)
  count_reads "$OUT_R1" > "$tmp_r1" &
  PID_R1=$!
  count_reads "$OUT_R2" > "$tmp_r2" &
  PID_R2=$!
  wait "$PID_R1" "$PID_R2"

  N_R1=$(tr -d ' ' < "$tmp_r1")
  N_R2=$(tr -d ' ' < "$tmp_r2")
  rm -f "$tmp_r1" "$tmp_r2"

  {
    echo "  R1 reads : $N_R1"
    echo "  R2 reads : $N_R2"
    if [[ "$N_R1" -ne "$N_R2" ]]; then
      echo "  WARNING: R1 y R2 tienen conteos distintos — posible mismatch."
    else
      echo "  OK: R1 y R2 pareados correctamente."
    fi
    echo "  Archivos:"
    echo "    $OUT_R1"
    echo "    $OUT_R2"
    echo ""
  } | tee -a "$OUTPUT_TXT"
done

dur=$SECONDS
echo "Proceso terminado en $((dur/3600))h $((dur%3600/60))m $((dur%60))s." | tee -a "$OUTPUT_TXT"
echo "Reporte guardado en: $OUTPUT_TXT"
