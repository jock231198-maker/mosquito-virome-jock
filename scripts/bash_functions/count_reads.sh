#!/bin/bash
# count_unmapped_reads.sh

SECONDS=0
eval "$(conda shell.bash hook)"
conda activate bow

input_dir="$1"           # Directorio con los BAMs o FASTQs
result_from="$2"         # Nombre del experimento
THREADS=6                # Núcleos para pigz (descompresión paralela)
outdir="$1"
OUTPUT_TXT="${outdir}/read_counts.txt"
> "$OUTPUT_TXT"

count_reads() {
  local f="$1"
  if [[ "$f" == *.gz ]]; then
    pigz -dc -p "$THREADS" "$f" | awk 'NR%4==1' | wc -l
  else
    awk 'NR%4==1' "$f" | wc -l
  fi
}

#Cabecera del reporte
{
  echo "=== Read Count Report: ${result_from} ==="
  echo "Generated: $(date)"
  echo ""
} | tee -a "$OUTPUT_TXT"

# ── Loop sobre los BAMs para reconstruir nombres de muestra ─────────────────
# ── Loop sobre los FASTQs R1 para reconstruir nombres de muestra ────────────
find "$input_dir" -type f -name "*R1*.fastq.gz" | while read -r fastq_file; do

  # Extrae el nombre de muestra quitando el sufijo _unmapped_R1.fastq.gz
  sample=$(basename "$fastq_file" _unmapped_R1.fastq.gz)
  OUT_R1="${outdir}/${sample}_unmapped_R1.fastq.gz"
  OUT_R2="${outdir}/${sample}_unmapped_R2.fastq.gz"

  # Verificar que los FASTQ existan
  if [[ ! -f "$OUT_R1" || ! -f "$OUT_R2" ]]; then
    echo "  [SKIP] $sample — archivos FASTQ no encontrados" | tee -a "$OUTPUT_TXT"
    continue
  fi

  echo "=== Contando: $sample ===" | tee -a "$OUTPUT_TXT"

  # Contar R1 y R2 en paralelo con & + wait
  count_reads "$OUT_R1" > /tmp/_r1_count &
  PID_R1=$!
  count_reads "$OUT_R2" > /tmp/_r2_count &
  PID_R2=$!

  wait "$PID_R1" "$PID_R2"

  N_R1=$(cat /tmp/_r1_count | tr -d ' ')
  N_R2=$(cat /tmp/_r2_count | tr -d ' ')

  {
    echo "  R1 reads : $N_R1"
    echo "  R2 reads : $N_R2"

    if [[ "$N_R1" -ne "$N_R2" ]]; then
      echo "  ⚠ WARNING: R1 y R2 tienen conteos distintos — posible mismatch."
    else
      echo "  ✓ R1 y R2 están pareados correctamente."
    fi

    echo "  Archivos:"
    echo "    $OUT_R1"
    echo "    $OUT_R2"
    echo ""
  } | tee -a "$OUTPUT_TXT"

done

rm -f /tmp/_r1_count /tmp/_r2_count

duration=$SECONDS
echo "Proceso terminado en $(($duration / 3600))h $(($duration % 3600 / 60))m $(($duration % 60))s." | tee -a "$OUTPUT_TXT"
echo "Reporte guardado en: $OUTPUT_TXT"

conda deactivate