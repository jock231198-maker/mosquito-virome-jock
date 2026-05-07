#!/bin/bash
# Define el archivo de salida
OUTPUT_TXT="read_counts.txt"
# Limpia o crea el archivo
> "$OUTPUT_TXT"
count_reads() {
  local f="$1"
  if [[ "$f" == *.gz ]]; then
    zcat "$f" | awk 'NR%4==1' | wc -l
  else
    awk 'NR%4==1' "$f" | wc -l
  fi
}
 
N_R1=$(count_reads "$OUT_R1")
N_R2=$(count_reads "$OUT_R2")
 {
  echo "=== Read Count Report ==="
  echo "  R1 reads written : $N_R1"
  echo "  R2 reads written : $N_R2"

  if [[ "$N_R1" -ne "$N_R2" ]]; then
    echo "  WARNING: R1 and R2 counts differ — check for name mismatches."
  fi

  if [[ "$N_R1" -ne "$N_READS" ]]; then
    echo "  NOTE: $((N_READS - N_R1)) read name(s) not found in FASTQ."
  fi

  echo ""
  echo "=== Done! Output files ==="
  echo "  $OUT_R1"
  echo "  $OUT_R2"
 }| tee "$OUTPUT_TXT"
 
echo ""
echo "=== Done! Output files ==="
echo "  $OUT_R1"
echo "  $OUT_R2"