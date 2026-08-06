#!/bin/bash
# ---------------------------------------------------------------------------
# genomad_worker.sh — identificación/clasificación viral con geNomad
#                     (un ensamblaje por elemento)
#
# Prep + envío:
#   source .../config.sh
#   find "$SCRATCH/spades" -type f -name "*_final.fasta" | sort > "$SCRATCH/fasta_list.txt"
#   N=$(wc -l < "$SCRATCH/fasta_list.txt")
#   bsub -J "genomad[1-$N]%4" \
#        -o "$LOGS_DIR/genomad.%J.%I.log" -e "$LOGS_DIR/genomad.%J.%I.err" \
#        -q normal -n 8 -M 32000 \
#        -R "select[mem>32000] rusage[mem=32000] span[hosts=1]" \
#        "./genomad_worker.sh $SCRATCH/fasta_list.txt"
#
# --splits intercambia memoria por tiempo: súbelo si te quedas sin RAM.
# --cleanup ya borra los intermedios de geNomad automáticamente.
# SALIDA (se conserva): $RESULTS_DIR/genomad/<sample>/
# ---------------------------------------------------------------------------
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

fasta_list="${1:?Uso: genomad_worker.sh <fasta_list.txt> [output_base] [genomad_db]}"
outbase="${2:-$RESULTS_DIR/genomad}"
db="${3:-$GENOMAD_DB}"

SPLITS="${SPLITS:-8}"
activate_env "$ENV_GENOMAD"

idx="${LSB_JOBINDEX:?Enviar como job array}"
fasta=$(sed -n "${idx}p" "$fasta_list")
[[ -n "$fasta" ]] || { echo "ERROR: linea $idx vacia en $fasta_list"; exit 1; }
sample=$(basename "$fasta" .fasta)
outdir="$outbase/$sample"
mkdir -p "$outdir"

echo "geNomad $sample (host: $(hostname))"
genomad end-to-end "$fasta" "$outdir" "$db" --cleanup --splits "$SPLITS"
echo "Done: $outdir"
