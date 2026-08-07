#!/bin/bash
# ---------------------------------------------------------------------------
# checkv_worker.sh — CheckV sobre un ensamblaje viral por elemento
#
# Prep + envío:
#   source .../config.sh
#   find "$SCRATCH/spades" -type f -name "*_final.fasta" | sort > "$SCRATCH/fasta_list.txt"
#   N=$(wc -l < "$SCRATCH/fasta_list.txt")
#   bsub -J "checkv[1-$N]%6" \
#        -o "$LOGS_DIR/checkv.%J.%I.log" -e "$LOGS_DIR/checkv.%J.%I.err" \
#        -q normal -n 8 -M 16000 \
#        -R "select[mem>16000] rusage[mem=16000] span[hosts=1]" \
#        "./checkv_worker.sh $SCRATCH/fasta_list.txt"
#
# SALIDA (se conserva): $RESULTS_DIR/checkv/<sample>/quality_summary.tsv, etc.
# Los .tmp/ internos de CheckV se borran al terminar.
# ---------------------------------------------------------------------------
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

fasta_list="${1:?Uso: checkv_worker.sh <fasta_list.txt> [output_base] [checkv_db]}"
outbase="${2:-$RESULTS_DIR/checkv}"
db="${3:-$CHECKV_DB}"

THREADS="${LSB_DJOB_NUMPROC:-8}"
activate_env "$ENV_CHECKV"

idx="${LSB_JOBINDEX:?Enviar como job array}"
fasta=$(sed -n "${idx}p" "$fasta_list")
[[ -n "$fasta" ]] || { echo "ERROR: linea $idx vacia en $fasta_list"; exit 1; }
sample=$(basename "$fasta" .fasta)
outdir="$outbase/$sample"
mkdir -p "$outdir"

echo "CheckV $sample (host: $(hostname))"
checkv end_to_end "$fasta" "$outdir" -d "$db" -t "$THREADS"

# tmp/ de CheckV: intermedios de prodigal/DIAMOND, no se necesitan
rm -rf "$outdir/tmp"

echo "Done: $outdir"
