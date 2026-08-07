#!/bin/bash
# ---------------------------------------------------------------------------
# checkv_worker.sh  — CheckV sobre un ensamblaje viral por elemento
#
# Prep + envío:
#   find spades -type f -name "*_final.fasta" | sort > fasta_list.txt
#   N=$(wc -l < fasta_list.txt)
#   bsub -J "checkv[1-$N]%6" -o "logs/checkv.%J.%I.log" -e "logs/checkv.%J.%I.err" \
#        -q normal -n 8 -M 16000 \
#        -R "select[mem>16000] rusage[mem=16000] span[hosts=1]" \
#        "./checkv_worker.sh fasta_list.txt checkv_out refs/checkv-db-v1.5"
# ---------------------------------------------------------------------------
set -euo pipefail

fasta_list="${1:?Uso: checkv_worker.sh <fasta_list.txt> <output_base> <checkv_db>}"
outbase="${2:?Falta output_base}"
CHECKV_DB="${3:?Falta ruta a la base de datos CheckV}"

THREADS="${LSB_DJOB_NUMPROC:-8}"
CONDA_ENV="${CONDA_ENV:-checkv}"

eval "$(conda shell.bash hook)"
conda activate "$CONDA_ENV"

idx="${LSB_JOBINDEX:?Enviar como job array}"
fasta=$(sed -n "${idx}p" "$fasta_list")
sample=$(basename "$fasta" .fasta)
outdir="$outbase/$sample"
mkdir -p "$outdir"

echo "CheckV $sample (host: $(hostname))"
checkv end_to_end "$fasta" "$outdir" -d "$CHECKV_DB" -t "$THREADS"
echo "Done: $outdir"
