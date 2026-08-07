#!/bin/bash
# ---------------------------------------------------------------------------
# genomad_worker.sh  — identificación/clasificación viral con geNomad
#                      (un ensamblaje por elemento)
#
# Prep + envío:
#   find spades -type f -name "*_final.fasta" | sort > fasta_list.txt
#   N=$(wc -l < fasta_list.txt)
#   bsub -J "genomad[1-$N]%4" -o "logs/genomad.%J.%I.log" -e "logs/genomad.%J.%I.err" \
#        -q normal -n 8 -M 32000 \
#        -R "select[mem>32000] rusage[mem=32000] span[hosts=1]" \
#        "./genomad_worker.sh fasta_list.txt genomad_out refs/genomad_db"
#
# --splits intercambia memoria por tiempo: súbelo si te quedas sin RAM.
# ---------------------------------------------------------------------------
set -euo pipefail

fasta_list="${1:?Uso: genomad_worker.sh <fasta_list.txt> <output_base> <genomad_db>}"
outbase="${2:?Falta output_base}"
GENOMAD_DB="${3:?Falta ruta a la base de datos geNomad}"

SPLITS="${SPLITS:-8}"
CONDA_ENV="${CONDA_ENV:-genomad}"

eval "$(conda shell.bash hook)"
conda activate "$CONDA_ENV"

idx="${LSB_JOBINDEX:?Enviar como job array}"
fasta=$(sed -n "${idx}p" "$fasta_list")
sample=$(basename "$fasta" .fasta)
outdir="$outbase/$sample"
mkdir -p "$outdir"

echo "geNomad $sample (host: $(hostname))"
genomad end-to-end "$fasta" "$outdir" "$GENOMAD_DB" --cleanup --splits "$SPLITS"
echo "Done: $outdir"
