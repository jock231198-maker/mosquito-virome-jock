#!/bin/bash
# ---------------------------------------------------------------------------
# fastqc_worker.sh — FastQC sobre UN fastq.gz por elemento de un job array LSF
#
# Prep + envío (desde un nodo de login):
#   source ~/jkvirome/jm79/mosquito-virome-jock/sanger_scripts_jobs/config.sh
#   mkdir -p "$LOGS_DIR" "$SCRATCH"
#   find "$DATA_DIR" -type f -name "*.fastq.gz" | sort > "$SCRATCH/filelist.txt"
#   N=$(wc -l < "$SCRATCH/filelist.txt")
#
#   bsub -J "fastqc[1-$N]%10" \
#        -o "$LOGS_DIR/fastqc.%J.%I.log" -e "$LOGS_DIR/fastqc.%J.%I.err" \
#        -q normal -n 4 -M 2000 \
#        -R "select[mem>2000] rusage[mem=2000] span[hosts=1]" \
#        "./fastqc_worker.sh $SCRATCH/filelist.txt"
#
# SALIDA (se conserva): $RESULTS_DIR/fastqc_raw/  -> .html + .zip por muestra
# ---------------------------------------------------------------------------
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

filelist="${1:?Uso: fastqc_worker.sh <filelist.txt> [output_dir]}"
outdir="${2:-$RESULTS_DIR/fastqc_raw}"

threads="${LSB_DJOB_NUMPROC:-4}"
activate_env "$ENV_FASTQC"

idx="${LSB_JOBINDEX:?Este script debe enviarse como job array}"
file=$(sed -n "${idx}p" "$filelist")
[[ -n "$file" ]] || { echo "ERROR: linea $idx vacia en $filelist"; exit 1; }

echo "Array element ${idx}: ${file}  (host: $(hostname))"
mkdir -p "$outdir"

fastqc -t "$threads" "$file" -o "$outdir"

echo "FastQC terminado para: $file"
