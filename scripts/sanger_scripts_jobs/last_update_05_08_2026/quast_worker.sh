#!/bin/bash
# ---------------------------------------------------------------------------
# quast_worker.sh — QC del ensamblaje con QUAST (un ensamblaje por elemento)
#
# Prep + envío:
#   source .../config.sh
#   find "$SCRATCH/spades" -type f -name "contigs.fasta" | sort > "$SCRATCH/assembly_list.txt"
#   N=$(wc -l < "$SCRATCH/assembly_list.txt")
#   bsub -J "quast[1-$N]%10" \
#        -o "$LOGS_DIR/quast.%J.%I.log" -e "$LOGS_DIR/quast.%J.%I.err" \
#        -q normal -n 8 -M 8000 \
#        -R "select[mem>8000] rusage[mem=8000] span[hosts=1]" \
#        "./quast_worker.sh $SCRATCH/assembly_list.txt"
#
# El nombre de muestra sale de la carpeta que contiene el contigs.fasta.
# SALIDA (se conserva): $RESULTS_DIR/quast/<sample>/report.html + report.tsv
# ---------------------------------------------------------------------------
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

assembly_list="${1:?Uso: quast_worker.sh <assembly_list.txt> [output_base]}"
outbase="${2:-$RESULTS_DIR/quast}"

THREADS="${LSB_DJOB_NUMPROC:-8}"
activate_env "$ENV_QUAST"

idx="${LSB_JOBINDEX:?Enviar como job array}"
fasta=$(sed -n "${idx}p" "$assembly_list")
[[ -n "$fasta" ]] || { echo "ERROR: linea $idx vacia en $assembly_list"; exit 1; }
sample=$(basename "$(dirname "$fasta")")   # nombre de la carpeta del sample
outdir="$outbase/$sample"
mkdir -p "$outdir"

echo "QUAST $sample (host: $(hostname))"
quast.py "$fasta" \
    -o "$outdir" \
    --rna-finding \
    --fast \
    --threads "$THREADS" \
    --min-contig 200
echo "Done: $outdir"
