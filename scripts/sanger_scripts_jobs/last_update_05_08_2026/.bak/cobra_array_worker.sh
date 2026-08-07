#!/bin/bash
# ---------------------------------------------------------------------------
# cobra_array_worker.sh — COBRA como job array gobernado por un MANIFIESTO
#
# Para entradas heterogéneas (rutas que no siguen un patrón regular) el array
# se alimenta de un manifiesto TSV: una muestra por línea, columnas separadas
# por TAB en este orden:
#
#   sample <TAB> query.fasta <TAB> assembly.fasta <TAB> map.bam <TAB> coverage.tsv <TAB> outdir
#
# Generar el manifiesto (ajusta al layout real):
#   source .../config.sh
#   for s in $(cat "$SCRATCH/samples.txt"); do
#     printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$s" \
#       "$SCRATCH/spades/$s/viral_query_contigs.fasta" \
#       "$SCRATCH/spades/$s/contigs.fasta" \
#       "$SCRATCH/spades_maps/${s}_sorted_mapped.bam" \
#       "$SCRATCH/spades_maps/${s}_coverage_for_cobra.tsv" \
#       "$RESULTS_DIR/cobra/$s"
#   done > "$SCRATCH/cobra_manifest.tsv"
#
# Envío:
#   N=$(grep -cve '^[[:space:]]*$' "$SCRATCH/cobra_manifest.tsv")
#   bsub -J "cobra[1-$N]%6" \
#        -o "$LOGS_DIR/cobra.%J.%I.log" -e "$LOGS_DIR/cobra.%J.%I.err" \
#        -q normal -n 8 -M 16000 \
#        -R "select[mem>16000] rusage[mem=16000] span[hosts=1]" \
#        "./cobra_array_worker.sh $SCRATCH/cobra_manifest.tsv"
# ---------------------------------------------------------------------------
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

manifest="${1:?Uso: cobra_array_worker.sh <manifest.tsv>}"

THREADS="${LSB_DJOB_NUMPROC:-8}"
activate_env "$ENV_COBRA"

idx="${LSB_JOBINDEX:?Enviar como job array}"

# Lee la línea N del manifiesto y reparte las columnas (separadas por TAB)
line=$(sed -n "${idx}p" "$manifest")
[[ -n "$line" ]] || { echo "ERROR: linea $idx vacia en $manifest"; exit 1; }
IFS=$'\t' read -r sample query assembly bam coverage outdir <<< "$line"
outdir="${outdir:-$RESULTS_DIR/cobra/$sample}"

echo "COBRA $sample (host: $(hostname))"
mkdir -p "$(dirname "$outdir")"

cobra-meta \
    -q "$query" \
    -f "$assembly" \
    -a metaspades \
    -mink 21 \
    -maxk 55 \
    -m "$bam" \
    -c "$coverage" \
    -t "$THREADS" \
    -o "$outdir"

echo "Done: $outdir"
