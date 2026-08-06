#!/bin/bash
# ---------------------------------------------------------------------------
# cobra_worker.sh — extender contigs con COBRA (una muestra por ejecución)
#
# Como las entradas vienen de varias carpetas distintas, este toma los 4
# archivos de forma explícita (más robusto que adivinar el patrón de nombres).
# Para varias muestras a la vez usa cobra_array_worker.sh con un manifiesto.
#
#   source .../config.sh
#   bsub -o "$LOGS_DIR/cobra.%J.log" -e "$LOGS_DIR/cobra.%J.err" \
#        -q normal -n 8 -M 16000 \
#        -R "select[mem>16000] rusage[mem=16000] span[hosts=1]" \
#        "./cobra_worker.sh M1_S1 \
#            $SCRATCH/spades/M1_S1/viral_query_contigs.fasta \
#            $SCRATCH/spades/M1_S1/contigs.fasta \
#            $SCRATCH/spades_maps/M1_S1_sorted_mapped.bam \
#            $SCRATCH/spades_maps/M1_S1_coverage_for_cobra.tsv"
#
# SALIDA (se conserva): $RESULTS_DIR/cobra/<sample>/
# ---------------------------------------------------------------------------
set -euo pipefail
SECONDS=0

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

sample="${1:?Uso: cobra_worker.sh <sample> <query.fasta> <assembly.fasta> <map.bam> <coverage.tsv> [outdir]}"
query="${2:?Falta query}"
assembly="${3:?Falta assembly}"
bam="${4:?Falta bam}"
coverage="${5:?Falta coverage}"
outdir="${6:-$RESULTS_DIR/cobra/$sample}"

THREADS="${LSB_DJOB_NUMPROC:-8}"
activate_env "$ENV_COBRA"

mkdir -p "$(dirname "$outdir")"
echo "COBRA $sample (host: $(hostname))"

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

dur=$SECONDS
echo "COBRA terminado: $outdir  ($((dur/60))m $((dur%60))s)"
