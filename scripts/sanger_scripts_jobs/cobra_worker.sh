#!/bin/bash
# ---------------------------------------------------------------------------
# cobra_worker.sh  — extender contigs con COBRA (una muestra por ejecución)
#
# Como las entradas vienen de varias carpetas distintas, este toma los 4
# archivos de forma explícita (más robusto que adivinar el patrón de nombres).
#
#   bsub -o "logs/cobra.%J.log" -e "logs/cobra.%J.err" \
#        -q normal -n 8 -M 16000 \
#        -R "select[mem>16000] rusage[mem=16000] span[hosts=1]" \
#        "./cobra_worker.sh M1_S1 \
#            spades/viral_query_contigs.fasta \
#            spades/Assembly/.../M1_S1/contigs.fasta \
#            bowtie2/spades_maps/m1_sorted_mapped.bam \
#            bowtie2/spades_maps/coverage_for_cobra.tsv \
#            spades/COBRA_M1_S1"
# ---------------------------------------------------------------------------
set -euo pipefail
SECONDS=0

sample="${1:?Uso: cobra_worker.sh <sample> <query.fasta> <assembly.fasta> <map.bam> <coverage.tsv> <outdir>}"
query="${2:?Falta query}"
assembly="${3:?Falta assembly}"
bam="${4:?Falta bam}"
coverage="${5:?Falta coverage}"
outdir="${6:?Falta outdir}"

THREADS="${LSB_DJOB_NUMPROC:-8}"
CONDA_ENV="${CONDA_ENV:-cobra}"

eval "$(conda shell.bash hook)"
conda activate "$CONDA_ENV"

mkdir -p "$(dirname "$outdir")"
echo "COBRA $sample (host: $(hostname))"

# Continuaciones de línea corregidas (sin el "\ #" que rompía el comando)
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
