#!/bin/bash
# ---------------------------------------------------------------------------
# bowtie_index.sh  — construir el índice bowtie2 del genoma huésped (one-shot)
#
#   bsub -o "logs/bt2idx.%J.log" -e "logs/bt2idx.%J.err" \
#        -q long -n 8 -M 16000 \
#        -R "select[mem>16000] rusage[mem=16000] span[hosts=1]" \
#        "./bowtie_index.sh refs/aedes_genome refs/bt2_index aedes_aegypti"
#
# Aedes aegypti (~1.3 Gbp) cabe en el índice normal. Si tu genoma supera
# ~4 Gbp añade --large-index a bowtie2-build.
# ---------------------------------------------------------------------------
set -euo pipefail
SECONDS=0

input_dir="${1:?Uso: bowtie_index.sh <input_dir_con_fna> <output_dir> <nombre>}"
outdir="${2:?Falta output_dir}"
result_from="${3:?Falta nombre}"

THREADS="${LSB_DJOB_NUMPROC:-8}"
CONDA_ENV="${CONDA_ENV:-bowtie2}"   # usa EXACTAMENTE el mismo nombre en todos los pasos

eval "$(conda shell.bash hook)"
conda activate "$CONDA_ENV"

mkdir -p "$outdir"   # (faltaba en tu versión original)

# bowtie2-build admite varios fasta separados por comas
fnas=$(ls "$input_dir"/*.fna | paste -sd, -)
echo "Indexando: $fnas"
bowtie2-build --threads "$THREADS" "$fnas" "$outdir/${result_from}_index"

dur=$SECONDS
echo "Índice listo en $outdir/${result_from}_index"
echo "Tiempo: $((dur/60)) min $((dur%60)) s."
