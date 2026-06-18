#!/bin/bash
# ---------------------------------------------------------------------------
# star_index.sh  — generar el índice STAR del genoma huésped (one-shot)
#
# MUCHA memoria: Aedes aegypti (~1.3 Gbp) necesita ~30-40 GB de RAM.
#   bsub -o "logs/staridx.%J.log" -e "logs/staridx.%J.err" \
#        -q hugemem -n 16 -M 45000 \
#        -R "select[mem>45000] rusage[mem=45000] span[hosts=1]" \
#        "./star_index.sh refs/aedes_genome/genome.fna refs/aedes_genome/genes.gtf refs/star_index"
# ---------------------------------------------------------------------------
set -euo pipefail
SECONDS=0

genome_fasta="${1:?Uso: star_index.sh <genome.fna> <genes.gtf> <genome_dir>}"
gtf_file="${2:?Falta gtf}"
genome_dir="${3:?Falta genome_dir}"

THREADS="${LSB_DJOB_NUMPROC:-16}"
CONDA_ENV="${CONDA_ENV:-star_env}"
SJDB_OVERHANG="${SJDB_OVERHANG:-99}"      # = (longitud de read) - 1
SA_NBASES="${SA_NBASES:-14}"              # 14 para genomas grandes; menor para virales

eval "$(conda shell.bash hook)"
conda activate "$CONDA_ENV"

mkdir -p "$genome_dir"
STAR --runThreadN "$THREADS" \
     --runMode genomeGenerate \
     --genomeDir "$genome_dir" \
     --genomeFastaFiles "$genome_fasta" \
     --sjdbGTFfile "$gtf_file" \
     --sjdbOverhang "$SJDB_OVERHANG" \
     --genomeSAindexNbases "$SA_NBASES"

dur=$SECONDS
echo "Índice STAR listo en $genome_dir"
echo "Tiempo: $((dur/60)) min $((dur%60)) s."
