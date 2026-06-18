#!/bin/bash
# ---------------------------------------------------------------------------
# star_align_worker.sh  — alinear pares al huésped con STAR (un sample por elemento)
#
# OJO memoria: STAR carga el índice del genoma en RAM (~30 GB para Aedes) en
# CADA elemento del array. Limita la concurrencia (%4) y reserva ~40 GB.
#
# Prep + envío:
#   ls trimmed/M57/*_R1_001_paired.fastq.gz | xargs -n1 basename | cut -d_ -f1-2 | sort -u > samples.txt
#   N=$(wc -l < samples.txt)
#   bsub -J "star[1-$N]%4" -o "logs/star.%J.%I.log" -e "logs/star.%J.%I.err" \
#        -q normal -n 8 -M 40000 \
#        -R "select[mem>40000] rusage[mem=40000] span[hosts=1]" \
#        "./star_align_worker.sh samples.txt trimmed/M57 refs/star_index aligned M57"
# ---------------------------------------------------------------------------
set -euo pipefail

samples="${1:?Uso: star_align_worker.sh <samples.txt> <input_dir> <genome_dir> <output_base> <result_from>}"
indir="${2:?Falta input_dir}"
genome_dir="${3:?Falta genome_dir}"
outbase="${4:?Falta output_base}"
result_from="${5:?Falta result_from}"

THREADS="${LSB_DJOB_NUMPROC:-8}"
CONDA_ENV="${CONDA_ENV:-star_env}"

eval "$(conda shell.bash hook)"
conda activate "$CONDA_ENV"

idx="${LSB_JOBINDEX:?Enviar como job array}"
sample=$(sed -n "${idx}p" "$samples")

R1="$indir/${sample}_R1_001_paired.fastq.gz"
R2="$indir/${sample}_R2_001_paired.fastq.gz"   # <- R2 distinto de R1 (bug corregido)
outdir="$outbase/$result_from/${sample}"
mkdir -p "$outdir"

echo "STAR align $sample (host: $(hostname))"
if [[ -f "$R1" && -f "$R2" ]]; then
  STAR --runThreadN "$THREADS" \
       --runMode alignReads \
       --genomeDir "$genome_dir" \
       --readFilesIn "$R1" "$R2" \
       --readFilesCommand zcat \
       --outSAMtype BAM Unsorted \
       --outFileNamePrefix "${outdir}/${sample}_"
  # Atajo opcional para virome: STAR puede escribir directamente los no-mapeados con
  #   --outReadsUnmapped Fastx
  # y te ahorras el paso extract_unmapped para la rama STAR.
  echo "Done: $sample"
else
  echo "ERROR: faltan R1/R2 para $sample"; exit 1
fi
