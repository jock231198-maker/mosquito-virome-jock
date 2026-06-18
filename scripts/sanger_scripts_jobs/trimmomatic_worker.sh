#!/bin/bash
# ---------------------------------------------------------------------------
# trimmomatic_worker.sh  — recorte de adaptadores (un par de muestras por elemento)
#
# Prep + envío (desde un nodo de login, en tu scratch):
#   ls cat_fastq/*_R1_001.fastq.gz | xargs -n1 basename | cut -d_ -f1-2 | sort -u > samples.txt
#   N=$(wc -l < samples.txt)
#   bsub -J "trim[1-$N]%10" -o "logs/trim.%J.%I.log" -e "logs/trim.%J.%I.err" \
#        -q normal -n 4 -M 20000 \
#        -R "select[mem>20000] rusage[mem=20000] span[hosts=1]" \
#        "./trimmomatic_worker.sh samples.txt cat_fastq trimmed M57"
#
# OJO: -Xmx (memoria Java) debe ser MENOR que el -M de LSF (deja margen al JVM).
# ---------------------------------------------------------------------------
set -euo pipefail

samples="${1:?Uso: trimmomatic_worker.sh <samples.txt> <input_dir> <output_base> <result_from>}"
indir="${2:?Falta input_dir}"
outbase="${3:?Falta output_base}"
result_from="${4:?Falta result_from}"

# Ruta COMPLETA al fasta de adaptadores (no relativa: en la farm hay que dar la ruta absoluta)
ADAPTERS="${ADAPTERS:-$HOME/jkvirome/refs/adapters/TruSeq3-PE-2.fa}"
JAVA_MEM="${JAVA_MEM:-16g}"
THREADS="${LSB_DJOB_NUMPROC:-4}"
CONDA_ENV="${CONDA_ENV:-trimmomatic}"

eval "$(conda shell.bash hook)"
conda activate "$CONDA_ENV"

idx="${LSB_JOBINDEX:?Este script debe enviarse como job array}"
sample=$(sed -n "${idx}p" "$samples")

R1="$indir/${sample}_R1_001.fastq.gz"
R2="$indir/${sample}_R2_001.fastq.gz"
outdir="$outbase/$result_from"
mkdir -p "$outdir"

echo "Trimming $sample (host: $(hostname))"
if [[ -f "$R1" && -f "$R2" ]]; then
  trimmomatic PE -Xmx${JAVA_MEM} -threads "$THREADS" -phred33 \
    "$R1" "$R2" \
    "$outdir/${sample}_R1_001_paired.fastq.gz" "$outdir/${sample}_R1_001_unpaired.fastq.gz" \
    "$outdir/${sample}_R2_001_paired.fastq.gz" "$outdir/${sample}_R2_001_unpaired.fastq.gz" \
    ILLUMINACLIP:"${ADAPTERS}":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:25 MINLEN:35
  echo "Done: $sample"
else
  echo "ERROR: faltan R1/R2 para $sample"; exit 1
fi
