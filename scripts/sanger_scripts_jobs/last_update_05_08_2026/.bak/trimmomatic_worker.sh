#!/bin/bash
# ---------------------------------------------------------------------------
# trimmomatic_worker.sh — recorte de adaptadores (una muestra por elemento)
#
# Prep + envío:
#   source .../config.sh
#   ls "$SCRATCH/cat_fastq"/*_R1_001.fastq.gz | xargs -n1 basename | cut -d_ -f1-2 \
#       | sort -u > "$SCRATCH/samples.txt"
#   N=$(wc -l < "$SCRATCH/samples.txt")
#   bsub -J "trim[1-$N]%10" \
#        -o "$LOGS_DIR/trim.%J.%I.log" -e "$LOGS_DIR/trim.%J.%I.err" \
#        -q normal -n 4 -M 20000 \
#        -R "select[mem>20000] rusage[mem=20000] span[hosts=1]" \
#        "./trimmomatic_worker.sh $SCRATCH/samples.txt $SCRATCH/cat_fastq"
#
# OJO: -Xmx (memoria Java) debe ser MENOR que el -M de LSF (deja margen al JVM).
#
# SALIDA (intermedia, en SCRATCH): $SCRATCH/trimmed/$RESULT_FROM/
#   *_paired.fastq.gz    -> se usan aguas abajo
#   *_unpaired.fastq.gz  -> BORRABLES (ver cleanup.sh)
# ---------------------------------------------------------------------------
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

samples="${1:?Uso: trimmomatic_worker.sh <samples.txt> <input_dir> [output_base] [result_from]}"
indir="${2:?Falta input_dir}"
outbase="${3:-$SCRATCH/trimmed}"
result_from="${4:-$RESULT_FROM}"

JAVA_MEM="${JAVA_MEM:-16g}"
THREADS="${LSB_DJOB_NUMPROC:-4}"

activate_env "$ENV_TRIMMOMATIC"

idx="${LSB_JOBINDEX:?Este script debe enviarse como job array}"
sample=$(sed -n "${idx}p" "$samples")
[[ -n "$sample" ]] || { echo "ERROR: linea $idx vacia en $samples"; exit 1; }

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
