#!/bin/bash
# ---------------------------------------------------------------------------
# star_align_worker.sh — alinear pares al huésped con STAR (una muestra por elemento)
#
# OJO memoria: STAR carga el índice del genoma en RAM (~30 GB para Aedes) en
# CADA elemento del array. Limita la concurrencia (%4) y reserva ~40 GB.
#
# Prep + envío:
#   source .../config.sh
#   ls "$SCRATCH/trimmed/$RESULT_FROM"/*_R1_001_paired.fastq.gz \
#       | xargs -n1 basename | cut -d_ -f1-2 | sort -u > "$SCRATCH/samples.txt"
#   N=$(wc -l < "$SCRATCH/samples.txt")
#   bsub -J "star[1-$N]%4" \
#        -o "$LOGS_DIR/star.%J.%I.log" -e "$LOGS_DIR/star.%J.%I.err" \
#        -q normal -n 8 -M 40000 \
#        -R "select[mem>40000] rusage[mem=40000] span[hosts=1]" \
#        "./star_align_worker.sh $SCRATCH/samples.txt $SCRATCH/trimmed/$RESULT_FROM"
#
# SALIDA:
#   $SCRATCH/aligned/...  BAM + _STARtmp (temporales BORRABLES)
#   $RESULTS_DIR/star_stats/<sample>_Log.final.out   se conserva
#
# --outReadsUnmapped Fastx hace que STAR escriba directamente los no-mapeados,
# ahorrándote el paso extract_unmapped en la rama STAR.
# ---------------------------------------------------------------------------
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

samples="${1:?Uso: star_align_worker.sh <samples.txt> <input_dir> [genome_dir] [output_base] [result_from]}"
indir="${2:?Falta input_dir}"
genome_dir="${3:-$STAR_INDEX}"
outbase="${4:-$SCRATCH/aligned}"
result_from="${5:-$RESULT_FROM}"

THREADS="${LSB_DJOB_NUMPROC:-8}"
activate_env "$ENV_STAR"

idx="${LSB_JOBINDEX:?Enviar como job array}"
sample=$(sed -n "${idx}p" "$samples")
[[ -n "$sample" ]] || { echo "ERROR: linea $idx vacia en $samples"; exit 1; }

R1="$indir/${sample}_R1_001_paired.fastq.gz"
R2="$indir/${sample}_R2_001_paired.fastq.gz"
outdir="$outbase/$result_from/${sample}"
statsdir="$RESULTS_DIR/star_stats"
mkdir -p "$outdir" "$statsdir"

echo "STAR align $sample (host: $(hostname))"
if [[ -f "$R1" && -f "$R2" ]]; then
  STAR --runThreadN "$THREADS" \
       --runMode alignReads \
       --genomeDir "$genome_dir" \
       --readFilesIn "$R1" "$R2" \
       --readFilesCommand zcat \
       --outSAMtype BAM Unsorted \
       --outReadsUnmapped Fastx \
       --outFileNamePrefix "${outdir}/${sample}_"

  # Log final -> resultado permanente ; _STARtmp -> temporal, se borra ya
  cp "${outdir}/${sample}_Log.final.out" "$statsdir/${sample}_Log.final.out" 2>/dev/null || true
  rm -rf "${outdir}/${sample}__STARtmp"

  echo "Done: $sample"
else
  echo "ERROR: faltan R1/R2 para $sample"; exit 1
fi
