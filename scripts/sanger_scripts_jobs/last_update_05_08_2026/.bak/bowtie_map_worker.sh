#!/bin/bash
# ---------------------------------------------------------------------------
# bowtie_map_worker.sh — mapear al huésped, salida BAM (una muestra por elemento)
#
# Prep + envío:
#   source .../config.sh
#   ls "$SCRATCH/trimmed/$RESULT_FROM"/*_R1_001_paired.fastq.gz \
#       | xargs -n1 basename | cut -d_ -f1-2 | sort -u > "$SCRATCH/samples.txt"
#   N=$(wc -l < "$SCRATCH/samples.txt")
#   bsub -J "btmap[1-$N]%10" \
#        -o "$LOGS_DIR/btmap.%J.%I.log" -e "$LOGS_DIR/btmap.%J.%I.err" \
#        -q normal -n 4 -M 8000 \
#        -R "select[mem>8000] rusage[mem=8000] span[hosts=1]" \
#        "./bowtie_map_worker.sh $SCRATCH/samples.txt $SCRATCH/trimmed/$RESULT_FROM"
#
# SALIDA:
#   $SCRATCH/mapped/$RESULT_FROM/<sample>.bam     intermedio (BORRABLE tras extraer unmapped)
#   $RESULTS_DIR/mapping_stats/<sample>.txt       se conserva (% de mapeo al huésped)
# ---------------------------------------------------------------------------
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

samples="${1:?Uso: bowtie_map_worker.sh <samples.txt> <input_dir> [index] [output_base] [result_from]}"
indir="${2:?Falta input_dir}"
index="${3:-$BT2_INDEX}"
outbase="${4:-$SCRATCH/mapped}"
result_from="${5:-$RESULT_FROM}"

THREADS="${LSB_DJOB_NUMPROC:-4}"
activate_env "$ENV_BOWTIE2"

idx="${LSB_JOBINDEX:?Enviar como job array}"
sample=$(sed -n "${idx}p" "$samples")
[[ -n "$sample" ]] || { echo "ERROR: linea $idx vacia en $samples"; exit 1; }

R1="$indir/${sample}_R1_001_paired.fastq.gz"
R2="$indir/${sample}_R2_001_paired.fastq.gz"
outdir="$outbase/$result_from"
statsdir="$RESULTS_DIR/mapping_stats"
mkdir -p "$outdir" "$statsdir"

echo "Mapping $sample (host: $(hostname))"
if [[ -f "$R1" && -f "$R2" ]]; then
  # El resumen de bowtie2 va a stderr -> lo guardamos como resultado permanente
  bowtie2 --threads "$THREADS" -x "$index" -1 "$R1" -2 "$R2" \
    2> "$statsdir/${sample}_bowtie2_summary.txt" \
    | samtools view -bS -@ "$THREADS" - > "$outdir/${sample}.bam"
  echo "Done: $outdir/${sample}.bam"
  echo "Stats: $statsdir/${sample}_bowtie2_summary.txt"
else
  echo "ERROR: faltan R1/R2 para $sample"; exit 1
fi
