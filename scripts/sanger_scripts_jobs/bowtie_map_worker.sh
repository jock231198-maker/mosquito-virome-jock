#!/bin/bash
# ---------------------------------------------------------------------------
# bowtie_map_worker.sh  — mapear al huésped, salida BAM (un sample por elemento)
#
# Prep + envío:
#   ls trimmed/M57/*_R1_001_paired.fastq.gz | xargs -n1 basename | cut -d_ -f1-2 | sort -u > samples.txt
#   N=$(wc -l < samples.txt)
#   bsub -J "btmap[1-$N]%10" -o "logs/btmap.%J.%I.log" -e "logs/btmap.%J.%I.err" \
#        -q normal -n 4 -M 8000 \
#        -R "select[mem>8000] rusage[mem=8000] span[hosts=1]" \
#        "./bowtie_map_worker.sh samples.txt trimmed/M57 refs/bt2_index/aedes_aegypti_index mapped M57"
# ---------------------------------------------------------------------------
set -euo pipefail

samples="${1:?Uso: bowtie_map_worker.sh <samples.txt> <input_dir> <index_prefix> <output_base> <result_from>}"
indir="${2:?Falta input_dir}"
index="${3:?Falta index_prefix}"
outbase="${4:?Falta output_base}"
result_from="${5:?Falta result_from}"

THREADS="${LSB_DJOB_NUMPROC:-4}"
CONDA_ENV="${CONDA_ENV:-bowtie2}"   # debe incluir bowtie2 Y samtools

eval "$(conda shell.bash hook)"
conda activate "$CONDA_ENV"

idx="${LSB_JOBINDEX:?Enviar como job array}"
sample=$(sed -n "${idx}p" "$samples")

R1="$indir/${sample}_R1_001_paired.fastq.gz"
R2="$indir/${sample}_R2_001_paired.fastq.gz"
outdir="$outbase/$result_from"
mkdir -p "$outdir"

echo "Mapping $sample (host: $(hostname))"
if [[ -f "$R1" && -f "$R2" ]]; then
  bowtie2 --threads "$THREADS" -x "$index" -1 "$R1" -2 "$R2" \
    | samtools view -bS -@ "$THREADS" - > "$outdir/${sample}.bam"
  echo "Done: $outdir/${sample}.bam"
else
  echo "ERROR: faltan R1/R2 para $sample"; exit 1
fi
