#!/bin/bash
# ---------------------------------------------------------------------------
# assembly_spades.sh  — ensamblar lecturas no mapeadas con rnaviralSPAdes
#                       (un sample por elemento)
#
# Prep + envío:
#   ls unmapped_fastq/*_unmapped_R1.fastq.gz | xargs -n1 basename \
#     | sed 's/_unmapped_R1.fastq.gz//' | sort -u > samples.txt
#   N=$(wc -l < samples.txt)
#   bsub -J "spades[1-$N]%4" -o "logs/spades.%J.%I.log" -e "logs/spades.%J.%I.err" \
#        -q long -n 8 -M 64000 \
#        -R "select[mem>64000] rusage[mem=64000] span[hosts=1]" \
#        "./assembly_spades.sh samples.txt unmapped_fastq assembly"
#
# OJO: el -m de SPAdes (límite de memoria, en GB) debe ser MENOR que el -M de LSF.
# ---------------------------------------------------------------------------
set -euo pipefail

samples="${1:?Uso: assembly_spades.sh <samples.txt> <input_dir> <output_base>}"
indir="${2:?Falta input_dir}"
outbase="${3:?Falta output_base}"

THREADS="${LSB_DJOB_NUMPROC:-8}"
MEM_GB="${MEM_GB:-60}"            # < reserva LSF (64000 MB ~ 64 GB)
CONDA_ENV="${CONDA_ENV:-spades}"

eval "$(conda shell.bash hook)"
conda activate "$CONDA_ENV"

idx="${LSB_JOBINDEX:?Enviar como job array}"
sample=$(sed -n "${idx}p" "$samples")

R1="$indir/${sample}_unmapped_R1.fastq.gz"
R2="$indir/${sample}_unmapped_R2.fastq.gz"
outdir="$outbase/${sample}"
mkdir -p "$outdir"

echo "rnaviralSPAdes $sample (host: $(hostname))"
if [[ -f "$R1" && -f "$R2" ]]; then
  # Tus reads no-mapeados salen como R1/R2 separados, así que usamos -1/-2
  # (no --12, que es para reads intercalados en un solo archivo).
  rnaviralspades.py -1 "$R1" -2 "$R2" -t "$THREADS" -m "$MEM_GB" -o "$outdir"
  echo "Done: contigs en $outdir"
else
  echo "ERROR: faltan R1/R2 para $sample"; exit 1
fi
