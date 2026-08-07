#!/bin/bash
# ---------------------------------------------------------------------------
# quast_worker.sh  — QC del ensamblaje con QUAST (un ensamblaje por elemento)
#
# Prep + envío:
#   find spades -type f -name "contigs.fasta" | sort > assembly_list.txt
#   N=$(wc -l < assembly_list.txt)
#   bsub -J "quast[1-$N]%10" -o "logs/quast.%J.%I.log" -e "logs/quast.%J.%I.err" \
#        -q normal -n 8 -M 8000 \
#        -R "select[mem>8000] rusage[mem=8000] span[hosts=1]" \
#        "./quast_worker.sh assembly_list.txt quast_out"
#
# Si el nombre del ensamblaje es siempre "contigs.fasta", el nombre de muestra
# se toma de la carpeta que lo contiene.
# ---------------------------------------------------------------------------
set -euo pipefail

assembly_list="${1:?Uso: quast_worker.sh <assembly_list.txt> <output_base>}"
outbase="${2:?Falta output_base}"

THREADS="${LSB_DJOB_NUMPROC:-8}"
CONDA_ENV="${CONDA_ENV:-quast}"

eval "$(conda shell.bash hook)"
conda activate "$CONDA_ENV"

idx="${LSB_JOBINDEX:?Enviar como job array}"
fasta=$(sed -n "${idx}p" "$assembly_list")
sample=$(basename "$(dirname "$fasta")")   # nombre de la carpeta del sample
outdir="$outbase/$sample"
mkdir -p "$outdir"

echo "QUAST $sample (host: $(hostname))"
quast.py "$fasta" \
    -o "$outdir" \
    --rna-finding \
    --fast \
    --threads "$THREADS" \
    --min-contig 200
echo "Done: $outdir"
