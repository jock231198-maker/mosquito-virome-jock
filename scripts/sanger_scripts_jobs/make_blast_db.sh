#!/bin/bash
# ---------------------------------------------------------------------------
# make_blast_db.sh  — construir base de datos BLAST de nucleótidos (one-shot)
#
#   bsub -o "logs/blastdb.%J.log" -e "logs/blastdb.%J.err" \
#        -q normal -n 2 -M 8000 \
#        -R "select[mem>8000] rusage[mem=8000] span[hosts=1]" \
#        "./make_blast_db.sh refs/viral_refseq_all.fasta refs/blast_db/viral_refseq_db"
# ---------------------------------------------------------------------------
set -euo pipefail

in_fasta="${1:?Uso: make_blast_db.sh <input.fasta> <db_prefix_salida>}"
db_prefix="${2:?Falta el prefijo de salida}"

CONDA_ENV="${CONDA_ENV:-blast}"

eval "$(conda shell.bash hook)"
conda activate "$CONDA_ENV"

mkdir -p "$(dirname "$db_prefix")"

makeblastdb \
  -in "$in_fasta" \
  -dbtype nucl \
  -title "Viral_RefSeq" \
  -out "$db_prefix" \
  -parse_seqids        # permite recuperar secuencias por ID más adelante

echo "Base de datos BLAST lista: $db_prefix"
