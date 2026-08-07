#!/bin/bash
# ---------------------------------------------------------------------------
# diamond_blastx_worker.sh  — DIAMOND blastx de contigs vs proteínas virales
#                             (un query por elemento)
#
# Prep + envío:
#   find . -type f -name "contigs_1000bp.fasta" | sort > query_list.txt
#   N=$(wc -l < query_list.txt)
#   bsub -J "dmnd[1-$N]%6" -o "logs/dmnd.%J.%I.log" -e "logs/dmnd.%J.%I.err" \
#        -q normal -n 8 -M 16000 \
#        -R "select[mem>16000] rusage[mem=16000] span[hosts=1]" \
#        "./diamond_blastx_worker.sh query_list.txt refs/diamond_db/viral_refseq_proteins.dmnd diamond_out"
# ---------------------------------------------------------------------------
set -euo pipefail

query_list="${1:?Uso: diamond_blastx_worker.sh <query_list.txt> <db.dmnd> <output_base>}"
DB="${2:?Falta la base de datos .dmnd}"
outbase="${3:?Falta output_base}"

THREADS="${LSB_DJOB_NUMPROC:-8}"
CONDA_ENV="${CONDA_ENV:-diamond}"

eval "$(conda shell.bash hook)"
conda activate "$CONDA_ENV"

idx="${LSB_JOBINDEX:?Enviar como job array}"
query=$(sed -n "${idx}p" "$query_list")
sample=$(basename "$(dirname "$query")")
mkdir -p "$outbase"
out="$outbase/${sample}_diamond_extended.tsv"

echo "DIAMOND blastx $sample (host: $(hostname))"
diamond blastx \
    -q "$query" \
    -d "$DB" \
    -o "$out" \
    -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle qcovhsp scovhsp \
    -e 1e-5 \
    --more-sensitive \
    -p "$THREADS"

# Añadir cabecera de columnas
hdr="qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstitle\tqcovhsp\tscovhsp"
tmp=$(mktemp)
{ echo -e "$hdr"; cat "$out"; } > "$tmp" && mv "$tmp" "$out"

echo "Done: $out"
