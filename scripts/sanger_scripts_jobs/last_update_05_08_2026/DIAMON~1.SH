#!/bin/bash
# ---------------------------------------------------------------------------
# diamond_blastx_worker.sh — DIAMOND blastx de contigs vs proteínas virales
#                            (un query por elemento)
#
# Prep + envío:
#   source .../config.sh
#   find "$SCRATCH/spades" -type f -name "contigs_1000bp.fasta" | sort > "$SCRATCH/query_list.txt"
#   N=$(wc -l < "$SCRATCH/query_list.txt")
#   bsub -J "dmnd[1-$N]%6" \
#        -o "$LOGS_DIR/dmnd.%J.%I.log" -e "$LOGS_DIR/dmnd.%J.%I.err" \
#        -q normal -n 8 -M 16000 \
#        -R "select[mem>16000] rusage[mem=16000] span[hosts=1]" \
#        "./diamond_blastx_worker.sh $SCRATCH/query_list.txt"
#
# SALIDA (se conserva): $RESULTS_DIR/diamond/<sample>_diamond_extended.tsv
# ---------------------------------------------------------------------------
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

query_list="${1:?Uso: diamond_blastx_worker.sh <query_list.txt> [db.dmnd] [output_base]}"
db="${2:-$DIAMOND_DB}"
outbase="${3:-$RESULTS_DIR/diamond}"

THREADS="${LSB_DJOB_NUMPROC:-8}"
activate_env "$ENV_DIAMOND"

idx="${LSB_JOBINDEX:?Enviar como job array}"
query=$(sed -n "${idx}p" "$query_list")
[[ -n "$query" ]] || { echo "ERROR: linea $idx vacia en $query_list"; exit 1; }
sample=$(basename "$(dirname "$query")")
mkdir -p "$outbase"
out="$outbase/${sample}_diamond_extended.tsv"

echo "DIAMOND blastx $sample (host: $(hostname))"
diamond blastx \
    -q "$query" \
    -d "$db" \
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
