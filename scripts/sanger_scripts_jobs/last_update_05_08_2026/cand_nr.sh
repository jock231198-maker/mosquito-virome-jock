#!/bin/bash
# cand_nr.sh — blastx de la vOTU candidata contra nr, sensibilidad maxima
#
# Uso:
#   bsub ... cand_nr.sh [fasta_consulta] [salida]
#
# Por defecto consulta la candidata k119_1227 y escribe en $RESULTS_DIR.
# Se lanza a cola porque `nr.dmnd` son 350 GB: cargarla lleva minutos y en el
# nodo de login el proceso se corta.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

query="${1:-$RESULTS_DIR/qc_control/candidata_k119_1227.fasta}"
out="${2:-$RESULTS_DIR/cand_nr.tsv}"

[[ -s "$query" ]] || { echo "ERROR: no existe o esta vacio $query"; exit 1; }
[[ -s "$REFS_DIR/diamond/nr.dmnd" ]] || { echo "ERROR: falta nr.dmnd"; exit 1; }

activate_env "${ENV_DIAMOND:-diamond_2.2.6}"
THREADS="${LSB_DJOB_NUMPROC:-8}"

echo "== consulta: $query  ($(grep -c '^>' "$query") secuencias)"
echo "== base:     nr.dmnd  ($(du -h "$REFS_DIR/diamond/nr.dmnd" | cut -f1))"
echo "== hilos:    $THREADS"
echo

diamond blastx \
  -d "$REFS_DIR/diamond/nr.dmnd" \
  -q "$query" \
  -o "$out.partial" \
  --threads "$THREADS" \
  --very-sensitive \
  --max-target-seqs 25 \
  --evalue 1e-3 \
  --block-size 6 --index-chunks 1 \
  --outfmt 6 qseqid pident length qcovhsp scovhsp evalue bitscore stitle

mv "$out.partial" "$out"

echo
echo "== resultado: $(wc -l < "$out") hits =="
if [[ -s "$out" ]]; then
  sort -k7 -rn "$out" | head -10 | cut -c1-170
else
  echo "  SIN HITS a evalue 1e-3"
fi
