#!/bin/bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
activate_env "${ENV_DIAMOND:-diamond_2.2.6}"
diamond blastx \
  -d "$REFS_DIR/diamond/nr.dmnd" \
  -q "$RESULTS_DIR/qc_control/votus" \
  -o "$RESULTS_DIR/votus_nr_arbiter.tsv" \
  --threads "${LSB_DJOB_NUMPROC:-8}" \
  --very-sensitive --max-target-seqs 25 --evalue 1e-3 \
  --block-size 6 --index-chunks 1 \
  --outfmt 6 qseqid sseqid pident length qlen evalue bitscore staxids sscinames qcovhsp stitle
echo "hits: $(wc -l < "$RESULTS_DIR/votus_nr_arbiter.tsv")"
