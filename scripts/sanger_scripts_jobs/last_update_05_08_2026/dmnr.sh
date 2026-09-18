#!/bin/bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
activate_env "diamond_2.2.6"
diamond blastx \
  -d "$REFS_DIR/diamond/nr.dmnd" \
  -q "$RESULTS_DIR/qc_control/votus" \
  -o "$RESULTS_DIR/diamond_votus_nr.tsv" \
  --threads "${LSB_DJOB_NUMPROC:-16}" \
  --max-target-seqs 5 --evalue 1e-5 --block-size 6 --index-chunks 1 \
  --outfmt 6 qseqid sseqid pident length evalue bitscore staxids sscinames sskingdoms stitle
