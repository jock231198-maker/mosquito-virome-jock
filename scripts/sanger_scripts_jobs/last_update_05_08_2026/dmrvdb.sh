#!/bin/bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
activate_env "diamond_2.2.6"
diamond blastx \
  -d "$REFS_DIR/diamond/rvdb.dmnd" \
  -q "$RESULTS_DIR/qc_control/votus" \
  -o "$RESULTS_DIR/diamond_votus_rvdb.tsv" \
  --threads "${LSB_DJOB_NUMPROC:-8}" \
  --max-target-seqs 5 --evalue 1e-5 \
  --outfmt 6 qseqid sseqid pident length evalue bitscore stitle
