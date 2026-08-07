#!/bin/bash
# ---------------------------------------------------------------------------
# multi_qc.sh — agrega los .zip de FastQC en un único reporte HTML.
#
# Uso:
#   ./multi_qc.sh [input_dir] [output_dir] [nombre_reporte]
#   ./multi_qc.sh                      # $RESULTS_DIR/fastqc_raw -> $RESULTS_DIR/multiqc
#
# Es ligero; se puede correr en un nodo de login o enviarlo con:
#   bsub -o "$LOGS_DIR/multiqc.%J.log" -e "$LOGS_DIR/multiqc.%J.err" \
#        -q normal -M 4000 -R "select[mem>4000] rusage[mem=4000]" \
#        "./multi_qc.sh"
#
# SALIDA (se conserva): $RESULTS_DIR/multiqc/<nombre_reporte>.html
# ---------------------------------------------------------------------------
set -euo pipefail
SECONDS=0

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

input_dir="${1:-$RESULTS_DIR/fastqc_raw}"
output_dir="${2:-$RESULTS_DIR/multiqc}"
report_name="${3:-multiqc_${RESULT_FROM}}"

activate_env "$ENV_MULTIQC"

mkdir -p "$output_dir"
echo "MultiQC sobre: $input_dir"

# --force sobrescribe un reporte previo con el mismo nombre
multiqc "$input_dir" -o "$output_dir" -n "$report_name" --force

dur=$SECONDS
echo "MultiQC terminado. Reporte: $output_dir/${report_name}.html"
echo "Tiempo: $((dur / 60)) min $((dur % 60)) s."
