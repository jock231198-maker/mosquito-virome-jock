#!/bin/bash
# ---------------------------------------------------------------------------
# fast_qc.sh — FastQC en SERIE sobre una carpeta entera.
#
# Versión secuencial: cómoda para lotes pequeños o para probar en un nodo de
# login. Para el lote completo de MERIDA usa fastqc_worker.sh como job array.
#
# Uso:
#   ./fast_qc.sh [input_dir] [output_dir] [threads]
#   ./fast_qc.sh                      # usa $DATA_DIR -> $RESULTS_DIR/fastqc_raw
# ---------------------------------------------------------------------------
set -euo pipefail
SECONDS=0

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

input_dir="${1:-$DATA_DIR}"
output_dir="${2:-$RESULTS_DIR/fastqc_raw}"
threads="${3:-${LSB_DJOB_NUMPROC:-4}}"

activate_env "$ENV_FASTQC"

echo "FastQC sobre: $input_dir"
echo "Salida en:    $output_dir"
mkdir -p "$output_dir"

n=0
while IFS= read -r file; do
    echo "Processing $file ..."
    fastqc -t "$threads" "$file" -o "$output_dir"
    n=$((n + 1))
done < <(find "$input_dir" -type f -name "*.fastq.gz" | sort)

[[ $n -gt 0 ]] || { echo "ERROR: no se encontraron .fastq.gz en $input_dir"; exit 1; }

echo "FastQC terminado: $n archivos. Resultados en $output_dir"
dur=$SECONDS
echo "Tiempo: $((dur / 60)) min $((dur % 60)) s."
