#!/bin/bash
# ---------------------------------------------------------------------------
# run_multiqc.sh
# Agrega los reportes de FastQC (y de otras herramientas) en un único HTML.
#
# Uso:
#   ./run_multiqc.sh <input_dir> <output_dir> [nombre_reporte]
#
# En la farm conviene enviarlo a LSF (es ligero, 1 CPU, poca memoria):
#   bsub -o logs/multiqc.%J.log -e logs/multiqc.%J.err \
#        -q normal -M 4000 -R "select[mem>4000] rusage[mem=4000]" \
#        "./run_multiqc.sh fastqc_results multiqc_report"
# ---------------------------------------------------------------------------
set -euo pipefail
SECONDS=0

input_dir="${1:?Uso: run_multiqc.sh <input_dir> <output_dir> [nombre_reporte]}"
output_dir="${2:?Falta output_dir}"
report_name="${3:-multiqc_report}"

# MultiQC está en el mismo entorno conda que FastQC
eval "$(conda shell.bash hook)"
conda activate Fast_QC

mkdir -p "$output_dir"
echo "Ejecutando MultiQC sobre $input_dir ..."

# --force sobrescribe un reporte previo con el mismo nombre
multiqc "$input_dir" -o "$output_dir" -n "$report_name" --force

dur=$SECONDS
echo "MultiQC terminado. Reporte: $output_dir/${report_name}.html"
echo "Tiempo: $((dur / 60)) min $((dur % 60)) s."
