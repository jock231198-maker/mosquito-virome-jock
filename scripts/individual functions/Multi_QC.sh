#!/bin/bash
# Descripción: Ejecuta MultiQC sobre los resultados de FastQC
# Definir rutas de entrada y salida
input_dir="$HOME/Desktop/FastQC_results"     # Carpeta con resultados de FastQC
output_dir="$HOME/Desktop/FastQC_results_MultiQC" # Carpeta donde se guardará el reporte de MultiQC
# Inicializa conda # Activar el entorno con MultiQC instalado
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate Multi_QC
multiqc "$input_dir"/*paired_fastqc.zip -o "$output_dir" # Ejecutar MultiQC sobre los archivos FastQC (solo los *paired_fastqc.zip)
conda deactivate # Desactivar el entorno conda
echo "MultiQC finalizado"
