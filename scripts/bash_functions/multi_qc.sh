#!/bin/bash
SECONDS=0
eval "$(conda shell.bash hook)"
conda activate Multi_QC # Inicializa conda # Activar el entorno con MultiQC instalado
# Descripción: Ejecuta MultiQC sobre los resultados de FastQC
# Definir rutas de entrada y salida
input_dir="$1"   # Carpeta con resultados de FastQC
result_from="$2" # Carpeta donde se guardará el reporte de MultiQC
output_dir=$HOME/Desktop/multiqc/${result_from}
mkdir -p "$output_dir" # Crear carpeta de salida si no existe

multiqc "$input_dir"/*_paired_fastqc.zip -o "$output_dir" # Ejecutar MultiQC sobre los archivos FastQC (solo los *fastqc.zip)

conda deactivate # Desactivar el entorno conda
echo "MultiQC finalizado"
duration=$SECONDS
echo "Process finished in $(($duration / 60)) minutes and $(($duration % 60)) seconds."






