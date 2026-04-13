#!/bin/bash
SECONDS=0
eval "$(conda shell.bash hook)"
conda activate Fast_QC # Inicializa conda # Activar el entorno con MultiQC instalado
input_dir="$1"  # Carpeta con archivos fastq.gz para FastQC
result_from="$2" # Nombre de la carpeta de resultados basado en la fuente de datos
threads="$3" # Número de hilos para FastQC

# Carpeta de salida para los resultados
output_dir=$HOME/Desktop/fast_qc/${result_from}
echo "Running FastQC on files in $input_dir..."
mkdir -p "$output_dir"
# Buscar solo en input_dir
find "$input_dir" -type f -name "*paired.fastq.gz" | while read -r file; do
    echo "Processing $file..."
    fastqc -t $threads "$file" -o "$output_dir"
done
echo "FastQC finished. Results are in $output_dir"
duration=$SECONDS
echo "Process finished in $(($duration / 60)) minutes and $(($duration % 60)) seconds."
