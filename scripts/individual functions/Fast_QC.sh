#!/bin/bash
RUN FASTQC
# Carpeta donde están tus FASTQ
input_dir="$HOME/Desktop/SeqwithTrimmAdaptersdata"
# Carpeta de salida para los resultados
output_dir="$HOME/Desktop/FastQC_results"
echo "Running FastQC on files in $input_dir..."
mkdir -p "$output_dir"
# Buscar solo en input_dir
find "$input_dir" -type f -name "*paired.fastq.gz" | while read -r file; do
    echo "Processing $file..."
    fastqc "$file" -o "$output_dir"
done
echo "FastQC finished. Results are in $output_dir"