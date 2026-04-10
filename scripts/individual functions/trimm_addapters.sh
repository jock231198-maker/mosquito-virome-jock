#!/bin/bash
SECONDS=0
eval "$(conda shell.bash hook)"
conda activate trimmomatic # Inicializa conda # Activar el entorno con Trimmomatic instalado
input_dir="$1"  # Carpeta con archivos fastq.gz para Trimmomatic
result_from="$2" # Nombre de la carpeta de resultados basado en la fuente de datos
output_dir=/Users/JK/Desktop/trimmomatic/adapter_trim/${result_from} # carpeta salida
mkdir -p "$output_dir" # Crear carpeta de salida si no existe
for R1 in "$input_dir"/*_R1_001.fastq.gz
do
BASENAME=$(basename "$R1" | cut -d'_' -f1-2)
R2="$input_dir/${BASENAME}_R2_001.fastq.gz"
if [[ -f "$R1" && -f "$R2" ]]; then
trimmomatic PE -threads 4 -phred33 \
-Xmx16g \
-threads 4 \
-phred33 \
"$R1" "$R2" \
"$output_dir/${BASENAME}_R1_001_paired.fastq.gz" \
"$output_dir/${BASENAME}_R1_001_unpaired.fastq.gz" \
"$output_dir/${BASENAME}_R2_001_paired.fastq.gz" \
"$output_dir/${BASENAME}_R2_001_unpaired.fastq.gz" \
ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:25 MINLEN:35 
else
echo "Archivo R1 y R2 faltante para $BASENAME. Saltando..."
fi
done
conda deactivate
duration=$SECONDS
echo "Process finished in $(($duration / 60)) minutes and $(($duration % 60)) seconds."
