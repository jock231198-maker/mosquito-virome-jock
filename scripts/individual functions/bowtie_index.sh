#!/bin/bash
SECONDS=0
eval "$(conda shell.bash hook)"
conda activate Bow # Inicializa conda # Activar el entorno con Bowtie2 instalado
input_dir="$1"  # Carpeta con archivos fna para construir el índice
result_from="$2" # Nombre de la o las especies a indexar 
output_dir="/Users/JK/Desktop/bowtie2/IndexBow/${result_from}" # Carpeta donde se guardará el índice
#Indexing
# Building a small index
#comandtobuild /path/to/genome/fna /path/to/index/directory
bowtie2-build --threads 8 "$input_dir"/*.fna "$output_dir/${result_from}_index"
conda deactivate
duration=$SECONDS
echo "Process finished in $(($duration / 60)) minutes and $(($duration % 60)) seconds."     

#bowtie2-build /Users/JK/Desktop/jkviromedata/aedes_ae_ref_genome/GCF_002204515.2_AaegL5.0_genomic.fna /Users/JK/Desktop/bowtie2/IndexBow/aedes_aegypti_smallindex
# Building a large index --- bowtie2-build --large-index example/reference/lambda_virus.fa example/index/lambda_virus
