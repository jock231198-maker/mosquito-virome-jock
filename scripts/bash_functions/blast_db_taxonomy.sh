#!/bin/bash 
output_dir="$1" # Carpeta de salida para los resultados

mkdir -p "$output_dir"
cd "$output_dir"
# TaxID de Viruses en NCBI es 10239
# Descargar el archivo de taxonomía
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xvzf taxdump.tar.gz