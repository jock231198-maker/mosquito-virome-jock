#!/bin/bash
output_dir="$1" # Carpeta de salida para los resultados 

mkdir -p "$output_dir"
cd "$output_dir"

# Obtener la lista de archivos y descargarlos uno a uno
curl -s "https://ftp.ncbi.nlm.nih.gov/genbank/" \
  | grep -oE 'gbvrl\d+\.seq\.gz' \
  | sort -u \
  | while read f; do
      echo "Descargando $f..."
      curl -O "https://ftp.ncbi.nlm.nih.gov/genbank/$f"
    done