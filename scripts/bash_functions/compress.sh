#!/bin/bash
# Descripción: concatena las lecturas de los lanes L001 y L002 para cada muestra (R1 y R2)
workdir="/Users/JK/Desktop/JKViromedata/Aa"
# Cambia al directorio de trabajo
cd "$workdir" || { echo "No se pudo acceder a $workdir"; exit 1; }

# Buscar muestras únicas basadas en el patrón de L001
for sample in $(ls *_L001_R1_001.fastq | cut -d'_' -f1-2 | sort -u); do
    echo "Concatenating sample: $sample"

    # Definir nombres de archivos
    r1_l001="${sample}_L001_R1_001.fastq"
    r1_l002="${sample}_L002_R1_001.fastq"
    r2_l001="${sample}_L001_R2_001.fastq"
    r2_l002="${sample}_L002_R2_001.fastq"
    r1_out="${sample}_R1.fastq"
    r2_out="${sample}_R2.fastq"

    # Concatenar R1 si ambos archivos existen
    if [[ -f "$r1_l001" && -f "$r1_l002" ]]; then
        cat "$r1_l001" "$r1_l002" > "$r1_out"
        echo "Created: $r1_out"
    else
        echo "WARNING: Missing files for $sample R1"
    fi

    # Concatenar R2 si ambos archivos existen
    if [[ -f "$r2_l001" && -f "$r2_l002" ]]; then
        cat "$r2_l001" "$r2_l002" > "$r2_out"
        echo "Created: $r2_out"
    else
        echo "WARNING: Missing files for $sample R2"
    fi
done
rm *_L001_*.fastq *_L002_*.fastq