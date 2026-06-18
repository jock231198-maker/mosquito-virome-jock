#!/bin/bash
# ---------------------------------------------------------------------------
# cat_samples.sh
# Concatena los lanes L001 + L002 (R1 y R2) por muestra.
#
# Uso:
#   ./cat_samples.sh <input_dir> <output_dir>
#
# Es una tarea de I/O puro. Para pocas muestras puedes correrlo en un nodo
# de login. Para muchas, envíalo a LSF:
#   bsub -o logs/cat.%J.log -e logs/cat.%J.err -q normal -M 1000 \
#        -R "select[mem>1000] rusage[mem=1000]" \
#        "./cat_samples.sh raw_fastq cat_fastq"
# ---------------------------------------------------------------------------
set -euo pipefail

# Raíz por argumento -> el mismo script sirve en tu Mac y en la farm
indir="${1:?Falta input_dir}"
outdir="${2:?Falta output_dir}"

mkdir -p "$outdir"
cd "$indir" || { echo "No se pudo acceder a $indir"; exit 1; }

# Lista de muestras únicas a partir de los nombres de archivo
for sample in $(ls *_L001_R1_001.fastq.gz 2>/dev/null | cut -d'_' -f1-2 | sort -u); do
    echo "Concatenando $sample"

    r1_l001="${sample}_L001_R1_001.fastq.gz"
    r1_l002="${sample}_L002_R1_001.fastq.gz"
    r2_l001="${sample}_L001_R2_001.fastq.gz"
    r2_l002="${sample}_L002_R2_001.fastq.gz"

    # outdir es constante; las salidas se nombran por muestra
    r1_out="${outdir}/${sample}_R1_001.fastq.gz"
    r2_out="${outdir}/${sample}_R2_001.fastq.gz"

    # R1
    if [[ -f "$r1_l001" && -f "$r1_l002" ]]; then
        cat "$r1_l001" "$r1_l002" > "$r1_out"
        echo "  Creado: $r1_out"
    else
        echo "  WARNING: faltan archivos R1 para $sample"
    fi

    # R2  (ahora sí escribe en $outdir, como R1)
    if [[ -f "$r2_l001" && -f "$r2_l002" ]]; then
        cat "$r2_l001" "$r2_l002" > "$r2_out"
        echo "  Creado: $r2_out"
    else
        echo "  WARNING: faltan archivos R2 para $sample"
    fi
done

echo "Concatenación terminada. Resultados en: $outdir"
