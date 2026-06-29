#!/bin/bash
# ---------------------------------------------------------------------------
# cobra_array_worker.sh  — COBRA como job array gobernado por un MANIFIESTO
#
# Para entradas heterogéneas (rutas que no siguen un patrón regular) el array
# se alimenta de un manifiesto TSV: una muestra por línea, columnas separadas
# por TAB en este orden:
#
#   sample <TAB> query.fasta <TAB> assembly.fasta <TAB> map.bam <TAB> coverage.tsv <TAB> outdir
#
# Ejemplo de cómo crear el manifiesto (ajusta las rutas a tu layout):
#   printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
#     M1_S1 \
#     spades/M1_S1/viral_query_contigs.fasta \
#     spades/Assembly/M1_S1/contigs.fasta \
#     bowtie2/spades_maps/M1_S1_sorted_mapped.bam \
#     bowtie2/spades_maps/M1_S1_coverage_for_cobra.tsv \
#     spades/COBRA_M1_S1 \
#     >> cobra_manifest.tsv
#   # ...una línea por muestra (o genéralas con un bucle)
#
# Envío:
#   N=$(grep -cve '^[[:space:]]*$' cobra_manifest.tsv)   # líneas no vacías
#   bsub -J "cobra[1-$N]%6" -o "logs/cobra.%J.%I.log" -e "logs/cobra.%J.%I.err" \
#        -q normal -n 8 -M 16000 \
#        -R "select[mem>16000] rusage[mem=16000] span[hosts=1]" \
#        "./cobra_array_worker.sh cobra_manifest.tsv"
# ---------------------------------------------------------------------------
set -euo pipefail

manifest="${1:?Uso: cobra_array_worker.sh <manifest.tsv>}"

THREADS="${LSB_DJOB_NUMPROC:-8}"
CONDA_ENV="${CONDA_ENV:-cobra}"

eval "$(conda shell.bash hook)"
conda activate "$CONDA_ENV"

idx="${LSB_JOBINDEX:?Enviar como job array}"

# Lee la línea N del manifiesto y reparte las columnas (separadas por TAB)
line=$(sed -n "${idx}p" "$manifest")
IFS=$'\t' read -r sample query assembly bam coverage outdir <<< "$line"

echo "COBRA $sample (host: $(hostname))"
mkdir -p "$(dirname "$outdir")"

cobra-meta \
    -q "$query" \
    -f "$assembly" \
    -a metaspades \
    -mink 21 \
    -maxk 55 \
    -m "$bam" \
    -c "$coverage" \
    -t "$THREADS" \
    -o "$outdir"

echo "Done: $outdir"
