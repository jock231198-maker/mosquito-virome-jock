#!/bin/bash
# ---------------------------------------------------------------------------
# diamond_make_db.sh  — construir la base de datos DIAMOND (one-shot)
#
#   bsub -o "logs/dmnddb.%J.log" -e "logs/dmnddb.%J.err" \
#        -q normal -n 4 -M 8000 \
#        -R "select[mem>8000] rusage[mem=8000] span[hosts=1]" \
#        "./diamond_make_db.sh refs/diamond_src refs/diamond_db/viral_refseq_proteins"
# ---------------------------------------------------------------------------
set -euo pipefail

src_dir="${1:?Uso: diamond_make_db.sh <dir_con_faa.gz> <db_prefix_salida>}"
db_prefix="${2:?Falta el prefijo de salida de la base de datos}"

THREADS="${LSB_DJOB_NUMPROC:-4}"
CONDA_ENV="${CONDA_ENV:-diamond}"

eval "$(conda shell.bash hook)"
conda activate "$CONDA_ENV"

cd "$src_dir"
mkdir -p "$(dirname "$db_prefix")"

# Combinar los dos faa comprimidos
cat viral.1.protein.faa.gz viral.2.protein.faa.gz > viral_refseq.faa.gz

# Construir la .dmnd  (DIAMOND lee el .gz directamente: antes apuntabas a .faa
# sin comprimir, que no existía)
diamond makedb \
    --in viral_refseq.faa.gz \
    --db "$db_prefix" \
    --threads "$THREADS"

echo "Base de datos DIAMOND lista: ${db_prefix}.dmnd"
