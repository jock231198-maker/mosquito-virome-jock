#!/bin/bash
# cov_worker.sh — amplitud de cobertura por vOTU, una muestra por elemento del array
#
# Uso:
#   bsub -J "cov[1-22]%8" ... cov_worker.sh <lista_muestras> <dir_bam> <dir_salida>
#
# Salida: <dir_salida>/<sample>.cov  con columnas
#   muestra  votu  longitud  bases_cubiertas  %cubierto  profundidad_media
#
# Idempotente: si el .cov ya existe y no está vacío, no rehace nada.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

list="${1:?falta la lista de muestras}"
bamdir="${2:?falta el directorio de BAM}"
outdir="${3:?falta el directorio de salida}"

: "${LSB_JOBINDEX:?este script necesita ejecutarse como array de LSF}"
sample="$(sed -n "${LSB_JOBINDEX}p" "$list")"
[[ -n "$sample" ]] || { echo "ERROR: indice $LSB_JOBINDEX fuera de rango en $list"; exit 1; }

bam="$bamdir/${sample}.bam"
out="$outdir/${sample}.cov"
mkdir -p "$outdir"

if [[ -s "$out" ]]; then
  echo "[$sample] ya existe $out — nada que hacer"; exit 0
fi
[[ -s "$bam" ]] || { echo "ERROR: no existe $bam"; exit 1; }

activate_env "${ENV_SAMTOOLS:-bowtie2_2.5.5}"
THREADS="${LSB_DJOB_NUMPROC:-4}"

# temporal en lustre, no en /tmp del nodo: los BAM rondan los GB
tmp="$SCRATCH/tmp_cov/${sample}_${LSB_JOBID:-manual}"
mkdir -p "$tmp"
trap 'rm -rf "$tmp"' EXIT

echo "[$sample] ordenando ($THREADS hilos)"
samtools sort -@ "$THREADS" -m 1G -T "$tmp/sort" -o "$tmp/${sample}.s.bam" "$bam"
samtools index -@ "$THREADS" "$tmp/${sample}.s.bam"

echo "[$sample] calculando cobertura"
samtools coverage "$tmp/${sample}.s.bam" \
  | awk -v S="$sample" 'NR>1{print S"\t"$1"\t"$3"\t"$5"\t"$6"\t"$7}' > "$out.partial"

mv "$out.partial" "$out"
echo "[$sample] listo: $(wc -l < "$out") vOTUs"
