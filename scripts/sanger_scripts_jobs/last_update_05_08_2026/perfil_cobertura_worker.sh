#!/bin/bash
# perfil_cobertura_worker.sh — perfil de profundidad a lo largo del genoma de las
# vOTUs de un BED, para una muestra (un elemento del array de LSF).
#
# Ordena e indexa el BAM de vOTUs de la muestra (temporal en Lustre), saca la
# profundidad por posicion con samtools depth -a (incluye las posiciones a 0) y la
# resume en ventanas de W pb para que el fichero sea pequeno.
#
# Uso (lo lanza perfil_cobertura.sh):
#   perfil_cobertura_worker.sh <lista_muestras> <bed> <dir_salida> [ventana_pb]
#
# Salida: <dir_salida>/<muestra>.perfil.tsv, sin cabecera:
#   muestra  votu  inicio  fin  profundidad_media
#
# Idempotente: si la salida ya existe y no esta vacia, no hace nada.

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

lista="${1:?falta la lista de muestras}"
bed="${2:?falta el BED de vOTUs}"
outdir="${3:?falta el directorio de salida}"
W="${4:-50}"

: "${LSB_JOBINDEX:?este script necesita ejecutarse como array de LSF}"
sample="$(sed -n "${LSB_JOBINDEX}p" "$lista")"
[[ -n "$sample" ]] || { echo "ERROR: indice $LSB_JOBINDEX fuera de rango en $lista"; exit 1; }

bam="$SCRATCH/mapped_votus/votus/${sample}.bam"
out="$outdir/${sample}.perfil.tsv"
mkdir -p "$outdir"

if [[ -s "$out" ]]; then echo "[$sample] ya existe $out"; exit 0; fi
[[ -s "$bam" ]] || { echo "ERROR: no existe $bam"; exit 1; }
[[ -s "$bed" ]] || { echo "ERROR: no existe $bed"; exit 1; }

activate_env "${ENV_SAMTOOLS:-bowtie2_2.5.5}"
THREADS="${LSB_DJOB_NUMPROC:-4}"
export LC_ALL=C

tmp="$SCRATCH/tmp_perfil/${sample}_${LSB_JOBID:-manual}"
mkdir -p "$tmp"
trap 'rm -rf "$tmp"' EXIT

echo "[$sample] ordenando ($THREADS hilos)"
samtools sort -@ "$THREADS" -m 1G -T "$tmp/sort" -o "$tmp/s.bam" "$bam"
samtools index -@ "$THREADS" "$tmp/s.bam"

echo "[$sample] profundidad en $(wc -l < "$bed") vOTUs, ventanas de $W pb"
samtools depth -a -b "$bed" "$tmp/s.bam" \
| awk -F'\t' -v OFS='\t' -v S="$sample" -v W="$W" '
    { w = int(($2 - 1) / W)
      if ($1 != v || w != cw) {
        if (v != "") print S, v, cw * W + 1, cw * W + n, sum / n
        v = $1; cw = w; sum = 0; n = 0
      }
      sum += $3; n++ }
    END { if (v != "") print S, v, cw * W + 1, cw * W + n, sum / n }' > "$out.partial"

mv "$out.partial" "$out"
echo "[$sample] listo: $(cut -f2 "$out" | sort -u | wc -l) vOTUs, $(wc -l < "$out") ventanas"
