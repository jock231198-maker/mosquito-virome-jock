#!/bin/bash
# kraken_reads_run.sh — Kraken2 + Bracken sobre lecturas limpias, las 22 muestras
#
# Uso:
#   bsub ... kraken_reads_run.sh <lista_muestras> <dir_lecturas> <dir_salida>
#
# NO es un array a proposito. La base de Kraken2 ocupa decenas de GB y se carga
# entera en memoria; con un array de 22 elementos se leerian esos GB desde Lustre
# veintidos veces. En un solo trabajo secuencial la base queda en la cache de
# pagina del nodo tras la primera muestra y las siguientes arrancan al instante.
#
# Idempotente: salta las muestras que ya tengan .report.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

list="${1:?falta la lista de muestras}"
indir="${2:?falta el directorio de lecturas}"
outdir="${3:?falta el directorio de salida}"

DB="${KRAKEN_DB:-$REFS_DIR/kraken2/pluspf}"
READLEN="${READLEN:-150}"        # para Bracken; ajustar si las lecturas no son 150 bp
LEVEL="${LEVEL:-S}"              # nivel de Bracken: S especie, G genero, F familia
MINHITS="${MINHITS:-10}"         # lecturas minimas para que Bracken estime una especie

mkdir -p "$outdir"
activate_env "${ENV_KRAKEN:-kraken2}"

THREADS="${LSB_DJOB_NUMPROC:-16}"

[[ -s "$DB/hash.k2d" ]] || { echo "ERROR: no encuentro la base en $DB"; exit 1; }
echo "== base: $DB  ($(du -h "$DB/hash.k2d" | cut -f1))"
echo "== hilos: $THREADS   lecturas de $READLEN bp"
echo

brk="$DB/database${READLEN}mers.kmer_distrib"
[[ -s "$brk" ]] || echo "AVISO: no hay $brk — se saltara Bracken"

n=0
while read -r sample; do
  [[ -n "$sample" ]] || continue
  n=$((n+1))
  R1="$indir/${sample}_R1_001_paired.fastq.gz"
  R2="$indir/${sample}_R2_001_paired.fastq.gz"
  rep="$outdir/${sample}.report"

  if [[ -s "$rep" ]]; then
    echo "[$n] $sample — ya hecho, saltando"; continue
  fi
  if [[ ! -f "$R1" || ! -f "$R2" ]]; then
    echo "[$n] $sample — FALTAN R1/R2, saltando"; continue
  fi

  echo "[$n] $sample — clasificando  $(date +%H:%M:%S)"
  kraken2 \
    --db "$DB" \
    --threads "$THREADS" \
    --paired "$R1" "$R2" \
    --gzip-compressed \
    --report "$rep.partial" \
    --report-minimizer-data \
    --output - \
    > /dev/null

  mv "$rep.partial" "$rep"
  echo "     report: $(wc -l < "$rep") lineas"

  # ---- Bracken: reestima abundancia corrigiendo el sesgo de Kraken2 ----------
  # Kraken2 deja muchas lecturas asignadas a nodos altos del arbol (genero,
  # familia) porque su k-mero no distingue mas. Bracken las redistribuye a
  # especie usando las distribuciones de k-meros de la propia base.
  if [[ -s "$brk" ]]; then
    bracken -d "$DB" -i "$rep" \
            -o "$outdir/${sample}.bracken" \
            -w "$outdir/${sample}.bracken.report" \
            -r "$READLEN" -l "$LEVEL" -t "$MINHITS" \
      || echo "     AVISO: bracken fallo en $sample (sigue)"
  fi

  echo "     hecho  $(date +%H:%M:%S)"
done < "$list"

# ---- resumen ----------------------------------------------------------------
echo
echo "===== RESUMEN: % de lecturas clasificadas por muestra ====="
for rep in "$outdir"/*.report; do
  [[ -s "$rep" ]] || continue
  s=$(basename "$rep" .report)
  # linea 'unclassified' = taxon 0; el resto es lo clasificado
  awk -v S="$s" -F'\t' '
    $NF ~ /^[ ]*unclassified/ { u=$1 }
    $(NF-2)=="R" || $(NF-2)=="root" { r=$1 }
    END { printf "  %-10s sin clasificar %6.2f%%\n", S, u }' "$rep"
done

echo
echo "Salidas en $outdir:"
echo "  <sample>.report          arbol de Kraken2 (con datos de minimizers)"
echo "  <sample>.bracken         abundancia reestimada a nivel $LEVEL"
echo "  <sample>.bracken.report  el mismo arbol, corregido por Bracken"
