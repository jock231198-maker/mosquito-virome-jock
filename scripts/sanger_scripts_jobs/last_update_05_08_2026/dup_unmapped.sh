#!/bin/bash
# ---------------------------------------------------------------------------
# dup_unmapped.sh — fraccion de secuencias unicas en unmapped_fastq
#
# QUE MIDE Y POR QUE
#   Es el conteo de la §9 del HANDOFF, empaquetado. Sobre las primeras N reads
#   de cada <sample>_unmapped_R1.fastq.gz cuenta cuantas secuencias distintas
#   hay. Una fraccion de unicas baja significa que unas pocas moleculas copan
#   la libreria; en MERIDA eso resulto ser rRNA 28S residual, y explicaba por
#   que las muestras con MAS material eran las que PEOR ensamblaban.
#
#   Referencia MERIDA (22 muestras, 400,000 reads):
#     19.5% M36_S17  27.0% M74_S21  29.9% M72_S20  35.0% M62_S18  37.0% M71_S19
#     ---- el resto entre 37.3% (M48_S5) y 68.1% (M79_S28) ----
#
#   Debajo de ~37% = zona mala. Por encima de ~60% = libreria sana.
#
# OJO CON LA COMPARACION ENTRE LOTES
#   Este numero depende de la profundidad y de la longitud de read. Las UVRI
#   vienen recortadas con otra herramienta y pueden ser mas cortas, lo que
#   sube artificialmente la colision de secuencias identicas y BAJA el % de
#   unicas sin que haya mas duplicacion biologica. El script imprime la
#   longitud media muestreada justo para que no te cuelen esa.
#
# Uso:
#   source .../env_poblaciones.sh
#   ./dup_unmapped.sh                        # $SCRATCH/unmapped_fastq
#   ./dup_unmapped.sh <dir> [n_reads]
#
# SALIDA: $RESULTS_DIR/qc_control/duplicacion_unmapped.tsv
#         $RESULTS_DIR/qc_control/top_secuencias_unmapped.fasta  (3 por muestra)
# ---------------------------------------------------------------------------
set -euo pipefail

: "${SCRATCH:?source env_poblaciones.sh primero}"
: "${RESULTS_DIR:?source env_poblaciones.sh primero}"

indir="${1:-$SCRATCH/unmapped_fastq}"
NREADS="${2:-400000}"

# OJO: la variable NO puede llamarse LINES.
# `LINES` es una variable reservada de bash: la reescribe con el numero de filas
# de la terminal. Con LINES=1600000 y una ventana de 48 filas, el `head -n
# "$LINES"` de abajo leia 48 lineas = 12 reads, y el script devolvia "100% de
# unicas" en TODAS las muestras sin quejarse. Lo delato el control: M36_S17 de
# MERIDA, que vale 19.5%, salia al 100%.
NLINES=$(( NREADS * 4 ))

outdir="$RESULTS_DIR/qc_control"
tsv="$outdir/duplicacion_unmapped.tsv"
fa="$outdir/top_secuencias_unmapped.fasta"
mkdir -p "$outdir"

# sort fuera de Lustre: son ficheros temporales pequenos y muchos accesos
TMP="${TMPDIR:-/tmp}/dup_$$"
mkdir -p "$TMP"
trap 'rm -rf "$TMP"' EXIT

shopt -s nullglob
files=("$indir"/*_unmapped_R1.fastq.gz)
shopt -u nullglob
(( ${#files[@]} )) || { echo "ERROR: 0 ficheros *_unmapped_R1.fastq.gz en $indir"; exit 1; }

printf 'sample\treads_muestreadas\tsecuencias_unicas\tpct_unicas\tlong_media\ttop1_frac\n' > "$tsv"
: > "$fa"

echo "Muestreando $NREADS reads por fichero en $indir"
echo

for f in "${files[@]}"; do
  sample=$(basename "$f" _unmapped_R1.fastq.gz)
  seqs="$TMP/$sample.seq"

  # El SIGPIPE de gzip al cerrarse head haria fallar el script con pipefail.
  # El `|| true` dentro del grupo lo absorbe sin ocultar un gzip corrupto de
  # verdad: eso lo caza el chequeo de "0 lineas" de abajo.
  { gzip -cd "$f" 2>/dev/null || true; } | head -n "$NLINES" | awk 'NR % 4 == 2' > "$seqs"

  n=$(wc -l < "$seqs")
  if (( n == 0 )); then
    printf '  %-16s SIN DATOS (fichero vacio o gzip corrupto)\n' "$sample"
    printf '%s\tNA\tNA\tNA\tNA\tNA\n' "$sample" >> "$tsv"
    continue
  fi
  # GUARDIA: si el muestreo sale mucho mas corto de lo pedido, el resultado no
  # vale y hay que verlo, no publicarlo. Un fichero puede tener legitimamente
  # menos reads que NREADS, pero entonces el aviso lo dice y no se calla.
  if (( n < NREADS )); then
    printf '  %-16s AVISO: solo %d reads muestreadas de %d pedidas\n' "$sample" "$n" "$NREADS"
    if (( n < NREADS / 10 )); then
      echo "         Con tan pocas reads el % de unicas no significa nada."
      echo "         Si el fichero SI tiene mas reads, algo trunca la tuberia."
      printf '%s\t%d\tNA\tNA\tNA\tNA\n' "$sample" "$n" >> "$tsv"
      rm -f "$seqs"; continue
    fi
  fi

  # Conteo de unicas y, de paso, las 3 mas frecuentes
  LC_ALL=C sort -T "$TMP" "$seqs" | uniq -c | LC_ALL=C sort -rn > "$TMP/$sample.cnt"
  u=$(wc -l < "$TMP/$sample.cnt")
  pct=$(awk -v u="$u" -v n="$n" 'BEGIN{printf "%.1f", 100*u/n}')
  lmed=$(awk '{s += length($0)} END{printf "%.0f", s/NR}' "$seqs")
  top1=$(awk -v n="$n" 'NR==1{printf "%.2f", 100*$1/n}' "$TMP/$sample.cnt")

  head -n 3 "$TMP/$sample.cnt" | awk -v s="$sample" '{printf ">%s_top%d_n%d\n%s\n", s, NR, $1, $2}' >> "$fa"

  printf '%s\t%d\t%d\t%s\t%s\t%s\n' "$sample" "$n" "$u" "$pct" "$lmed" "$top1" >> "$tsv"

  flag=""
  awk -v p="$pct" 'BEGIN{exit !(p < 37.3)}' && flag="  <-- zona S17-S21"
  printf '  %-16s %s%% unicas   (%s bp medios, top1 = %s%% de las reads)%s\n' \
         "$sample" "$pct" "$lmed" "$top1" "$flag"

  rm -f "$seqs" "$TMP/$sample.cnt"
done

echo
echo "TSV   : $tsv"
echo "FASTA : $fa   (entrada del test de SortMeRNA de la §11)"
