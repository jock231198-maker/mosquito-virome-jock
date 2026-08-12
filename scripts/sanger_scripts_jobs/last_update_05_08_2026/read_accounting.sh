#!/bin/bash
# ---------------------------------------------------------------------------
# read_accounting.sh — LA TABLA. Cuantas reads entran y cuantas sobreviven
#                      a cada etapa, por muestra.
#
# Es el control de calidad mas util del pipeline: una muestra que pierde el 95%
# en el trimming, o que mapea el 99% al huesped, o que se queda con 3000 reads
# virales, se ve de un vistazo. Sin esta tabla esos problemas aparecen tres
# semanas despues, en el arbol filogenetico.
#
# NO recuenta los FASTQ (seria lentisimo): reaprovecha numeros que las
# herramientas YA calcularon.
#
#   raw_pairs      <- raw_pairs_baseline.tsv (verificado por compress.sh)
#                     y si no existe, de los .zip de FastQC
#   trimmed_pairs  <- log de Trimmomatic ("Both Surviving")
#   bt_pairs       <- summary de bowtie2 ("N reads; of these:")
#   host_pct       <- summary de bowtie2 ("overall alignment rate")
#   unmapped_pairs <- extract_unmapped_worker
#   contigs        <- se cuentan del FASTA (barato)
#
# Uso:
#   ./read_accounting.sh
#   ./read_accounting.sh --csv
#   ./read_accounting.sh --exclude M41_S4          # muestras fallidas
# ---------------------------------------------------------------------------
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

SEP=$'\t'; EXT="tsv"; EXCLUDE=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --csv)     SEP=","; EXT="csv" ;;
    --exclude) EXCLUDE="${2:?--exclude necesita una lista separada por comas}"; shift ;;
    -h|--help) sed -n '2,30p' "$0"; exit 0 ;;
    *) echo "Opcion desconocida: $1"; exit 1 ;;
  esac
  shift
done

OUTDIR="$RESULTS_DIR/qc_control"
mkdir -p "$OUTDIR"
OUT="$OUTDIR/read_accounting.$EXT"
BASELINE="$OUTDIR/raw_pairs_baseline.tsv"

# --- helpers ----------------------------------------------------------------

# Reads crudas del .zip de FastQC. SUMA todos los ficheros de la muestra: con
# lanes hay dos por read y quedarse con el primero da la mitad de la cifra.
# El patron '_*R1*' cubre los dos layouts (M11_S3_R1_001 y M11_S3_L001_R1_001).
reads_from_fastqc() {
  local sample="$1" total=0 n=0 v
  while IFS= read -r zip; do
    [[ -n "$zip" ]] || continue
    v=$(unzip -p "$zip" '*/fastqc_data.txt' 2>/dev/null \
        | awk -F'\t' '/^Total Sequences/ {print $2; exit}' | tr -d '[:space:]')
    if [[ "$v" =~ ^[0-9]+$ ]]; then total=$((total + v)); n=$((n + 1)); fi
  done < <(find "$RESULTS_DIR/fastqc_raw" -name "${sample}_*R1*_fastqc.zip" 2>/dev/null | sort)
  (( n > 0 )) && echo "$total" || echo ""
}

# Prioriza la linea base verificada por compress.sh sobre lo que diga FastQC.
reads_raw() {
  local sample="$1" v=""
  if [[ -f "$BASELINE" ]]; then
    v=$(awk -F'\t' -v s="$sample" '$1 == s {print $2; exit}' "$BASELINE")
  fi
  [[ -n "$v" ]] || v=$(reads_from_fastqc "$sample")
  echo "$v"
}

# Pares que sobreviven al trimming.
#
# Dos trampas resueltas aqui:
#  1. Trimmomatic escribe sus estadisticas a STDERR, asi que con -o/-e separados
#     acaban en trim.*.err y no en trim.*.log. Se miran los dos.
#  2. Entre el "Trimming <sample>" del worker y la linea "Both Surviving" hay
#     ~7 lineas de Trimmomatic. Cualquier ventana -A fija se queda corta; como
#     el array escribe un fichero por elemento, no hace falta ventana: se
#     localiza el fichero de esa muestra y se busca en el entero.
reads_from_trimlog() {
  local sample="$1" f e
  f=$(grep -l "Trimming ${sample} " "$LOGS_DIR"/trim.*.log 2>/dev/null \
      | xargs -r ls -t 2>/dev/null | head -1)
  [[ -n "$f" ]] || return
  e="${f%.log}.err"
  { cat "$f"; [[ -f "$e" ]] && cat "$e"; } 2>/dev/null \
    | sed -n 's/.*Both Surviving: \([0-9]*\).*/\1/p' | head -1
}

# Del summary de bowtie2: total de pares procesados y % de alineamiento global.
bowtie_total() { sed -n 's/^\([0-9]*\) reads; of these:.*/\1/p' "$1" 2>/dev/null | head -1; }
bowtie_rate()  { sed -n 's/^\([0-9.]*\)% overall alignment rate.*/\1/p' "$1" 2>/dev/null | head -1; }

# Pares no mapeados: los guarda extract_unmapped_worker.
# Tolera que el .tsv lleve cabecera: coge el primer valor numerico de la col 2.
unmapped_pairs() {
  local f="$RESULTS_DIR/unmapped_counts/${1}_unmapped_pairs.tsv"
  [[ -f "$f" ]] || return
  awk -F'\t' '$2 ~ /^[0-9]+$/ {print $2; exit}' "$f"
}

# Numero de contigs (contar '>' es barato). Ruta directa: un find con comodines
# confunde M1_S1 con M1_S10.
contig_count() {
  local fa="$SCRATCH/spades/$1/contigs.fasta"
  [[ -f "$fa" ]] || return
  grep -c '^>' "$fa"
}

# Porcentaje seguro (evita division por cero y NA).
pct() {
  local num="$1" den="$2"
  [[ "$num" =~ ^[0-9]+$ && "$den" =~ ^[0-9]+$ && "$den" -gt 0 ]] || { echo "NA"; return; }
  awk -v n="$num" -v d="$den" 'BEGIN {printf "%.2f", 100*n/d}'
}

# --- muestras ---------------------------------------------------------------
# Orden de preferencia: la lista que se paso a los job arrays, luego cat_fastq,
# luego los crudos. Asi la tabla cubre exactamente lo que se proceso.
SRC=""
if [[ -f "$SCRATCH/samples.txt" ]]; then
  mapfile -t SAMPLES < <(grep -ve '^[[:space:]]*$' "$SCRATCH/samples.txt" | sort -u)
  SRC="samples.txt"
elif [[ -d "$SCRATCH/cat_fastq" ]]; then
  mapfile -t SAMPLES < <(find "$SCRATCH/cat_fastq" -name "*_R1_*.fastq.gz" -printf '%f\n' 2>/dev/null \
                         | cut -d_ -f1-2 | sort -u)
  SRC="cat_fastq"
else
  mapfile -t SAMPLES < <(find "$DATA_DIR" -name "*_R1_*.fastq.gz" -printf '%f\n' 2>/dev/null \
                         | cut -d_ -f1-2 | sort -u)
  SRC="DATA_DIR"
fi

if [[ ${#SAMPLES[@]} -eq 0 ]]; then
  echo "ERROR: no se detectaron muestras (fuente intentada: $SRC)"; exit 1
fi

# Exclusiones
if [[ -n "$EXCLUDE" ]]; then
  IFS=',' read -ra EXC <<< "$EXCLUDE"
  KEEP=()
  for s in "${SAMPLES[@]}"; do
    skip=0
    for x in "${EXC[@]}"; do [[ "$s" == "$x" ]] && skip=1; done
    (( skip )) || KEEP+=("$s")
  done
  SAMPLES=("${KEEP[@]}")
fi

# --- construir la tabla -----------------------------------------------------
{
printf '%s\n' "sample${SEP}raw_pairs${SEP}trimmed_pairs${SEP}pct_trim${SEP}bt_pairs${SEP}host_pct${SEP}unmapped_pairs${SEP}pct_unmapped${SEP}contigs"

for s in "${SAMPLES[@]}"; do
  raw=$(reads_raw "$s");                   raw=${raw:-NA}
  trim=$(reads_from_trimlog "$s");         trim=${trim:-NA}
  bsum="$RESULTS_DIR/mapping_stats/${s}_bowtie2_summary.txt"
  btot=$(bowtie_total "$bsum");            btot=${btot:-NA}
  rate=$(bowtie_rate "$bsum");             rate=${rate:-NA}
  unm=$(unmapped_pairs "$s");              unm=${unm:-NA}
  ctg=$(contig_count "$s");                ctg=${ctg:-NA}

  p_trim=$(pct "$trim" "$raw")
  p_unm=$(pct "$unm" "$trim")

  printf '%s\n' "${s}${SEP}${raw}${SEP}${trim}${SEP}${p_trim}${SEP}${btot}${SEP}${rate}${SEP}${unm}${SEP}${p_unm}${SEP}${ctg}"
done
} > "$OUT"

# --- mostrar y avisar -------------------------------------------------------
echo "==========================================================="
echo "READ ACCOUNTING   $(date '+%Y-%m-%d %H:%M')"
echo "==========================================================="
echo "muestras : ${#SAMPLES[@]}  (fuente: $SRC)"
if [[ -f "$BASELINE" ]]; then
  echo "raw_pairs: $BASELINE  (verificado por compress.sh)"
else
  echo "raw_pairs: .zip de FastQC, sumando lanes  (no hay raw_pairs_baseline.tsv)"
fi
[[ -n "$EXCLUDE" ]] && echo "excluidas: $EXCLUDE"
echo
column -t -s"$SEP" "$OUT" 2>/dev/null || cat "$OUT"
echo
echo "Tabla: $OUT"
echo

echo "== Avisos automaticos"
warned=0
while IFS="$SEP" read -r s raw trim p_trim btot rate unm p_unm ctg; do
  [[ "$s" == "sample" ]] && continue

  # Supervivencia al trimming baja -> adaptadores mal, o calidad mala
  if [[ "$p_trim" != "NA" ]] && awk -v v="$p_trim" 'BEGIN {exit !(v < 70)}'; then
    echo "  [!] $s: solo ${p_trim}% sobrevive al trimming (revisa adaptadores/calidad)"; warned=1
  fi
  # Muy pocas reads crudas -> muestra fallida en secuenciacion
  if [[ "$raw" =~ ^[0-9]+$ ]] && (( raw < 1000000 )); then
    echo "  [!] $s: solo $raw pares crudos (muestra fallida?)"; warned=1
  fi
  # Descuadre entre lo que salio del trimming y lo que vio bowtie2
  if [[ "$trim" =~ ^[0-9]+$ && "$btot" =~ ^[0-9]+$ ]] && (( trim != btot )); then
    echo "  [!] $s: trimming dio $trim pares pero bowtie2 proceso $btot (descuadre)"; warned=1
  fi
  # Casi todo mapea al huesped -> poco material viral
  if [[ "$rate" != "NA" ]] && awk -v v="$rate" 'BEGIN {exit !(v > 99)}'; then
    echo "  [!] $s: ${rate}% mapea al huesped (queda muy poco no-mapeado)"; warned=1
  fi
  # Pocas reads no mapeadas -> ensamblaje sera pobre
  if [[ "$unm" =~ ^[0-9]+$ ]] && (( unm < 10000 )); then
    echo "  [!] $s: solo $unm pares no mapeados (ensamblaje poco fiable)"; warned=1
  fi
  # Muy pocos contigs
  if [[ "$ctg" =~ ^[0-9]+$ ]] && (( ctg < 100 )); then
    echo "  [!] $s: solo $ctg contigs"; warned=1
  fi
done < "$OUT"
[[ $warned -eq 0 ]] && echo "  ninguno: todas las muestras dentro de rango"
echo

# Carpetas que nadie crea y que este script necesita leer
[[ -d "$RESULTS_DIR/unmapped_counts" ]] || \
  echo "  [nota] no existe $RESULTS_DIR/unmapped_counts: la columna unmapped_pairs saldra NA"

echo "== Columnas NA"
echo "  Un NA no es un error: significa que ese paso aun no ha corrido, o que"
echo "  el fichero de donde se lee el numero no esta donde se espera."