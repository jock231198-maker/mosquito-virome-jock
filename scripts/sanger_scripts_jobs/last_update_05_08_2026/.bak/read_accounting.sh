#!/bin/bash
# ---------------------------------------------------------------------------
# read_accounting.sh — LA TABLA. Cuantas reads entran y cuantas sobreviven
#                      a cada etapa, por muestra.
#
# Es el control de calidad mas util de todo el pipeline: una muestra que
# pierde el 95% en el trimming, o que mapea el 99% al huesped, o que se queda
# con 3000 reads virales, se ve de un vistazo. Sin esta tabla esos problemas
# aparecen tres semanas despues, en el arbol filogenetico.
#
# NO recuenta los FASTQ (seria lentisimo): reaprovecha numeros que las
# herramientas YA calcularon:
#   - reads crudas      <- fastqc_data.txt dentro del .zip de FastQC
#   - tras trimming     <- log de LSF de Trimmomatic ("Input Read Pairs...")
#   - mapeadas/no map.  <- summary de bowtie2 que guarda bowtie_map_worker
#   - contigs           <- se cuentan del FASTA (barato)
#
# Uso:
#   ./read_accounting.sh
#   ./read_accounting.sh --csv        # separado por comas en vez de tabs
# ---------------------------------------------------------------------------
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

SEP=$'\t'; EXT="tsv"
[[ "${1:-}" == "--csv" ]] && { SEP=","; EXT="csv"; }

OUTDIR="$RESULTS_DIR/qc_control"
mkdir -p "$OUTDIR"
OUT="$OUTDIR/read_accounting.$EXT"

# --- helpers ----------------------------------------------------------------

# Reads crudas desde el .zip de FastQC (campo "Total Sequences")
reads_from_fastqc() {
  local sample="$1" zip
  zip=$(find "$RESULTS_DIR/fastqc_raw" -name "${sample}_R1*_fastqc.zip" 2>/dev/null | head -1)
  [[ -n "$zip" ]] || { echo "NA"; return; }
  # -p extrae a stdout; el fichero vive en <base>/fastqc_data.txt
  unzip -p "$zip" '*/fastqc_data.txt' 2>/dev/null \
    | awk -F'\t' '/^Total Sequences/ {print $2; exit}' \
    | tr -d '[:space:]' || echo "NA"
}

# Pares que sobreviven al trimming, desde el log de LSF de Trimmomatic.
# Trimmomatic imprime: "Input Read Pairs: 1000 Both Surviving: 950 (95.00%) ..."
reads_from_trimlog() {
  local sample="$1"
  grep -h -A2 "Trimming ${sample} " "$LOGS_DIR"/trim.*.log 2>/dev/null \
    | grep -m1 "Both Surviving" \
    | sed -n 's/.*Both Surviving: \([0-9]*\).*/\1/p'
}

# Del summary de bowtie2: total de pares y % de alineamiento global
bowtie_total()  { sed -n 's/^\([0-9]*\) reads; of these:.*/\1/p' "$1" 2>/dev/null | head -1; }
bowtie_rate()   { sed -n 's/^\([0-9.]*\)% overall alignment rate.*/\1/p' "$1" 2>/dev/null | head -1; }

# Pares no mapeados: los guarda extract_unmapped_worker
unmapped_pairs() {
  # OJO: dos lineas a proposito. `local a="$1" b="${a}"` NO funciona: bash
  # expande todos los argumentos de `local` antes de ejecutarlo.
  local sample="$1"
  local f="$RESULTS_DIR/unmapped_counts/${sample}_unmapped_pairs.tsv"
  [[ -f "$f" ]] && cut -f2 "$f" || echo "NA"
}

# Numero de contigs de un ensamblaje (contar '>' es barato)
contig_count() {
  local sample="$1" fa
  fa=$(find "$SCRATCH/spades" -path "*${sample}*" -name "contigs.fasta" 2>/dev/null | head -1)
  [[ -n "$fa" ]] && grep -c '^>' "$fa" || echo "NA"
}

# porcentaje seguro (evita division por cero y NA)
pct() {
  local num="$1" den="$2"
  [[ "$num" =~ ^[0-9]+$ && "$den" =~ ^[0-9]+$ && "$den" -gt 0 ]] || { echo "NA"; return; }
  awk -v n="$num" -v d="$den" 'BEGIN{printf "%.2f", 100*n/d}'
}

# --- muestras ---------------------------------------------------------------
mapfile -t SAMPLES < <(find "$DATA_DIR" -name "*_R1_*.fastq.gz" -printf '%f\n' 2>/dev/null \
                       | cut -d_ -f1-2 | sort -u)

if [[ ${#SAMPLES[@]} -eq 0 ]]; then
  echo "ERROR: no se detectaron muestras en $DATA_DIR"; exit 1
fi

# --- construir la tabla -----------------------------------------------------
{
printf '%s\n' "sample${SEP}raw_pairs${SEP}trimmed_pairs${SEP}pct_trim_survive${SEP}mapped_host_pct${SEP}unmapped_pairs${SEP}pct_unmapped${SEP}contigs"

for s in "${SAMPLES[@]}"; do
  raw=$(reads_from_fastqc "$s");           raw=${raw:-NA}
  trim=$(reads_from_trimlog "$s");         trim=${trim:-NA}
  bsum="$RESULTS_DIR/mapping_stats/${s}_bowtie2_summary.txt"
  rate=$(bowtie_rate "$bsum");             rate=${rate:-NA}
  unm=$(unmapped_pairs "$s");              unm=${unm:-NA}
  ctg=$(contig_count "$s");                ctg=${ctg:-NA}

  p_trim=$(pct "$trim" "$raw")
  p_unm=$(pct "$unm" "$trim")

  printf '%s\n' "${s}${SEP}${raw}${SEP}${trim}${SEP}${p_trim}${SEP}${rate}${SEP}${unm}${SEP}${p_unm}${SEP}${ctg}"
done
} > "$OUT"

# --- mostrar y avisar -------------------------------------------------------
echo "==========================================================="
echo "READ ACCOUNTING   $(date)"
echo "==========================================================="
column -t -s"$SEP" "$OUT" 2>/dev/null || cat "$OUT"
echo
echo "Tabla: $OUT"
echo

echo "== Avisos automaticos"
warned=0
while IFS="$SEP" read -r s raw trim p_trim rate unm p_unm ctg; do
  [[ "$s" == "sample" ]] && continue

  # Supervivencia al trimming baja -> adaptadores mal, o calidad mala
  if [[ "$p_trim" != "NA" ]] && awk -v v="$p_trim" 'BEGIN{exit !(v<70)}'; then
    echo "  [!] $s: solo ${p_trim}% sobrevive al trimming (revisa adaptadores/calidad)"
    warned=1
  fi
  # Muy pocas reads crudas -> muestra fallida en secuenciacion
  if [[ "$raw" =~ ^[0-9]+$ ]] && (( raw < 1000000 )); then
    echo "  [!] $s: solo $raw pares crudos (muestra con poca profundidad)"
    warned=1
  fi
  # Casi todo mapea al huesped -> poco material viral
  if [[ "$rate" != "NA" ]] && awk -v v="$rate" 'BEGIN{exit !(v>99)}'; then
    echo "  [!] $s: ${rate}% mapea al huesped (queda muy poco no-mapeado)"
    warned=1
  fi
  # Pocas reads no mapeadas -> ensamblaje sera pobre
  if [[ "$unm" =~ ^[0-9]+$ ]] && (( unm < 10000 )); then
    echo "  [!] $s: solo $unm pares no mapeados (ensamblaje poco fiable)"
    warned=1
  fi
  # Muy pocos contigs
  if [[ "$ctg" =~ ^[0-9]+$ ]] && (( ctg < 100 )); then
    echo "  [!] $s: solo $ctg contigs"
    warned=1
  fi
done < "$OUT"
[[ $warned -eq 0 ]] && echo "  ninguno: todas las muestras dentro de rango"
echo

echo "== Columnas NA"
echo "  Un NA no es un error: significa que ese paso aun no ha corrido, o que"
echo "  el fichero de donde se lee el numero no esta donde se espera."
