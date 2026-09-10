#!/bin/bash
# ---------------------------------------------------------------------------
# collect_recovery.sh — reúne el recuento de lecturas de todas las etapas en
# UN solo fichero largo, listo para graficar en R.
#
# No recuenta FASTQ. Reaprovecha lo que las herramientas ya escribieron:
# el fastqc_data.txt dentro de los .zip, el "Both Surviving" de Trimmomatic,
# el summary de bowtie2, el Log.final.out de STAR y los .tsv de no-mapeados.
#
# SALIDAS  ($RESULTS_DIR/qc_control/)
#   read_recovery_long.tsv     sample · mapper · step · step_order · pairs
#   read_recovery_wide.tsv     una fila por muestra, todas las columnas crudas
#   read_recovery_metrics.tsv  tasas por muestra y mapeador
#
# USO
#   ./collect_recovery.sh --check     # qué encuentra y qué falta. EMPEZAR AQUI.
#   ./collect_recovery.sh             # escribe los tres TSV
#
# Rutas que puedes forzar por variable de entorno si no las adivina:
#   FQ_RAW_DIR  FQ_NOPOLYG_DIR  FQ_TRIMMED_DIR  STAR_DIR  BT_STATS_DIR
#   UNMAPPED_COUNTS_DIR  SAMPLES
# ---------------------------------------------------------------------------
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

CHECK=0
[[ "${1:-}" == "--check" ]] && CHECK=1

OUT="$RESULTS_DIR/qc_control"
mkdir -p "$OUT"

LONG="$OUT/read_recovery_long.tsv"
WIDE="$OUT/read_recovery_wide.tsv"
METRICS="$OUT/read_recovery_metrics.tsv"

# --- localización de fuentes -----------------------------------------------
# Cada etapa acepta varias carpetas candidatas; se usa la primera que exista.
first_dir() {
  local d
  for d in "$@"; do [[ -d "$d" ]] && { echo "$d"; return; }; done
}

FQ_RAW_DIR="${FQ_RAW_DIR:-$(first_dir "$RESULTS_DIR/fastqc_raw" "$RESULTS_DIR/fastqc/raw")}"
FQ_NOPOLYG_DIR="${FQ_NOPOLYG_DIR:-$(first_dir "$RESULTS_DIR/fastqc_nopolyg" "$RESULTS_DIR/fastqc/nopolyg" "$RESULTS_DIR/fastqc_polyg")}"
FQ_TRIMMED_DIR="${FQ_TRIMMED_DIR:-$(first_dir "$RESULTS_DIR/fastqc_trimmed" "$RESULTS_DIR/fastqc/trimmed")}"
BT_STATS_DIR="${BT_STATS_DIR:-$(first_dir "$RESULTS_DIR/mapping_stats")}"
UNMAPPED_COUNTS_DIR="${UNMAPPED_COUNTS_DIR:-$(first_dir "$RESULTS_DIR/unmapped_counts")}"
STAR_DIR="${STAR_DIR:-$(first_dir "$SCRATCH/aligned" "$SCRATCH/star" "$SCRATCH/aligned_star")}"
BASELINE="$RESULTS_DIR/qc_control/raw_pairs_baseline.tsv"

# --- lista de muestras ------------------------------------------------------
SAMPLES="${SAMPLES:-}"
if [[ -z "$SAMPLES" ]]; then
  for c in "$SCRATCH/asm_samples.txt" "$SCRATCH/samples.txt"; do
    [[ -s "$c" ]] && { SAMPLES="$c"; break; }
  done
fi
if [[ -z "$SAMPLES" || ! -s "$SAMPLES" ]]; then
  echo "ERROR: no encuentro la lista de muestras. Pasa SAMPLES=/ruta/lista.txt" >&2
  exit 1
fi
mapfile -t SAMPLE_LIST < <(grep -v '^[[:space:]]*$' "$SAMPLES" | sort -u)
N=${#SAMPLE_LIST[@]}

# ---------------------------------------------------------------------------
# EXTRACTORES
# ---------------------------------------------------------------------------

# Pares desde los .zip de FastQC de un directorio.
# Suma los lanes: el patrón "${s}_R1*" NO matchea "M11_S3_L001_R1_001_fastqc.zip",
# que es exactamente el bug que dejaba raw_pairs en NA. Se usa "${s}_*R1*".
pairs_from_fastqc() {
  local dir="$1" s="$2" total=0 n=0 z v
  [[ -n "$dir" && -d "$dir" ]] || return
  while IFS= read -r z; do
    v=$(unzip -p "$z" '*/fastqc_data.txt' 2>/dev/null \
        | awk -F'\t' '/^Total Sequences/ {print $2; exit}')
    [[ "$v" =~ ^[0-9]+$ ]] || continue
    total=$((total + v)); n=$((n + 1))
  done < <(find "$dir" -maxdepth 1 -name "${s}_*R1*_fastqc.zip" 2>/dev/null | sort)
  (( n > 0 )) && echo "$total"
}

# Pares crudos: primero la línea base, si no los .zip de FastQC.
pairs_raw() {
  local s="$1" v=""
  if [[ -s "$BASELINE" ]]; then
    v=$(awk -F'\t' -v s="$s" '$1 == s {print $2; exit}' "$BASELINE")
  fi
  [[ "$v" =~ ^[0-9]+$ ]] && { echo "$v"; return; }
  pairs_from_fastqc "$FQ_RAW_DIR" "$s"
}

# Pares tras fastp (recorte de polyG).
pairs_fastp() { pairs_from_fastqc "$FQ_NOPOLYG_DIR" "$1"; }

# Pares tras Trimmomatic. Fuente autoritativa: "Both Surviving" del log.
# Trimmomatic escribe a STDERR, así que con -o/-e separados cae en el .err:
# se miran los dos. Se coge el intento MÁS RECIENTE por mtime.
pairs_trim() {
  local s="$1" f e v
  f=$(grep -l "Trimming ${s} " "$LOGS_DIR"/trim.*.log 2>/dev/null \
      | xargs -r ls -t 2>/dev/null | head -1)
  if [[ -n "$f" ]]; then
    e="${f%.log}.err"
    v=$({ cat "$f"; [[ -f "$e" ]] && cat "$e"; } 2>/dev/null \
        | sed -n 's/.*Both Surviving: \([0-9]*\).*/\1/p' | head -1)
    [[ "$v" =~ ^[0-9]+$ ]] && { echo "$v"; return; }
  fi
  pairs_from_fastqc "$FQ_TRIMMED_DIR" "$s"   # respaldo
}

# --- bowtie2 ---------------------------------------------------------------
bt_summary() {
  local s="$1" f
  for f in "$BT_STATS_DIR/${s}_bowtie2_summary.txt" "$BT_STATS_DIR/${s}.txt"; do
    [[ -s "$f" ]] && { echo "$f"; return; }
  done
}
bt_total()  { sed -n 's/^\([0-9]*\) reads; of these:.*/\1/p' "$1" 2>/dev/null | head -1; }
bt_rate()   { sed -n 's/^\([0-9.]*\)% overall alignment rate.*/\1/p' "$1" 2>/dev/null | head -1; }

# Pares NO mapeados de bowtie2 = lo que --un-conc-gz escribió de verdad.
# Ojo: son parejas NO CONCORDANTES, algo más que "ambos extremos sin mapear"
# (~1.4 % de diferencia medida). Es la definición que alimenta al ensamblaje,
# así que es la que se grafica; el overall alignment rate va aparte, en
# read_recovery_metrics.tsv, porque es POR LECTURA y no por pareja.
bt_unmapped() {
  local f="$UNMAPPED_COUNTS_DIR/${1}_unmapped_pairs.tsv"
  [[ -s "$f" ]] || return
  awk -F'\t' '$2 ~ /^[0-9]+$/ {print $2; exit}' "$f"
}

# --- STAR ------------------------------------------------------------------
# STAR nombra su salida según --outFileNamePrefix, y hay dos formas habituales:
#   a) prefijo en el fichero      aligned/M1_S1_Log.final.out
#   b) un directorio por muestra  aligned/M1_S1/Log.final.out   <- basename sin muestra
# Se prueban las dos, y luego un barrido por ruta. El ancla [/_.] evita que
# M1_S1 matchee un hipotético M1_S10, que es la confusión de siempre.
star_log() {
  local s="$1" f
  f=$(find "$STAR_DIR" -maxdepth 3 \
        \( -name "${s}_Log.final.out" -o -name "${s}Log.final.out" \
           -o -name "${s}_*Log.final.out" \) 2>/dev/null \
      | xargs -r ls -t 2>/dev/null | head -1)
  [[ -n "$f" ]] && { echo "$f"; return; }

  f=$(find "$STAR_DIR" -maxdepth 3 -path "*/${s}/*" -name "*Log.final.out" 2>/dev/null \
      | xargs -r ls -t 2>/dev/null | head -1)
  [[ -n "$f" ]] && { echo "$f"; return; }

  find "$STAR_DIR" -maxdepth 4 -name "*Log.final.out" 2>/dev/null \
    | grep -E "/${s}[/_.]" | xargs -r ls -t 2>/dev/null | head -1
}
star_field() {   # $1=log  $2=etiqueta exacta antes del |
  awk -F'|' -v k="$2" '
    { gsub(/^[ \t]+|[ \t]+$/, "", $1)
      if ($1 == k) { gsub(/[ \t%]/, "", $2); print $2; exit } }' "$1" 2>/dev/null
}

# ---------------------------------------------------------------------------
# MODO --check
# ---------------------------------------------------------------------------
if (( CHECK )); then
  printf '\n== Fuentes localizadas ==\n'
  printf '  %-22s %s\n' "muestras"        "$SAMPLES ($N)"
  printf '  %-22s %s\n' "línea base"      "${BASELINE}$( [[ -s $BASELINE ]] && echo '' || echo '   [no existe]')"
  printf '  %-22s %s\n' "FastQC crudos"   "${FQ_RAW_DIR:-[NO ENCONTRADO]}"
  printf '  %-22s %s\n' "FastQC nopolyg"  "${FQ_NOPOLYG_DIR:-[NO ENCONTRADO]}"
  printf '  %-22s %s\n' "FastQC trimmed"  "${FQ_TRIMMED_DIR:-[NO ENCONTRADO]}"
  printf '  %-22s %s\n' "logs Trimmomatic" "$LOGS_DIR/trim.*.log|.err"
  printf '  %-22s %s\n' "summary bowtie2" "${BT_STATS_DIR:-[NO ENCONTRADO]}"
  printf '  %-22s %s\n' "no-mapeados bt2" "${UNMAPPED_COUNTS_DIR:-[NO ENCONTRADO]}"
  printf '  %-22s %s\n' "STAR"            "${STAR_DIR:-[NO ENCONTRADO]}"

  printf '\n== Cobertura por fuente (muestras con dato / %d) ==\n' "$N"
  declare -A hit=()
  for s in "${SAMPLE_LIST[@]}"; do
    [[ -n "$(pairs_raw   "$s")" ]] && hit[raw]=$(( ${hit[raw]:-0} + 1 ))
    [[ -n "$(pairs_fastp "$s")" ]] && hit[fastp]=$(( ${hit[fastp]:-0} + 1 ))
    [[ -n "$(pairs_trim  "$s")" ]] && hit[trim]=$(( ${hit[trim]:-0} + 1 ))
    [[ -n "$(bt_summary  "$s")" ]] && hit[bt]=$(( ${hit[bt]:-0} + 1 ))
    [[ -n "$(bt_unmapped "$s")" ]] && hit[btun]=$(( ${hit[btun]:-0} + 1 ))
    [[ -n "$(star_log    "$s")" ]] && hit[star]=$(( ${hit[star]:-0} + 1 ))
  done
  printf '  %-28s %3d\n' "crudas"                "${hit[raw]:-0}"
  printf '  %-28s %3d\n' "tras fastp"            "${hit[fastp]:-0}"
  printf '  %-28s %3d\n' "tras Trimmomatic"      "${hit[trim]:-0}"
  printf '  %-28s %3d\n' "summary de bowtie2"    "${hit[bt]:-0}"
  printf '  %-28s %3d\n' "no-mapeados de bowtie2" "${hit[btun]:-0}"
  printf '  %-28s %3d\n' "Log.final.out de STAR" "${hit[star]:-0}"
  # Si STAR sale a 0, no dejarlo en "no encontrado": enseñar qué hay de verdad.
  if [[ "${hit[star]:-0}" -eq 0 ]]; then
    printf '\n-- STAR a 0: esto es lo que hay realmente --\n'
    if [[ -n "${STAR_DIR:-}" && -d "$STAR_DIR" ]]; then
      printf '  contenido de %s (primeras 10 entradas):\n' "$STAR_DIR"
      ls -1 "$STAR_DIR" 2>/dev/null | head -10 | sed 's/^/    /'
      [[ -z "$(ls -A "$STAR_DIR" 2>/dev/null)" ]] && printf '    (vacío)\n'
    else
      printf '  no hay directorio candidato de STAR\n'
    fi
    printf '  Log.final.out en todo $SCRATCH:\n'
    find "$SCRATCH" -maxdepth 4 -name "*Log.final.out" 2>/dev/null | head -5 | sed 's/^/    /'
    if ! find "$SCRATCH" -maxdepth 4 -name "*Log.final.out" 2>/dev/null | grep -q .; then
      printf '    ninguno. Los Log.final.out de STAR ya no están en scratch.\n'
      printf '    Sigue adelante sin STAR: el script de R saca los mapeadores del\n'
      printf '    propio TSV, así que con solo bowtie2 hace un panel en vez de dos.\n'
    fi
  fi

  printf '\nSi alguna sale a 0, forza su ruta por variable de entorno y repite.\n'
  printf 'Ej:  STAR_DIR=$SCRATCH/otra_carpeta ./collect_recovery.sh --check\n\n'
  exit 0
fi

# ---------------------------------------------------------------------------
# RECOLECCIÓN
# ---------------------------------------------------------------------------
printf 'sample\tmapper\tstep\tstep_order\tpairs\n' > "$LONG"
printf 'sample\traw_pairs\tfastp_pairs\ttrim_pairs\tbt2_total\tbt2_mapped\tbt2_unmapped\tbt2_overall_rate\tstar_input\tstar_mapped\tstar_unmapped\tstar_toomany_pct\n' > "$WIDE"
printf 'sample\tmapper\tpct_survive_fastp\tpct_survive_trim\tpct_survive_total\tpct_mapped_of_trim\tpct_unmapped_of_trim\tpct_unmapped_of_raw\n' > "$METRICS"

na() { [[ -n "${1:-}" ]] && echo "$1" || echo "NA"; }
pct() {  # $1 num  $2 den
  [[ "${1:-}" =~ ^[0-9]+$ && "${2:-}" =~ ^[0-9]+$ && "$2" -gt 0 ]] || { echo "NA"; return; }
  awk -v n="$1" -v d="$2" 'BEGIN{printf "%.4f", 100*n/d}'
}
emit() {  # sample mapper step order pairs
  [[ "${5:-}" =~ ^[0-9]+$ ]] || return
  printf '%s\t%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$4" "$5" >> "$LONG"
}

for s in "${SAMPLE_LIST[@]}"; do
  raw=$(pairs_raw "$s")
  fp=$(pairs_fastp "$s")
  tr=$(pairs_trim  "$s")

  # --- bowtie2
  btf=$(bt_summary "$s"); bt_tot=""; bt_r=""
  [[ -n "$btf" ]] && { bt_tot=$(bt_total "$btf"); bt_r=$(bt_rate "$btf"); }
  bt_un=$(bt_unmapped "$s")
  # base de pares que entró al mapeo: el total del summary, y si no, lo recortado
  bt_base="${bt_tot:-$tr}"
  bt_map=""
  [[ "$bt_base" =~ ^[0-9]+$ && "$bt_un" =~ ^[0-9]+$ ]] && bt_map=$((bt_base - bt_un))

  # --- STAR
  slog=$(star_log "$s"); st_in=""; st_u=""; st_m=""; st_many=""
  if [[ -n "$slog" ]]; then
    st_in=$(star_field "$slog" "Number of input reads")
    uq=$(star_field "$slog" "Uniquely mapped reads number")
    mu=$(star_field "$slog" "Number of reads mapped to multiple loci")
    st_many=$(star_field "$slog" "% of reads mapped to too many loci")
    [[ "$uq" =~ ^[0-9]+$ && "$mu" =~ ^[0-9]+$ ]] && st_m=$((uq + mu))
    [[ "$st_in" =~ ^[0-9]+$ && "$st_m" =~ ^[0-9]+$ ]] && st_u=$((st_in - st_m))
  fi

  # --- formato largo. Las etapas de limpieza son comunes a los dos mapeadores:
  # se marcan "shared" y el script de R las replica en cada panel.
  emit "$s" shared   "Raw"                1 "$raw"
  emit "$s" shared   "After fastp"        2 "$fp"
  emit "$s" shared   "After Trimmomatic"  3 "$tr"
  emit "$s" bowtie2  "Mapped to host"     4 "$bt_map"
  emit "$s" bowtie2  "Unmapped"           5 "$bt_un"
  emit "$s" STAR     "Mapped to host"     4 "$st_m"
  emit "$s" STAR     "Unmapped"           5 "$st_u"

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$s" \
    "$(na "$raw")" "$(na "$fp")" "$(na "$tr")" \
    "$(na "$bt_tot")" "$(na "$bt_map")" "$(na "$bt_un")" "$(na "$bt_r")" \
    "$(na "$st_in")" "$(na "$st_m")" "$(na "$st_u")" "$(na "$st_many")" >> "$WIDE"

  printf '%s\tbowtie2\t%s\t%s\t%s\t%s\t%s\t%s\n' "$s" \
    "$(pct "$fp" "$raw")" "$(pct "$tr" "$fp")" "$(pct "$tr" "$raw")" \
    "$(pct "$bt_map" "$bt_base")" "$(pct "$bt_un" "$bt_base")" "$(pct "$bt_un" "$raw")" >> "$METRICS"
  printf '%s\tSTAR\t%s\t%s\t%s\t%s\t%s\t%s\n' "$s" \
    "$(pct "$fp" "$raw")" "$(pct "$tr" "$fp")" "$(pct "$tr" "$raw")" \
    "$(pct "$st_m" "$st_in")" "$(pct "$st_u" "$st_in")" "$(pct "$st_u" "$raw")" >> "$METRICS"
done

echo "Muestras procesadas: $N"
echo "  $LONG      ($(($(wc -l < "$LONG") - 1)) filas)"
echo "  $WIDE"
echo "  $METRICS"
echo
echo "Cópialos a tu equipo y pásaselos al script de R:"
echo "  scp <usuario>@<login>:$LONG ."
