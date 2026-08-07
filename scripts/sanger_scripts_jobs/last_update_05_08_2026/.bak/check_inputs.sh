#!/bin/bash
# ---------------------------------------------------------------------------
# check_inputs.sh — valida los FASTQ crudos ANTES de gastar CPU en ellos.
#
# Comprueba:
#   1. Que el .gz no esta truncado          (gzip -t)
#   2. Que cada R1 tiene su R2              (emparejamiento)
#   3. Que el numero de lineas es multiplo de 4  (--deep)
#   4. Que R1 y R2 tienen el MISMO numero de reads (--deep)
#   5. Ficheros sospechosamente pequenos
#   6. Genera un manifiesto md5 para detectar corrupcion futura
#
# Uso:
#   ./check_inputs.sh                # rapido: gzip -t + emparejamiento
#   ./check_inputs.sh --deep         # ademas cuenta reads (lento, usa bsub)
#   ./check_inputs.sh --md5          # ademas genera/verifica md5
#   ./check_inputs.sh /otra/carpeta  # otra carpeta de entrada
#
# En la farm, --deep sobre muchas muestras conviene enviarlo:
#   bsub -o "$LOGS_DIR/check.%J.log" -q normal -n 4 -M 4000 \
#        -R "select[mem>4000] rusage[mem=4000]" "./check_inputs.sh --deep"
# ---------------------------------------------------------------------------
set -uo pipefail   # sin -e: queremos seguir aunque un fichero falle

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

DEEP=0; DOMD5=0; INDIR=""
for a in "$@"; do
  case "$a" in
    --deep) DEEP=1 ;;
    --md5)  DOMD5=1 ;;
    -*) echo "Opcion desconocida: $a"; exit 1 ;;
    *)  INDIR="$a" ;;
  esac
done
INDIR="${INDIR:-$DATA_DIR}"

REPORT_DIR="$RESULTS_DIR/qc_control"
mkdir -p "$REPORT_DIR"
REPORT="$REPORT_DIR/check_inputs_$(date +%Y%m%d_%H%M).txt"

# Umbral: por debajo de esto casi seguro es un fichero fallido (10 MB)
MIN_BYTES="${MIN_BYTES:-10000000}"

ERRORS=0
WARNS=0
err()  { echo "  [ERROR] $*"; ERRORS=$((ERRORS+1)); }
warn() { echo "  [warn ] $*"; WARNS=$((WARNS+1)); }
ok()   { echo "  [ok   ] $*"; }

{
echo "==========================================================="
echo "check_inputs.sh   $(date)"
echo "Carpeta: $INDIR"
echo "==========================================================="
echo

mapfile -t FILES < <(find "$INDIR" -type f -name "*.fastq.gz" | sort)
echo "Ficheros .fastq.gz encontrados: ${#FILES[@]}"
if [[ ${#FILES[@]} -eq 0 ]]; then
  echo "[ERROR] no hay nada que comprobar"; exit 1
fi
echo

# --- 1. Integridad del gzip -------------------------------------------------
echo "== 1. Integridad del gzip"
for f in "${FILES[@]}"; do
  if gzip -t "$f" 2>/dev/null; then
    :
  else
    err "gzip corrupto o truncado: $f"
  fi
done
[[ $ERRORS -eq 0 ]] && ok "los ${#FILES[@]} ficheros descomprimen bien"
echo

# --- 2. Tamanos -------------------------------------------------------------
echo "== 2. Tamanos de fichero"
for f in "${FILES[@]}"; do
  sz=$(stat -c %s "$f")
  if [[ $sz -lt $MIN_BYTES ]]; then
    warn "$(basename "$f"): solo $(numfmt --to=iec-i --suffix=B "$sz")"
  fi
done
echo "  total: $(du -sh "$INDIR" 2>/dev/null | cut -f1)"
echo

# --- 3. Emparejamiento R1 / R2 ---------------------------------------------
echo "== 3. Emparejamiento R1/R2"
unpaired=0
for f in "${FILES[@]}"; do
  b=$(basename "$f")
  case "$b" in
    *_R1_*) mate="${f/_R1_/_R2_}" ;;
    *_R2_*) continue ;;               # ya se comprueba desde el R1
    *) warn "nombre sin _R1_/_R2_: $b"; continue ;;
  esac
  if [[ ! -f "$mate" ]]; then
    err "falta el R2 de: $b"
    unpaired=$((unpaired+1))
  fi
done
[[ $unpaired -eq 0 ]] && ok "todos los R1 tienen su R2"
echo

# --- 4. Muestras detectadas -------------------------------------------------
echo "== 4. Muestras detectadas"
mapfile -t SAMPLES < <(printf '%s\n' "${FILES[@]}" | xargs -n1 basename \
                       | grep '_R1_' | cut -d_ -f1-2 | sort -u)
echo "  ${#SAMPLES[@]} muestras:"
printf '    %s\n' "${SAMPLES[@]}"
echo

# --- 5. Conteo de reads (--deep) -------------------------------------------
if [[ $DEEP -eq 1 ]]; then
  echo "== 5. Conteo de reads y simetria R1/R2  (lento)"
  printf '%s\t%s\t%s\t%s\n' "sample" "reads_R1" "reads_R2" "estado" \
    > "$REPORT_DIR/raw_read_counts.tsv"
  for s in "${SAMPLES[@]}"; do
    r1=$(find "$INDIR" -name "${s}_R1_*.fastq.gz" | head -1)
    r2="${r1/_R1_/_R2_}"
    [[ -f "$r1" && -f "$r2" ]] || { err "$s: falta R1 o R2"; continue; }

    l1=$(zcat "$r1" | wc -l)
    l2=$(zcat "$r2" | wc -l)
    n1=$((l1 / 4)); n2=$((l2 / 4))

    estado="OK"
    if (( l1 % 4 != 0 )); then err "$s R1: lineas no multiplo de 4 ($l1) -> TRUNCADO"; estado="TRUNCADO_R1"; fi
    if (( l2 % 4 != 0 )); then err "$s R2: lineas no multiplo de 4 ($l2) -> TRUNCADO"; estado="TRUNCADO_R2"; fi
    if (( n1 != n2 ));      then err "$s: R1=$n1 != R2=$n2 -> DESEMPAREJADO"; estado="DESEMPAREJADO"; fi
    [[ "$estado" == "OK" ]] && ok "$s: $n1 pares"

    printf '%s\t%s\t%s\t%s\n' "$s" "$n1" "$n2" "$estado" \
      >> "$REPORT_DIR/raw_read_counts.tsv"
  done
  echo "  tabla: $REPORT_DIR/raw_read_counts.tsv"
else
  echo "== 5. Conteo de reads: OMITIDO (usa --deep)"
fi
echo

# --- 6. Manifiesto md5 ------------------------------------------------------
if [[ $DOMD5 -eq 1 ]]; then
  echo "== 6. Manifiesto md5"
  MD5="$REPORT_DIR/raw_md5.txt"
  if [[ -f "$MD5" ]]; then
    echo "  existe: verificando contra $MD5"
    bad=$( cd "$INDIR" && md5sum -c "$MD5" 2>&1 | grep -v ': OK$' )
    if [[ -n "$bad" ]]; then
      err "ficheros que NO coinciden con el manifiesto:"
      printf '          %s\n' "$bad"
    else
      ok "todos los md5 coinciden"
    fi
  else
    echo "  generando (tarda: lee todos los ficheros)..."
    (cd "$INDIR" && find . -name "*.fastq.gz" | sort | xargs md5sum) > "$MD5"
    ok "creado: $MD5"
  fi
else
  echo "== 6. Manifiesto md5: OMITIDO (usa --md5)"
fi
echo

echo "==========================================================="
echo "RESUMEN: $ERRORS errores, $WARNS avisos"
if [[ $ERRORS -gt 0 ]]; then
  echo "NO SIGAS hasta resolverlos: un fastq truncado corrompe todo aguas abajo."
else
  echo "Entradas validadas. Puedes lanzar FastQC."
fi
echo "==========================================================="

# Sale con 1 si hubo errores. OJO: este bloque corre en un subshell por el
# pipe a tee, asi que el codigo se recupera fuera con PIPESTATUS.
exit $(( ERRORS > 0 ? 1 : 0 ))

} 2>&1 | tee "$REPORT"

rc=${PIPESTATUS[0]}
echo
echo "Informe guardado en: $REPORT"
exit "$rc"
