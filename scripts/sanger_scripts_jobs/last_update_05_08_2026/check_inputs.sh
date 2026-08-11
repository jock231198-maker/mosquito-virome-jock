#!/bin/bash
# ---------------------------------------------------------------------------
# check_inputs.sh — valida FASTQ ANTES de gastar CPU en ellos.
#
# Funciona sobre los DOS layouts del proyecto:
#   crudos    $DATA_DIR            M11_S3_L001_R1_001.fastq.gz   (con lane)
#   cat_fastq $SCRATCH/cat_fastq   M11_S3_R1_001.fastq.gz        (sin lane)
#
# Comprueba:
#   1. Integridad del gzip                             (gzip -t, en paralelo)
#   2. Tamanos: umbral absoluto + desviacion vs mediana
#   3. Emparejamiento R1/R2 en disco
#   4. Muestras detectadas; con lanes, que ninguna cojee un lane
#   5. CRUCE R1/R2 por ID del primer read              <-- barato y critico
#   6. --deep   conteo de reads, multiplo de 4, simetria, ID del ultimo read
#   7. --vs-raw el concatenado cuadra byte a byte con los crudos
#   8. --md5    manifiesto para detectar corrupcion futura
#
# Uso:
#   ./check_inputs.sh                                   # crudos, modo rapido
#   ./check_inputs.sh "$SCRATCH/cat_fastq"              # concatenados
#   ./check_inputs.sh --vs-raw "$SCRATCH/cat_fastq"     # + cuadre con crudos
#   ./check_inputs.sh --deep "$SCRATCH/cat_fastq"       # + conteo de reads
#   ./check_inputs.sh --md5 --expect 23 "$SCRATCH/cat_fastq"
#
# El modo rapido descomprime solo la primera cabecera de cada fichero: son
# segundos. --deep descomprime entero: eso va SIEMPRE por bsub.
#
#   bsub -J chkdeep -o "$LOGS_DIR/chk.%J.log" -e "$LOGS_DIR/chk.%J.err" \
#        -q normal -n 8 -M 4000 -R "select[mem>4000] rusage[mem=4000] span[hosts=1]" \
#        "$PWD/check_inputs.sh --deep $SCRATCH/cat_fastq"
# ---------------------------------------------------------------------------
set -uo pipefail   # sin -e: queremos seguir aunque un fichero falle

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# --- Modo interno: escanea UN fichero de una sola pasada ---------------------
# Se autoinvoca desde xargs. Va antes de sourcear config.sh: no necesita nada.
# Salida: fichero <TAB> n_lineas <TAB> primera_cabecera <TAB> ultima_cabecera
if [[ "${1:-}" == "--scan-one" ]]; then
  gzip -cd -- "$2" 2>/dev/null | awk -v F="$2" '
    NR==1        {a=$1}
    (NR-1)%4==0  {b=$1}
    END          {printf "%s\t%d\t%s\t%s\n", F, NR, a, b}'
  exit 0
fi

source "$SCRIPT_DIR/config.sh"

# --- Argumentos -------------------------------------------------------------
DEEP=0; DOMD5=0; VSRAW=0; EXPECT=""; INDIR=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --deep)   DEEP=1 ;;
    --md5)    DOMD5=1 ;;
    --vs-raw) VSRAW=1 ;;
    --expect) EXPECT="${2:?--expect necesita un numero}"; shift ;;
    -h|--help) sed -n '2,40p' "$0"; exit 0 ;;
    -*)  echo "Opcion desconocida: $1"; exit 1 ;;
    *)   INDIR="$1" ;;
  esac
  shift
done
INDIR="${INDIR:-$DATA_DIR}"
INDIR="${INDIR%/}"

NPROC="${LSB_DJOB_NUMPROC:-4}"
MIN_BYTES="${MIN_BYTES:-10000000}"        # 10 MB: por debajo casi seguro es fallo
SAMPLE_FIELDS="${SAMPLE_FIELDS:-2}"       # M11_S3 -> 2 campos separados por '_'

TAG=$(basename "$INDIR")
REPORT_DIR="$RESULTS_DIR/qc_control"
mkdir -p "$REPORT_DIR"
REPORT="$REPORT_DIR/check_inputs_${TAG}_$(date +%Y%m%d_%H%M).txt"

TMP=$(mktemp -d "${TMPDIR:-/tmp}/chkin.XXXXXX")
trap 'rm -rf "$TMP"' EXIT

ERRORS=0; WARNS=0
err()  { echo "  [ERROR] $*"; ERRORS=$((ERRORS+1)); }
warn() { echo "  [warn ] $*"; WARNS=$((WARNS+1)); }
ok()   { echo "  [ok   ] $*"; }

hum() { numfmt --to=iec-i --suffix=B "$1" 2>/dev/null || echo "${1}B"; }

# Quita '@', se queda con el primer campo y descarta el sufijo /1 o /2.
# Casava 1.8+:  @A00.. :1234:5678 1:N:0:ATCACG   -> A00..:1234:5678
# Formato viejo: @READID/1                        -> READID
norm_id() {
  local s="${1#@}"
  s="${s%% *}"; s="${s%/1}"; s="${s%/2}"
  printf '%s' "$s"
}

# Mate de un R1: se sustituye SOLO en el basename, nunca en la ruta.
mate_of() {
  local f="$1" d b
  d=$(dirname "$f"); b=$(basename "$f")
  printf '%s/%s' "$d" "${b/_R1_/_R2_}"
}

sample_of() { basename "$1" | cut -d_ -f1-"$SAMPLE_FIELDS"; }

# Agrega el escaneo de TODOS los ficheros de una muestra+read (varios si hay
# lanes) ordenados por nombre. Devuelve: n_lineas <TAB> primer_id <TAB>
# ultimo_id <TAB> n_ficheros.  Con lanes hay que SUMAR, no coger el primero.
agg_scan() {   # agg_scan <sample> <R1|R2> <fichero_scan>
  sort "$3" | awk -F'\t' -v s="/$1_" -v r="_$2_" '
    index($1, s) && index($1, r) { n += $2; if (c == 0) f = $3; l = $4; c++ }
    END { printf "%d\t%s\t%s\t%d\n", n, f, l, c }'
}

{
echo "==========================================================="
echo "check_inputs.sh   $(date)"
echo "Carpeta : $INDIR"
echo "Modo    : rapido${DEEP:+}$( ((DEEP)) && echo " +deep")$( ((VSRAW)) && echo " +vs-raw")$( ((DOMD5)) && echo " +md5")   (${NPROC} hilos)"
echo "==========================================================="
echo

# --- 0. Inventario y layout -------------------------------------------------
echo "== 0. Inventario"
mapfile -t FILES < <(find "$INDIR" -maxdepth 1 -type f -name "*.fastq.gz" | sort)
if [[ ${#FILES[@]} -eq 0 ]]; then
  echo "  [ERROR] no hay .fastq.gz en $INDIR"; exit 1
fi

mapfile -t R1S < <(printf '%s\n' "${FILES[@]}" | grep '_R1_' || true)
mapfile -t R2S < <(printf '%s\n' "${FILES[@]}" | grep '_R2_' || true)

LANES=0
printf '%s\n' "${FILES[@]}" | grep -q '_L00[0-9]_' && LANES=1

echo "  ficheros .fastq.gz : ${#FILES[@]}   (R1: ${#R1S[@]}, R2: ${#R2S[@]})"
echo "  layout             : $( ((LANES)) && echo "CON lanes (crudos)" || echo "SIN lanes (concatenado)")"
echo "  tamano total       : $(du -sh "$INDIR" 2>/dev/null | cut -f1)"
huerfanos=$(( ${#FILES[@]} - ${#R1S[@]} - ${#R2S[@]} ))
(( huerfanos > 0 )) && warn "$huerfanos ficheros sin _R1_/_R2_ en el nombre"
echo

# --- 1. Integridad del gzip -------------------------------------------------
echo "== 1. Integridad del gzip  (gzip -t, $NPROC en paralelo)"
printf '%s\n' "${FILES[@]}" \
  | xargs -P "$NPROC" -I{} bash -c 'gzip -t "$1" 2>/dev/null || echo "$1"' _ {} \
  > "$TMP/badgz" 2>/dev/null
if [[ -s "$TMP/badgz" ]]; then
  while IFS= read -r f; do err "gzip corrupto o truncado: $f"; done < "$TMP/badgz"
else
  ok "los ${#FILES[@]} ficheros descomprimen bien"
fi
echo

# --- 2. Tamanos -------------------------------------------------------------
echo "== 2. Tamanos"
stat -c '%s %n' "${FILES[@]}" > "$TMP/sizes"
mapfile -t SORTED < <(cut -d' ' -f1 "$TMP/sizes" | sort -n)
MED=${SORTED[$(( ${#SORTED[@]} / 2 ))]}
MAXI=${SORTED[$(( ${#SORTED[@]} - 1 ))]}
echo "  mediana: $(hum "$MED")   min: $(hum "${SORTED[0]}")   max: $(hum "$MAXI")"
while read -r sz name; do
  if (( sz < MIN_BYTES )); then
    err "$(basename "$name"): solo $(hum "$sz")  (umbral $(hum "$MIN_BYTES"))"
  elif (( sz * 2 < MED )); then
    warn "$(basename "$name"): $(hum "$sz"), menos de la mitad de la mediana"
  fi
done < "$TMP/sizes"
echo

# --- 3. Emparejamiento R1 / R2 ---------------------------------------------
echo "== 3. Emparejamiento R1/R2 en disco"
unpaired=0
for f in ${R1S[@]+"${R1S[@]}"}; do
  m=$(mate_of "$f")
  [[ -f "$m" ]] || { err "falta el R2 de: $(basename "$f")"; unpaired=$((unpaired+1)); }
done
if (( ${#R1S[@]} != ${#R2S[@]} )); then
  err "hay ${#R1S[@]} R1 y ${#R2S[@]} R2: no cuadran"
elif (( unpaired == 0 )); then
  ok "los ${#R1S[@]} R1 tienen su R2"
fi
echo

# --- 4. Muestras ------------------------------------------------------------
echo "== 4. Muestras detectadas"
mapfile -t SAMPLES < <(printf '%s\n' ${R1S[@]+"${R1S[@]}"} | while IFS= read -r f; do sample_of "$f"; done | sort -u)
echo "  ${#SAMPLES[@]} muestras"
printf '    %s\n' ${SAMPLES[@]+"${SAMPLES[@]}"} | paste -d' ' - - - - 2>/dev/null || printf '    %s\n' "${SAMPLES[@]}"

if [[ -n "$EXPECT" ]]; then
  if (( ${#SAMPLES[@]} == EXPECT )); then
    ok "coincide con lo esperado ($EXPECT)"
  else
    err "esperabas $EXPECT muestras y hay ${#SAMPLES[@]}"
  fi
fi

# Con lanes: que ninguna muestra cojee un lane
if (( LANES )); then
  echo "  -- completitud de lanes:"
  nlanes_ref=""; lanes_mal=0
  for s in "${SAMPLES[@]}"; do
    n1=$(printf '%s\n' "${R1S[@]}" | grep -c "/${s}_L00[0-9]_R1_")
    n2=$(printf '%s\n' "${R2S[@]}" | grep -c "/${s}_L00[0-9]_R2_")
    [[ -z "$nlanes_ref" ]] && nlanes_ref=$n1
    (( n1 != n2 ))         && { err  "$s: $n1 lanes en R1 pero $n2 en R2"; lanes_mal=1; }
    (( n1 != nlanes_ref )) && { warn "$s: $n1 lanes (el resto tiene $nlanes_ref)"; lanes_mal=1; }
  done
  (( lanes_mal == 0 )) && ok "todas con $nlanes_ref lanes en R1 y en R2"
fi
echo

# --- 5. CRUCE R1/R2 ---------------------------------------------------------
# El emparejamiento FASTQ es POSICIONAL: registro i de R1 <-> registro i de R2.
# Si se concatenan los lanes en distinto orden en R1 y en R2, el numero de
# reads sigue siendo identico y NINGUN chequeo de conteo lo detecta: pasa
# Trimmomatic, pasa bowtie2, y produce un mapeo basura en silencio.
# Comparar el ID del PRIMER read cuesta milisegundos y lo caza.
echo "== 5. Cruce R1/R2 (ID del primer read)"
cross=0
for f in ${R1S[@]+"${R1S[@]}"}; do
  m=$(mate_of "$f"); [[ -f "$m" ]] || continue
  h1=$( { gzip -cd -- "$f" 2>/dev/null || true; } | head -1 )
  h2=$( { gzip -cd -- "$m" 2>/dev/null || true; } | head -1 )
  i1=$(norm_id "$h1"); i2=$(norm_id "$h2")
  if [[ -z "$i1" || -z "$i2" ]]; then
    err "$(basename "$f"): no se pudo leer la primera cabecera"
    cross=$((cross+1))
  elif [[ "$i1" != "$i2" ]]; then
    err "$(sample_of "$f"): R1 y R2 CRUZADOS o desalineados"
    echo "           R1: $i1"
    echo "           R2: $i2"
    cross=$((cross+1))
  fi
done
if (( cross == 0 )); then
  ok "primer read alineado en las ${#R1S[@]} parejas"
else
  echo "  NO SIGAS: un R1/R2 cruzado da un mapeo basura sin lanzar ningun error."
  echo "  Rehaz el concatenado ordenando los lanes con 'sort', identico en R1 y R2."
fi
echo

# --- 6. Conteo de reads (--deep) -------------------------------------------
if (( DEEP )); then
  echo "== 6. Conteo de reads, multiplo de 4 y simetria  (una pasada por fichero)"
  printf '%s\n' "${FILES[@]}" \
    | xargs -P "$NPROC" -n1 "$SCRIPT_DIR/$(basename "$0")" --scan-one \
    > "$TMP/scan" 2>/dev/null

  TSV="$REPORT_DIR/read_counts_${TAG}.tsv"
  printf 'sample\tficheros_R1\treads_R1\treads_R2\tprimer_id_ok\tultimo_id_ok\testado\n' > "$TSV"

  for s in "${SAMPLES[@]}"; do
    IFS=$'\t' read -r l1 f1 u1 c1 < <(agg_scan "$s" R1 "$TMP/scan")
    IFS=$'\t' read -r l2 f2 u2 c2 < <(agg_scan "$s" R2 "$TMP/scan")

    if (( c1 == 0 || c2 == 0 )); then
      err "$s: no se pudo escanear R1 o R2"; continue
    fi
    (( c1 != c2 )) && err "$s: $c1 ficheros R1 vs $c2 R2"
    n1=$((l1 / 4)); n2=$((l2 / 4))
    estado="OK"
    (( l1 % 4 != 0 )) && { err "$s R1: $l1 lineas, no multiplo de 4 -> TRUNCADO"; estado="TRUNCADO_R1"; }
    (( l2 % 4 != 0 )) && { err "$s R2: $l2 lineas, no multiplo de 4 -> TRUNCADO"; estado="TRUNCADO_R2"; }
    (( n1 != n2 ))    && { err "$s: R1=$n1 pares != R2=$n2 -> DESEMPAREJADO";      estado="DESEMPAREJADO"; }

    pid="si"; uid="si"
    [[ "$(norm_id "$f1")" == "$(norm_id "$f2")" ]] || { pid="NO"; estado="CRUZADO"; }
    [[ "$(norm_id "$u1")" == "$(norm_id "$u2")" ]] || {
      uid="NO"
      err "$s: el ULTIMO read no coincide entre R1 y R2 -> CRUZADO"
      estado="CRUZADO"
    }
    [[ "$estado" == "OK" ]] && ok "$s: $n1 pares en $c1 fichero(s), extremos alineados"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$s" "$c1" "$n1" "$n2" "$pid" "$uid" "$estado" >> "$TSV"
  done

  total=$(awk -F'\t' 'NR>1 {s+=$3} END {print s+0}' "$TSV")
  echo "  total: $total pares en ${#SAMPLES[@]} muestras"
  echo "  tabla: $TSV"
else
  echo "== 6. Conteo de reads: OMITIDO (usa --deep, y lanzalo por bsub)"
fi
echo

# --- 7. Cuadre con los crudos (--vs-raw) ------------------------------------
# 'cat a.gz b.gz > c.gz' es concatenacion de miembros gzip: el resultado pesa
# EXACTAMENTE la suma de los dos. Si cuadra al byte, el concatenado esta
# completo y no se recomprimio nada. Cuesta cero: solo lee metadatos.
if (( VSRAW )); then
  echo "== 7. Cuadre byte a byte contra $DATA_DIR"
  if (( LANES )); then
    warn "--vs-raw se usa sobre la carpeta CONCATENADA, no sobre los crudos"
  elif [[ ! -d "$DATA_DIR" ]]; then
    err "no existe DATA_DIR: $DATA_DIR"
  else
    for s in "${SAMPLES[@]}"; do
      for r in R1 R2; do
        cat_f="$INDIR/${s}_${r}_001.fastq.gz"
        [[ -f "$cat_f" ]] || { err "$s: falta $cat_f"; continue; }
        got=$(stat -c %s "$cat_f")
        mapfile -t parts < <(find "$DATA_DIR" -maxdepth 1 -type f \
                             -name "${s}_L*_${r}_*.fastq.gz" | sort)
        if (( ${#parts[@]} == 0 )); then
          warn "$s $r: sin crudos correspondientes en DATA_DIR"; continue
        fi
        want=$(stat -c %s "${parts[@]}" | awk '{s+=$1} END {print s}')
        if (( got == want )); then
          ok "$s $r: ${#parts[@]} lanes -> $(hum "$got") exacto"
        else
          err "$s $r: concatenado $(hum "$got") vs crudos $(hum "$want")  (dif $((got-want)) B)"
        fi
      done
    done
    nraw=$(find "$DATA_DIR" -maxdepth 1 -type f -name "*.fastq.gz" | wc -l)
    echo "  crudos en DATA_DIR: $nraw ficheros"
    (( nraw == 0 )) && err "DATA_DIR VACIO: los crudos eran la unica copia"
  fi
else
  echo "== 7. Cuadre con crudos: OMITIDO (usa --vs-raw)"
fi
echo

# --- 8. Manifiesto md5 ------------------------------------------------------
if (( DOMD5 )); then
  echo "== 8. Manifiesto md5"
  MD5="$REPORT_DIR/md5_${TAG}.txt"
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
    echo "  generando (lee todos los ficheros, tarda)..."
    ( cd "$INDIR" && printf '%s\n' "${FILES[@]}" | xargs -n1 basename \
        | xargs -P "$NPROC" -n1 md5sum ) | sort -k2 > "$MD5"
    ok "creado: $MD5"
  fi
else
  echo "== 8. Manifiesto md5: OMITIDO (usa --md5)"
fi
echo

# --- Resumen ----------------------------------------------------------------
echo "==========================================================="
echo "RESUMEN: $ERRORS errores, $WARNS avisos"
if (( ERRORS > 0 )); then
  echo "NO SIGAS hasta resolverlos: un fastq truncado o cruzado corrompe"
  echo "todo aguas abajo y no vuelve a dar la cara hasta el arbol filogenetico."
else
  echo "Entradas validadas."
fi
echo "==========================================================="

exit $(( ERRORS > 0 ? 1 : 0 ))

} 2>&1 | tee "$REPORT"

rc=${PIPESTATUS[0]}
echo
echo "Informe guardado en: $REPORT"
exit "$rc"