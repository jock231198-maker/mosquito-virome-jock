#!/bin/bash
# ---------------------------------------------------------------------------
# prep_lote.sh — staging del lote activo. No copia nada: symlinks.
#
#   source env_lote.sh MEXICO_INSP && ./prep_lote.sh
#   source env_lote.sh UGANDA_UVRI && ./prep_lote.sh
#
# QUE HACE
#   Pone los ficheros del lote con el nombre EXACTO que esperan los workers, en
#   el sitio por el que ese lote entra al pipeline, y escribe $SCRATCH/samples.txt
#   (mismo nombre y mismo formato que en MERIDA, para que todo lo demas encaje).
#
#   LOTE_KIND=crudo      -> $SCRATCH/cat_fastq/<s>_R1_001.fastq.gz
#                           entra por polyg_worker.sh
#   LOTE_KIND=recortado  -> $SCRATCH/trimmed/$RESULT_FROM/<s>_R1_001_paired.fastq.gz
#                           entra por bowtie_map_worker_fast.sh
#
#   Tres normalizaciones de nombre, todas necesarias:
#     _L001 fuera           los workers no lo esperan (y aqui no hay lanes que unir)
#     .fq.gz -> .fastq.gz   los workers buscan .fastq.gz literal
#     $LOTE_PREFIX delante  un nombre que empieza por digito (2_S3) se rompe en R
#                           y se coerciona a numero en awk
#
# Opciones:
#   --count        cuenta reads de cada R1 (lento: ~1 min/GB)
#   --forzar-polyg manda un lote 'recortado' a cat_fastq, para pasarlo por fastp
# ---------------------------------------------------------------------------
set -euo pipefail

: "${LOTE:?source env_lote.sh <LOTE> primero}"
: "${SCRATCH:?}" ; : "${RESULTS_DIR:?}" ; : "${RESULT_FROM:?}"
: "${LOTE_DATA:?}" ; : "${LOTE_KIND:?}" ; : "${LOTE_GLOB:?}"
LOTE_PREFIX="${LOTE_PREFIX:-}"

DO_COUNT=0; FORZAR_POLYG=0
for a in "$@"; do
  case "$a" in
    --count)         DO_COUNT=1 ;;
    --forzar-polyg)  FORZAR_POLYG=1 ;;
    *) echo "opcion desconocida: $a"; exit 2 ;;
  esac
done

kind="$LOTE_KIND"
nota=""
if (( FORZAR_POLYG )) && [[ "$kind" != "crudo" ]]; then
  kind="crudo"; nota="  -> FORZADO a crudo: pasara por fastp"
fi

if [[ "$kind" == "crudo" ]]; then
  DEST="$SCRATCH/cat_fastq"; SUF1="_R1_001.fastq.gz";        SUF2="_R2_001.fastq.gz"
  ENTRADA="polyg_worker.sh"
else
  DEST="$SCRATCH/trimmed/$RESULT_FROM"; SUF1="_R1_001_paired.fastq.gz"; SUF2="_R2_001_paired.fastq.gz"
  ENTRADA="bowtie_map_worker_fast.sh"
fi

QC_DIR="$RESULTS_DIR/qc_control"
MAP="$QC_DIR/staging_map.tsv"
mkdir -p "$DEST" "$QC_DIR" "$LOGS_DIR"
printf 'lote\tsample\tnombre_original\tread\tdestino\torigen\n' > "$MAP"

# Cabecera del primer read sin que el SIGPIPE de gzip tumbe el script
first_header() { { gzip -cd "$1" || true; } | head -n 1; }

# Solo las maquinas de dos colores generan colas de polyG
instrumento() {
  local id="${1#@}"; id="${id%%:*}"
  case "$id" in
    M0*|M1*|M2*|M3*|M4*|M5*|M6*|M7*|M8*|M9*) echo "$id  MiSeq (4 colores, SIN polyG)" ;;
    NB*|NS*)                                  echo "$id  NextSeq 500/550 (2 colores, CON polyG)" ;;
    VH*|VL*)                                  echo "$id  NextSeq 1000/2000 (2 colores, CON polyG)" ;;
    A0*|A1*)                                  echo "$id  NovaSeq 6000 (2 colores, CON polyG)" ;;
    LH*)                                      echo "$id  NovaSeq X (2 colores, CON polyG)" ;;
    D0*|SN*|K0*|J0*)                          echo "$id  HiSeq (4 colores, SIN polyG)" ;;
    "")                                       echo "(cabecera vacia o fichero ilegible)" ;;
    *)                                        echo "$id  no reconocido - miralo a mano" ;;
  esac
}

echo
echo "=== $LOTE  ($LOTE_KIND)$nota"
echo "    origen  : $LOTE_DATA"
echo "    destino : $DEST"
echo "    entra por: $ENTRADA"
[[ -d "$LOTE_DATA" ]] || { echo "ERROR: no existe $LOTE_DATA"; exit 1; }

shopt -s nullglob
files=("$LOTE_DATA"/$LOTE_GLOB)
shopt -u nullglob
(( ${#files[@]} )) || { echo "ERROR: 0 ficheros '$LOTE_GLOB' en $LOTE_DATA"; exit 1; }

: > "$SCRATCH/samples.txt"
# La lista de ficheros la escribe el staging, no un `find` posterior.
# MOTIVO: lo que hay en $DEST son SYMLINKS, y `find -type f` NO los encuentra
# (son -type l). Un `find -type f` sobre este directorio devuelve 0 y parece
# que el staging fallo cuando esta perfecto. Lo mismo vale para cualquier otro
# chequeo que uses aqui: necesita `find -L` o `-type l`.
: > "$SCRATCH/filelist.txt"
echo

for r1 in "${files[@]}"; do
  b=$(basename "$r1")
  orig=$(echo "$b" | cut -d_ -f1-2)
  sample="${LOTE_PREFIX}${orig}"

  # El R2 se deriva del nombre del R1, sea cual sea la convencion del lote
  case "$b" in
    *_val_1.fq.gz)    r2="${r1/_R1_001_val_1.fq.gz/_R2_001_val_2.fq.gz}" ;;
    *_R1_001.fastq.gz) r2="${r1/_R1_001.fastq.gz/_R2_001.fastq.gz}" ;;
    *) echo "ERROR: no se deducir el R2 de $b"; exit 1 ;;
  esac
  [[ -f "$r2" ]] || { echo "ERROR: $orig no tiene R2 ($r2)"; exit 1; }

  for pair in "$r1:$SUF1:R1" "$r2:$SUF2:R2"; do
    src="${pair%%:*}"; rest="${pair#*:}"; suf="${rest%%:*}"; rd="${rest##*:}"
    dst="$DEST/${sample}${suf}"
    [[ -r "$src" ]] || { echo "ERROR: no puedo leer $src"; exit 1; }
    ln -sfn "$src" "$dst"
    [[ -r "$dst" ]] || { echo "ERROR: el symlink $dst no resuelve"; exit 1; }
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$LOTE" "$sample" "$orig" "$rd" "$dst" "$src" >> "$MAP"
  done

  echo "$sample" >> "$SCRATCH/samples.txt"
  printf '%s\n%s\n' "$DEST/${sample}${SUF1}" "$DEST/${sample}${SUF2}" >> "$SCRATCH/filelist.txt"
  if [[ "$sample" == "$orig" ]]; then
    printf '  %-14s %8s  %s\n' "$sample" \
      "$(du -h --apparent-size "$r1" | cut -f1)" "$(instrumento "$(first_header "$r1")")"
  else
    printf '  %-14s (era %-7s) %8s  %s\n' "$sample" "$orig" \
      "$(du -h --apparent-size "$r1" | cut -f1)" "$(instrumento "$(first_header "$r1")")"
  fi
  (( DO_COUNT )) && echo "                 reads R1: $(( $(gzip -cd "$r1" | wc -l) / 4 ))"
done

sort -u -o "$SCRATCH/samples.txt" "$SCRATCH/samples.txt"

# --- Emparejamiento R1/R2: el chequeo de la §5 de check_inputs.sh -----------
# Un cruce de parejas conserva el numero de reads, asi que contar no lo detecta.
# Hay que comparar el ID del primer read. Cuesta milisegundos: gzip solo abre
# el primer bloque.
echo
echo "=== Emparejamiento R1/R2 (ID del primer read)"
fallos=0
while read -r s; do
  id1=$(first_header "$DEST/${s}${SUF1}"); id1="${id1%% *}"
  id2=$(first_header "$DEST/${s}${SUF2}"); id2="${id2%% *}"
  if [[ -n "$id1" && "$id1" == "$id2" ]]; then
    printf '  OK    %-14s %s\n' "$s" "$id1"
  else
    printf '  FALLO %-14s R1=%s  R2=%s\n' "$s" "$id1" "$id2"; fallos=$((fallos+1))
  fi
done < "$SCRATCH/samples.txt"

echo
sort -o "$SCRATCH/filelist.txt" "$SCRATCH/filelist.txt"
echo "  samples.txt : $(wc -l < "$SCRATCH/samples.txt") muestras -> $SCRATCH/samples.txt"
echo "  filelist.txt: $(wc -l < "$SCRATCH/filelist.txt") ficheros -> $SCRATCH/filelist.txt  (para FastQC)"
echo "  trazabilidad: $MAP"
if (( fallos )); then
  echo
  echo "!! $fallos parejas con IDs distintos. NO lances nada hasta resolverlo."
  exit 1
fi
echo
echo "Staging de $LOTE listo."