#!/bin/bash
# ---------------------------------------------------------------------------
# build_sortmerna_db.sh — descarga la base de rRNA de SortMeRNA y CONSTRUYE EL
#                         INDICE UNA SOLA VEZ, para que el array no lo repita 22.
#
# Mismo espiritu que build_bt2_index.sh y build_star_index.sh: idempotente,
# verifica al final, y no es un job array.
#
# POR QUE UN SCRIPT APARTE
#   SortMeRNA indexa la referencia en la primera ejecucion y la reutiliza si
#   encuentra el indice. Si lanzas las 22 muestras de golpe sin indice previo,
#   las 22 intentan construirlo A LA VEZ en el mismo sitio: carrera de escritura,
#   indice corrupto, y 22 veces el mismo trabajo. Se construye antes, una vez.
#
# SALIDA:
#   $REFS_DIR/sortmerna/smr_v4.3_default_db.fasta   (y las otras tres)
#   $REFS_DIR/sortmerna/idx/                        indice compartido, solo lectura
#
# Uso:
#   ./build_sortmerna_db.sh              # base 'default'
#   ./build_sortmerna_db.sh sensitive    # 2x mas lento, +0.008% de exactitud
#
# Envio (no es array, y la descarga necesita red: nodo de login o cola normal):
#   bsub -J smrdb -o "$LOGS_DIR/smrdb.%J.log" -e "$LOGS_DIR/smrdb.%J.err" \
#        -q normal -n 4 -M 8000 \
#        -R "select[mem>8000] rusage[mem=8000] span[hosts=1]" \
#        "$SCRIPTS_DIR/build_sortmerna_db.sh"
# ---------------------------------------------------------------------------
set -euo pipefail
SECONDS=0

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

variant="${1:-default}"
SMR_DIR="${SMR_DIR:-$REFS_DIR/sortmerna}"
DB_URL="${DB_URL:-https://github.com/sortmerna/sortmerna/releases/download/v4.3.4/database.tar.gz}"
THREADS="${LSB_DJOB_NUMPROC:-4}"

case "$variant" in
  fast|default|sensitive) db="$SMR_DIR/smr_v4.3_${variant}_db.fasta" ;;
  *) echo "ERROR: variante '$variant' no valida (fast|default|sensitive)"; exit 1 ;;
esac

mkdir -p "$SMR_DIR"

activate_env "$ENV_SORTMERNA"
command -v sortmerna >/dev/null || { echo "ERROR: sortmerna no esta en el PATH"; exit 1; }
echo "sortmerna: $(sortmerna --version 2>&1 | head -1)"

# --- 1. La base ------------------------------------------------------------
if [[ -s "$db" ]]; then
  echo "[ok] la base ya esta: $db"
else
  echo "== descargando database.tar.gz =="
  tgz="$SMR_DIR/database.tar.gz"
  [[ -s "$tgz" ]] || wget -O "$tgz" "$DB_URL"
  tar -xzf "$tgz" -C "$SMR_DIR"
  # el tarball puede traer las fasta sueltas o dentro de un subdirectorio
  if [[ ! -s "$db" ]]; then
    found=$(find "$SMR_DIR" -name "smr_v4.3_${variant}_db.fasta" | head -1)
    [[ -n "$found" ]] || { echo "ERROR: no aparece smr_v4.3_${variant}_db.fasta tras extraer"; ls -R "$SMR_DIR" | head -40; exit 1; }
    [[ "$found" == "$db" ]] || cp "$found" "$db"
  fi
fi

nseq=$(grep -c '^>' "$db")
echo "[ok] $db  ($nseq secuencias de referencia)"
(( nseq > 1000 )) || { echo "ERROR: la base parece truncada"; exit 1; }

# --- 2. El indice ----------------------------------------------------------
IDX="$SMR_DIR/idx_${variant}"
if [[ -d "$IDX" ]] && (( $(find "$IDX" -type f -size +0 | wc -l) >= 8 )); then
  echo "[ok] el indice ya esta: $IDX  ($(find "$IDX" -type f | wc -l) ficheros)"
else
  echo "== construyendo el indice (una sola vez) =="
  rm -rf "$IDX"; mkdir -p "$IDX"

  # SortMeRNA indexa como efecto secundario de una ejecucion, no hay subcomando.
  # Se le da un puñado de lecturas de mentira: lo unico que interesa es el idx.
  tmp="${TMPDIR:-/tmp}/smrbuild_$$"
  mkdir -p "$tmp/kvdb"
  printf '@r1\nACGTACGTACGTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n' \
    > "$tmp/dummy.fastq"

  sortmerna \
      --ref "$db" \
      --reads "$tmp/dummy.fastq" \
      --workdir "$tmp" \
      --idx-dir "$IDX" \
      --kvdb "$tmp/kvdb" \
      --threads "$THREADS" \
      --num_alignments 1 \
      -v

  rm -rf "$tmp"
fi

# --- 3. Verificacion -------------------------------------------------------
echo
echo "==========================================================="
nfiles=$(find "$IDX" -type f -size +0 | wc -l)
size=$(du -sh "$IDX" | cut -f1)
if (( nfiles >= 8 )); then
  echo "[ok] indice: $IDX   ($nfiles ficheros, $size)"
else
  echo "[FALLO] el indice solo tiene $nfiles ficheros no vacios. Revisa el log."
  ls -la "$IDX"; exit 1
fi
echo "Para el worker:  SMR_DB=$db  SMR_IDX=$IDX"
echo "Done  ($((SECONDS/60))m $((SECONDS%60))s)"
