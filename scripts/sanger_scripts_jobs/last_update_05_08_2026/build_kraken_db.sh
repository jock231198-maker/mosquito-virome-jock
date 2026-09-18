#!/bin/bash
# build_kraken_db.sh — descarga una base Kraken2 pre-construida (incluye Bracken)
#
# Uso:
#   K2_URL="https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_YYYYMMDD.tar.gz" \
#   bsub ... build_kraken_db.sh
#
# La URL NO se fija aquí a proposito: las versiones cambian y una URL caducada
# deja un fichero de 0 bytes sin avisar. Sacala de:
#   curl -s https://benlangmead.github.io/aws-indexes/k2 \
#     | grep -oE 'https://genome-idx[^"]*k2_pluspf_2[^"]*\.tar\.gz' | sort -u | tail -3
#
# El paquete trae la base de Kraken2 + las de Bracken (100, 150 y 200 meros),
# asi que no hay que construir nada.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

: "${K2_URL:?define K2_URL con la URL del .tar.gz (ver cabecera)}"
DB="${KRAKEN_DB:-$REFS_DIR/kraken2/pluspf}"
mkdir -p "$DB"
cd "$DB"

tarball="$(basename "$K2_URL")"

# ---- comprobacion de espacio -------------------------------------------------
echo "== espacio en el grupo =="
lfs quota -h -g team222 /lustre/scratch126 2>/dev/null | tail -2 || true
echo

# ---- si ya esta, no rehacer --------------------------------------------------
if [[ -s "$DB/hash.k2d" && -s "$DB/taxo.k2d" && -s "$DB/opts.k2d" ]]; then
  echo "== la base ya existe en $DB =="
  ls -lh "$DB"/*.k2d
  exit 0
fi

# ---- descarga ----------------------------------------------------------------
echo "== descargando $tarball =="
wget -c -O "$tarball" "$K2_URL"

# un .tar.gz de 0 bytes es el sintoma de una URL mala: fallar aqui, no despues
if [[ ! -s "$tarball" ]]; then
  echo "ERROR: $tarball quedo vacio. La URL es incorrecta o caduco."
  rm -f "$tarball"; exit 1
fi
echo "   descargado: $(du -h "$tarball" | cut -f1)"

echo "== extrayendo =="
tar -xzf "$tarball"

# ---- verificacion ------------------------------------------------------------
echo
echo "== ficheros de la base =="
ls -lh "$DB"
echo
for f in hash.k2d opts.k2d taxo.k2d; do
  [[ -s "$DB/$f" ]] && echo "  OK    $f  ($(du -h "$DB/$f" | cut -f1))" \
                    || { echo "  FALTA $f"; exit 1; }
done
echo
echo "== bases de Bracken disponibles =="
ls -1 "$DB"/database*mers.kmer_distrib 2>/dev/null || echo "  ninguna (Bracken no se podra usar)"

echo
echo "Base lista en: $DB"
echo "Memoria que hara falta al clasificar: ~$(du -sh "$DB/hash.k2d" | cut -f1) + margen"
echo
echo "Puedes borrar el tarball: rm $DB/$tarball"
