#!/bin/bash
# congelar_hits.sh — congela las tablas de hits virales por contig en
# results/analisis_R/datos_hits_vN, con MANIFEST de md5 (mismo contrato que
# congelar_datos.sh: una version congelada no se toca; si cambia, otra version).
#
# Se corre EN LA FARM, desde el clon del repo:
#   bash results/analisis_R/congelar_hits.sh          # datos_hits_v1
#   bash results/analisis_R/congelar_hits.sh v2     # v2+ incluye blastn.tsv y annot_por_votu.tsv
# Despues: git add/commit/push, y en el Mac: git pull && Rscript hits_tablas.R
#
# Origen: $RESULTS_DIR/union_annot/ (generado a partir de
# votu_annot/all_final_viraldb20260925/annot_por_votu.tsv + geNomad all_final).

set -euo pipefail
V="${1:-v1}"

AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAIZ="$(git -C "$AQUI" rev-parse --show-toplevel)"
SCRIPTS_DIR="${SCRIPTS_DIR:-$RAIZ/scripts/sanger_scripts_jobs/last_update_05_08_2026}"
source "$SCRIPTS_DIR/config.sh" >/dev/null 2>&1
[[ "${RESULT_FROM:-}" == "MERIDA" ]] || { echo "ERROR: el lote activo es '${RESULT_FROM:-?}', no MERIDA"; exit 1; }

O="$RESULTS_DIR/union_annot"
A="$RESULTS_DIR/votu_annot/all_final_viraldb20260925"
GS="$RESULTS_DIR/genomad_all/all_final/all_final_summary/all_final_virus_summary.tsv"
DEST="$AQUI/datos_hits_$V"
[[ -e "$DEST" ]] && { echo "ERROR: $DEST ya existe y esta congelado. Usa otra version."; exit 1; }

# --- comprobaciones: que sea la version de las tablas que esperamos -------------
F=(hits_por_contig.tsv hits_por_muestra.tsv candidatos_por_virus.tsv)
for f in "${F[@]}"; do [[ -s "$O/$f" ]] || { echo "ERROR: falta $O/$f"; exit 1; }; done
# Resultados crudos (desde v2): blastn de los 6,004 contigs y la anotacion por contig
FA=(blastn.tsv annot_por_votu.tsv)
for f in "${FA[@]}"; do [[ -s "$A/$f" ]] || { echo "ERROR: falta $A/$f"; exit 1; }; done
nc=$(head -1 "$A/blastn.tsv" | awk -F'\t' '{print NF}')
[[ $nc -eq 11 ]] || { echo "ERROR: blastn.tsv tiene $nc columnas, se esperaban 11"; exit 1; }

cab=$(head -1 "$O/hits_por_contig.tsv")
for c in muestra contig nombre nivel familia realm grupo en_genomad bn_pident dm_pident; do
  grep -qw "$c" <<<"$cab" || { echo "ERROR: hits_por_contig.tsv sin columna '$c' (version antigua?)"; exit 1; }
done
[[ $(head -1 "$O/hits_por_muestra.tsv" | cut -f1) == muestra ]] || { echo "ERROR: encabezado fuera de sitio en hits_por_muestra.tsv"; exit 1; }
[[ $(head -1 "$O/candidatos_por_virus.tsv" | cut -f1) == nombre ]] || { echo "ERROR: encabezado fuera de sitio en candidatos_por_virus.tsv"; exit 1; }

n_union=$(tail -n +2 "$A/annot_por_votu.tsv" | wc -l)
n_gen=$(tail -n +2 "$GS" | wc -l)
[[ $n_union -eq 6004 ]] || { echo "ERROR: annot_por_votu.tsv tiene $n_union contigs, se esperaban 6004"; exit 1; }
[[ $n_gen -eq 1144 ]]   || { echo "ERROR: geNomad tiene $n_gen virales, se esperaban 1144"; exit 1; }
n_mu=$(tail -n +2 "$O/hits_por_contig.tsv" | cut -f1 | sort -u | wc -l)
[[ $n_mu -eq 22 ]] || echo "AVISO: hits_por_contig.tsv tiene $n_mu muestras (se esperaban 22)"

# --- copia + MANIFEST ------------------------------------------------------------
TMPD="$DEST.tmp.$$"; trap 'rm -rf "$TMPD"' EXIT
mkdir -p "$TMPD"
cp "${F[@]/#/$O/}" "${FA[@]/#/$A/}" "$TMPD/"
{
  echo "# datos_hits_$V — congelado $(date '+%F %T') por ${USER:-$(id -un)} en $(hostname)"
  echo "# origen: $O"
  echo "# anotacion: $A/annot_por_votu.tsv ($n_union contigs)"
  echo "# geNomad: $GS ($n_gen virales)"
  echo "# commit del repo: $(git -C "$RAIZ" rev-parse --short HEAD)"
  echo -e "# md5\tlineas\tfichero"
  for f in "${F[@]}" "${FA[@]}"; do echo -e "$(md5sum "$TMPD/$f" | cut -d' ' -f1)\t$(wc -l < "$TMPD/$f")\t$f"; done
} > "$TMPD/MANIFEST.txt"
mv "$TMPD" "$DEST"   # solo aparece datos_hits_vN si todo salio bien

cat "$DEST/MANIFEST.txt"
echo
echo "Siguiente:"
echo "  git -C \"$RAIZ\" add results/analisis_R/datos_hits_$V && git -C \"$RAIZ\" commit -m 'datos_hits_$V congelados' && git -C \"$RAIZ\" push"
