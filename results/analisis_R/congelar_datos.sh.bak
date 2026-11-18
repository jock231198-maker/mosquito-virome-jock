#!/bin/bash
# congelar_datos.sh — copia de la farm al repo las tablas que usan las figuras,
# en una carpeta versionada (datos_v1, datos_v2, ...) con un MANIFEST de md5.
#
# Se corre EN LA FARM, desde el clon del repo. Despues: git add/commit/push,
# y en el Mac: git pull. Las figuras leen solo de datos_vN, nunca de la farm.
#
# Uso:
#   bash results/analisis_R/congelar_datos.sh            # crea datos_v1
#   bash results/analisis_R/congelar_datos.sh v2         # nueva version
#   bash results/analisis_R/congelar_datos.sh v1 --force # rehace v1 (evitar)
#
# Una version congelada no se toca: si cambia un dato, se crea la siguiente.

set -euo pipefail
V="${1:-v1}"
FORCE="${2:-}"

AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAIZ="$(git -C "$AQUI" rev-parse --show-toplevel)"
SCRIPTS_DIR="${SCRIPTS_DIR:-$RAIZ/scripts/sanger_scripts_jobs/last_update_05_08_2026}"
source "$SCRIPTS_DIR/config.sh" >/dev/null 2>&1

[[ "${RESULT_FROM:-}" == "MERIDA" ]] || { echo "ERROR: el lote activo es '${RESULT_FROM:-?}', no MERIDA (usa jkv)"; exit 1; }

DEST="$AQUI/datos_$V"
if [[ -e "$DEST" && "$FORCE" != "--force" ]]; then
  echo "ERROR: $DEST ya existe y esta congelado. Usa otra version (v2...) o --force."
  exit 1
fi
rm -rf "$DEST"
mkdir -p "$DEST/kraken2_nohost"

R="$RESULTS_DIR"
Q="$R/qc_control"

# destino                        origen                                                         lineas esperadas
LISTA=(
  "read_accounting.tsv           $Q/read_accounting.tsv                                          23"
  "contigs_union.tsv             $Q/contigs_union.tsv.cov                                        6005"
  "genomad_virus_summary.tsv     $R/genomad_all/all_final/all_final_summary/all_final_virus_summary.tsv 1145"
  "checkv_viral.tsv              $R/checkv_viral/all_final_virus.fna/quality_summary.tsv         1145"
  "votus.clstr                   $Q/votus.clstr                                                  -"
  "votu_coverage.tsv             $Q/votu_coverage.tsv                                            16456"
  "votus_en_genoma.tsv           $Q/votus_en_genoma.tsv                                          749"
  "TABLA_FINAL_votus_arb.tsv     $Q/TABLA_FINAL_votus_arb.tsv                                    749"
  "veredicto_contigs.tsv         $R/contig_annot/contigs_union_MERIDA_nr_jr46/veredicto_por_votu.tsv -"
)

M="$DEST/MANIFEST.tsv"
{
  echo "# congelado: $(date '+%F %T')  host: $(hostname)  lote: $RESULT_FROM"
  echo "# scripts: $(git -C "$SCRIPTS_DIR" log -1 --format='%h %s' 2>/dev/null)"
  printf 'archivo\torigen\tlineas\tmd5\tmodificado_origen\n'
} > "$M"

problemas=0
anotar() {  # anotar <destino_relativo> <origen>
  local d="$1" o="$2"
  printf '%s\t%s\t%s\t%s\t%s\n' "$d" "$o" "$(wc -l < "$DEST/$d")" \
    "$(md5sum "$DEST/$d" | cut -d' ' -f1)" "$(date -r "$o" '+%F %T')" >> "$M"
}

for fila in "${LISTA[@]}"; do
  read -r dst src esp <<< "$fila"
  if [[ ! -s "$src" ]]; then
    printf '  FALTA  %-28s %s\n' "$dst" "$src"; problemas=$((problemas + 1)); continue
  fi
  cp "$src" "$DEST/$dst"
  n=$(wc -l < "$DEST/$dst")
  flag="ok"
  [[ "$esp" != "-" && "$n" != "$esp" ]] && { flag="DIF (esperado $esp)"; problemas=$((problemas + 1)); }
  printf '  %-6s %-28s %7s lin\n' "$flag" "$dst" "$n"
  anotar "$dst" "$src"
done

nk=0
for f in "$R"/kraken2_nohost/*.bracken; do
  [[ -s "$f" ]] || continue
  cp "$f" "$DEST/kraken2_nohost/"
  anotar "kraken2_nohost/$(basename "$f")" "$f"
  nk=$((nk + 1))
done
printf '  %-6s %-28s %7s ficheros\n' "$([[ $nk == 22 ]] && echo ok || echo DIF)" "kraken2_nohost/*.bracken" "$nk"
[[ $nk == 22 ]] || problemas=$((problemas + 1))

echo
echo "  $(du -sh "$DEST" | cut -f1)  en $DEST"
if (( problemas > 0 )); then
  echo "  ATENCION: $problemas problema(s). Revisa antes de hacer commit."
  exit 1
fi
echo "  listo. Siguiente:"
echo "    git -C \"$RAIZ\" add results/analisis_R/datos_$V"
echo "    git -C \"$RAIZ\" commit -m \"datos_$V congelados para figuras\" && git -C \"$RAIZ\" push"
