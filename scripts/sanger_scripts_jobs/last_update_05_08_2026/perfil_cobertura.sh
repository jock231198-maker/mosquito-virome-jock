#!/bin/bash
# perfil_cobertura.sh — perfiles de cobertura de las vOTUs virales y de los EVE.
#
# Dos pasos, en la farm, con el lote MERIDA activo:
#
#   bash perfil_cobertura.sh lanzar           # BED + array de LSF (22 muestras)
#   bash perfil_cobertura.sh recoger [v1]     # une las salidas y las congela en el repo
#
# "lanzar" escribe el BED con las vOTUs virales (es_virus = si) y los EVE (en el
# genoma y candidatos), para poder comparar la forma de la cobertura de unos y
# otros. "recoger" comprueba que estan las 22 muestras, las une y copia el resultado
# a results/analisis_R/datos_cobertura_vN/ con su MANIFEST. No toca datos_vN.

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh" >/dev/null 2>&1
export LC_ALL=C

[[ "${RESULT_FROM:-}" == "MERIDA" ]] || { echo "ERROR: lote activo '${RESULT_FROM:-?}', no MERIDA (usa jkv)"; exit 1; }

W=50
TABLA="$RESULTS_DIR/qc_control/TABLA_FINAL_votus_arb.tsv"
BASE="$RESULTS_DIR/qc_control/perfil_cobertura"
BED="$BASE/votus_perfil.bed"
OUT="$BASE/por_muestra"
LISTA="$SCRATCH/asm_samples.txt"

modo="${1:-}"
case "$modo" in
lanzar)
  mkdir -p "$BASE" "$OUT"
  [[ -s "$TABLA" ]] || { echo "ERROR: no existe $TABLA"; exit 1; }
  # columnas por nombre, no por posicion
  awk -F'\t' -v OFS='\t' 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next }
       $c["es_virus"] == "si" || $c["clase_final"] ~ /^EVE_/ { print $c["votu"], 0, $c["longitud"] }' \
       "$TABLA" | sort -k1,1 > "$BED"
  echo "BED: $(wc -l < "$BED") vOTUs  ($BED)"
  awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next }
       $c["es_virus"] == "si" { v++ } $c["clase_final"] ~ /^EVE_/ { e++ }
       END { printf "  virales: %d   EVE: %d\n", v, e }' "$TABLA"

  N=$(wc -l < "$LISTA")
  bsub -J "perfil[1-$N]%8" \
       -o "$LOGS_DIR/perfil.%J.%I.log" -e "$LOGS_DIR/perfil.%J.%I.err" \
       -q normal -n 4 -M 16000 \
       -R "select[mem>16000] rusage[mem=16000] span[hosts=1]" \
       "$SCRIPT_DIR/perfil_cobertura_worker.sh $LISTA $BED $OUT $W"
  echo "cuando termine (bjobs -J perfil):  bash $0 recoger v1"
  ;;

recoger)
  V="${2:-v1}"
  N=$(wc -l < "$LISTA")
  n=$(ls "$OUT"/*.perfil.tsv 2>/dev/null | wc -l)
  if (( n != N )); then
    echo "ERROR: hay $n de $N perfiles. Faltan:"
    while read -r s; do [[ -s "$OUT/$s.perfil.tsv" ]] || echo "  $s"; done < "$LISTA"
    exit 1
  fi
  todo="$BASE/perfiles_cobertura.tsv"
  { printf 'muestra\tvotu\tinicio\tfin\tprofundidad\n'; cat "$OUT"/*.perfil.tsv; } > "$todo"
  echo "unidos: $(($(wc -l < "$todo") - 1)) ventanas, $(cut -f2 "$todo" | tail -n +2 | sort -u | wc -l) vOTUs"

  RAIZ="$(git -C "$SCRIPT_DIR" rev-parse --show-toplevel)"
  DEST="$RAIZ/results/analisis_R/datos_cobertura_$V"
  [[ -e "$DEST" ]] && { echo "ERROR: $DEST ya existe (congelado). Usa otra version."; exit 1; }
  mkdir -p "$DEST"
  gzip -c "$todo" > "$DEST/perfiles_cobertura.tsv.gz"
  cp "$BED" "$DEST/votus_perfil.bed"
  {
    echo "# congelado: $(date '+%F %T')  host: $(hostname)  lote: $RESULT_FROM  ventana: $W pb"
    echo "# scripts: $(git -C "$SCRIPT_DIR" log -1 --format='%h %s' 2>/dev/null)"
    printf 'archivo\torigen\tlineas\tmd5\n'
    printf 'perfiles_cobertura.tsv.gz\t%s\t%s\t%s\n' "$todo" "$(wc -l < "$todo")" "$(md5sum "$DEST/perfiles_cobertura.tsv.gz" | cut -d' ' -f1)"
    printf 'votus_perfil.bed\t%s\t%s\t%s\n' "$BED" "$(wc -l < "$BED")" "$(md5sum "$DEST/votus_perfil.bed" | cut -d' ' -f1)"
  } > "$DEST/MANIFEST.tsv"
  echo "  $(du -sh "$DEST" | cut -f1) en $DEST"
  echo "  siguiente: git add results/analisis_R/datos_cobertura_$V && git commit -m \"perfiles de cobertura $V\" && git push"
  ;;

*)
  echo "uso: bash $0 lanzar | recoger [v1]"; exit 1 ;;
esac
