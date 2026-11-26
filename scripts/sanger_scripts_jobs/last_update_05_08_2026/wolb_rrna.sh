#!/usr/bin/env bash
# =============================================================================
# wolb_rrna.sh -- cierra la cuestion de la Wolbachia
#
# El blast de los contigs ya dijo que no hay genoma de Wolbachia en ninguna
# muestra: amplitud maxima 0.32 %, hits en 20 de 22 muestras (no en las 4 que
# marcaba Kraken2), y la cepa del mejor hit salta entre supergrupos sin orden.
# Los bloques de bp se repiten identicos en decenas de cepas no relacionadas
# (712, 1527, 2773, 4118), que son tamanos de 16S y 23S.
#
# Esto lo prueba: corre barrnap sobre los contigs que pegaron y separa rRNA de
# no-rRNA. Si todo es rRNA, la senal de Kraken2 es un artefacto de no haber
# depletado rRNA y la cuestion queda cerrada. Si queda sequence no-rRNA
# concentrada en M20/M21, hay algo que mirar.
#
#   bash wolb_rrna.sh
#
# Necesita wolb_rrna.py en el mismo directorio y el hits.tsv que dejo
# wolbachia_contigs.sh run.
# =============================================================================
set -euo pipefail

RAMA="${RAMA:-min500}"
CPUS="${CPUS:-4}"
BARRNAP_BIN=/software/treeoflife/conda/users/envs/team222/jm79/barrnap/bin

for c in ./config.sh ../config.sh "$HOME/jkvirome/jm79/mosquito-virome-jock/scripts/sanger_scripts_jobs/last_update_05_08_2026/config.sh"; do
  if [[ -f "$c" ]]; then source "$c"; break; fi
done
: "${RESULTS_DIR:=/lustre/scratch126/tol/teams/lawniczak/users/jm79/results/MERIDA}"

if [[ "$RAMA" == "min500" ]]; then
  UNION="$RESULTS_DIR/contig_annot_min500/contigs_union_MERIDA.fasta"
else
  UNION="$RESULTS_DIR/contig_annot/contigs_union_MERIDA.fasta"
fi
OUT="$RESULTS_DIR/wolbachia_$RAMA"
AQUI="$(cd "$(dirname "$0")" && pwd)"

export PATH="$BARRNAP_BIN:$PATH"

msg() { printf '  %s\n' "$*"; }

[[ -s "$OUT/hits.tsv" ]] || { echo "ERROR: falta $OUT/hits.tsv -- corre antes wolbachia_contigs.sh run"; exit 1; }
[[ -s "$UNION" ]]        || { echo "ERROR: falta $UNION"; exit 1; }
command -v barrnap >/dev/null 2>&1 || { echo "ERROR: no encuentro barrnap en $BARRNAP_BIN"; exit 1; }

echo "== 1. extraer los contigs con hit a Wolbachia =="
python3 -I "$AQUI/wolb_rrna.py" extraer \
    "$OUT/hits.tsv" "$UNION" "$OUT/wolb_contigs.fasta"

echo
echo "== 2. barrnap: buscar genes de rRNA en esos contigs =="
barrnap --kingdom bac --threads "$CPUS" --quiet \
        "$OUT/wolb_contigs.fasta" > "$OUT/barrnap.gff" 2> "$OUT/barrnap.log" || {
  echo "barrnap fallo, mira $OUT/barrnap.log"; exit 1; }
msg "genes de rRNA anotados: $(grep -cv '^#' "$OUT/barrnap.gff" || true)"

echo
echo "== 3. rRNA frente a no-rRNA =="
python3 -I "$AQUI/wolb_rrna.py" resumir \
    "$OUT/hits.tsv" "$OUT/barrnap.gff" "$OUT/wolb" | tee "$OUT/resumen_rrna.txt"

echo
msg "escrito: $OUT/resumen_rrna.txt y $OUT/wolb_clasificacion.tsv"
