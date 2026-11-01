#!/bin/bash
# ---------------------------------------------------------------------------
# congelar_diamond.sh — congela en el repo las tablas de DIAMOND sobre los
#   6,004 contigs de la unión (MERIDA), para hacer figuras y tablas en local.
#
#   Se corre EN LA FARM, desde el clon del repo, con el lote MERIDA activo.
#   Después: git add/commit/push, y en casa: git pull + Rscript diamond_tablas.R
#
#   Uso:
#     bash results/analisis_R/diamond_R/congelar_diamond.sh check        # solo comprueba
#     bash results/analisis_R/diamond_R/congelar_diamond.sh v1           # crea datos_diamond_v1
#     bash results/analisis_R/diamond_R/congelar_diamond.sh v1 --force   # rehace v1 (evitar)
#
#   'vN' carga el taxdump (~2 GB de RAM, 1-3 min): mejor dentro de un bsub -Is.
#   Una versión congelada no se toca: si cambia un dato, se crea la siguiente.
#
#   Qué congela (todo sobre los MISMOS 6,004 contigs, >= 1,000 pb):
#     hits_viraldb.tsv.gz   base viral propia (viral_custom/20260925)  --very-sensitive -e 1e-3 -k 25
#     hits_rvdb.tsv.gz      RVDB v32.0                                  --very-sensitive, e <= 1e-3
#     hits_nr.tsv.gz        nr de jr46 (nr_tax.dmnd)                    --sensitive -e 1e-5 -k 25
#     veredicto_contigs.tsv árbitro de nr, un veredicto por contig (arbiter_nr.py)
#     annot_viraldb.tsv     anotación con la base viral (blastn + DIAMOND, familia/realm)
#     + MANIFEST de cada corrida, resumen_veredicto.txt y cruce_votus/ (si existe)
#   A nr y a la base viral se les añade grupo/realm/familia/genero por hit (diamond_taxa.py,
#   autocontenido) y se verifica que esos grupos reproducen el veredicto del árbitro.
# ---------------------------------------------------------------------------
set -euo pipefail
MODO="${1:-check}"
FORCE="${2:-}"

AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAIZ="$(git -C "$AQUI" rev-parse --show-toplevel)"
SCRIPTS_DIR="${SCRIPTS_DIR:-$RAIZ/scripts/sanger_scripts_jobs/last_update_05_08_2026}"
source "$SCRIPTS_DIR/config.sh" >/dev/null 2>&1

die() { echo "ERROR: $*" >&2; exit 1; }
[[ "${RESULT_FROM:-}" == "MERIDA" ]] || die "el lote activo es '${RESULT_FROM:-?}', no MERIDA (usa jkv)"
[[ "$MODO" == "check" || "$MODO" =~ ^v[0-9]+$ ]] || die "uso: $0 check | vN [--force]"

R="$RESULTS_DIR"
CA="$R/contig_annot"
NRD="$CA/contigs_union_MERIDA_nr_jr46"
VAD="$CA/contigs_union_MERIDA_viraldb20260925"
TAXDUMP="${TAXDUMP:-$REFS_DIR/viral_custom/20260925/taxonomy}"
ESPERADO="${ESPERADO:-6004}"

COLS_VIRALDB="qseqid sseqid pident length qlen slen qcovhsp scovhsp evalue bitscore staxids sscinames stitle"
COLS_RVDB="qseqid sseqid pident length qlen slen qcovhsp scovhsp evalue bitscore stitle"
COLS_NR="qseqid sseqid pident length qlen evalue bitscore staxids sscinames qcovhsp stitle"

# nombre   origen                                  columnas  col_evalue  e_max
HITS=(
  "viraldb  $R/diamond_TODOS_vcustom.tsv            13        9           1e-3"
  "rvdb     $R/diamond_TODOS_rvdb.tsv               11        9           1e-3"
  "nr       $NRD/diamond.tsv                        11        6           1e-5"
)
VER="$NRD/veredicto_por_votu.tsv"
ANN="$VAD/annot_por_votu.tsv"

# ----------------------------------------------------------------- comprobaciones
problemas=0
ok()  { printf '  [ok]    %s\n' "$*"; }
mal() { printf '  [MAL]   %s\n' "$*"; problemas=$((problemas + 1)); }
TMPC="$(mktemp -d)"; trap 'rm -rf "$TMPC"' EXIT

echo "== fuentes"
for f in "$VER" "$ANN" "$NRD/MANIFEST.txt" "$TAXDUMP/nodes.dmp" "$TAXDUMP/names.dmp"; do
  [[ -s "$f" ]] && ok "$f" || mal "falta o vacío: $f"
done
for fila in "${HITS[@]}"; do read -r _ src _ _ _ <<< "$fila"; [[ -s "$src" ]] && ok "$src" || mal "falta: $src"; done
(( problemas == 0 )) || die "faltan ficheros; no sigo"

echo "== contigs de referencia (veredicto del árbitro)"
tail -n +2 "$VER" | cut -f1 | sort -u > "$TMPC/ver.ids"
nv=$(wc -l < "$TMPC/ver.ids")
[[ "$nv" -eq "$ESPERADO" ]] && ok "veredicto: $nv contigs" || mal "veredicto: $nv contigs, se esperaban $ESPERADO"
[[ "$(head -1 "$VER" | awk -F'\t' '{print $1, $NF}')" == "votu veredicto" ]] \
  && ok "veredicto: cabecera votu ... veredicto" || mal "veredicto: cabecera inesperada"
tail -n +2 "$ANN" | cut -f1 | sort -u > "$TMPC/ann.ids"
cmp -s "$TMPC/ver.ids" "$TMPC/ann.ids" && ok "annot_viraldb: mismos $nv contigs" \
  || mal "annot_viraldb: los IDs no coinciden con el veredicto ($(wc -l < "$TMPC/ann.ids"))"
grep -qE '__(meta|trinity-rescate)__' "$TMPC/ver.ids" && mal "hay contigs de spades_meta o trinity_rescate" \
  || ok "sin spades_meta ni trinity_rescate"
ns=$(cut -d_ -f1-2 "$TMPC/ver.ids" | sort -u | wc -l)
[[ "$ns" -eq 22 ]] && ok "22 muestras" || mal "$ns muestras (se esperaban 22)"

echo "== tablas de hits"
for fila in "${HITS[@]}"; do
  read -r nom src ncol cev emax <<< "$fila"
  awk -F'\t' -v n="$ncol" -v c="$cev" -v e="$emax" -v nom="$nom" '
    NF != n { bad++ } $1 == "qseqid" { hdr++ }
    { if ($c + 0 > me) me = $c + 0; if (NR == 1 || $5 + 0 < ml) ml = $5 + 0 }
    END { printf "%s\t%d\t%d\t%d\t%g\t%d\t%g\n", nom, NR, bad + 0, hdr + 0, me, ml, e }' "$src" > "$TMPC/$nom.stat"
  read -r _ nl bad hdr me ml e < "$TMPC/$nom.stat"
  cut -f1 "$src" | sort -u > "$TMPC/$nom.ids"
  nq=$(wc -l < "$TMPC/$nom.ids")
  fuera=$(comm -23 "$TMPC/$nom.ids" "$TMPC/ver.ids" | wc -l)
  [[ "$bad" -eq 0 && "$hdr" -eq 0 ]] && ok "$nom: $nl hits, $ncol columnas, sin cabecera" \
    || mal "$nom: $bad líneas sin $ncol columnas, $hdr cabeceras"
  [[ "$fuera" -eq 0 ]] && ok "$nom: $nq contigs con hit, todos entre los $nv" \
    || mal "$nom: $fuera contigs que NO están en la unión (¿otra corrida, p.ej. union_min500?)"
  awk -v a="$me" -v b="$e" 'BEGIN { exit !(a <= b) }' && ok "$nom: e-value máximo $me (<= $e)" \
    || mal "$nom: e-value máximo $me > $e: no es la corrida documentada"
  [[ "$ml" -ge 1000 ]] && ok "$nom: qlen mínima $ml" || mal "$nom: qlen mínima $ml (< 1000)"
done

if [[ "$MODO" == "check" ]]; then
  echo
  (( problemas == 0 )) && echo "TODO BIEN. Para congelar: bash $0 v1   (mejor dentro de bsub -Is)" \
    || echo "$problemas problemas: no congeles hasta resolverlos"
  exit $(( problemas > 0 ))
fi
(( problemas == 0 )) || die "$problemas problemas: no se congela nada"

# ----------------------------------------------------------------- congelar
DEST="$AQUI/datos_diamond_$MODO"
if [[ -e "$DEST" && "$FORCE" != "--force" ]]; then
  die "$DEST ya existe y está congelado. Usa otra versión (v2...) o --force."
fi
TMPD="$AQUI/.datos_diamond_$MODO.tmp"; rm -rf "$TMPD"; mkdir -p "$TMPD"
trap 'rm -rf "$TMPC" "$TMPD"' EXIT

PY="${PYTHON:-python3}"
echo "== taxonomía por hit ($TAXDUMP)"
echo " viraldb"; "$PY" "$AQUI/diamond_taxa.py" --taxdump "$TAXDUMP" \
  --cols "$COLS_VIRALDB" --in "$R/diamond_TODOS_vcustom.tsv" --out "$TMPD/hits_viraldb.tsv.gz"
echo " nr";      "$PY" "$AQUI/diamond_taxa.py" --taxdump "$TAXDUMP" \
  --cols "$COLS_NR" --in "$NRD/diamond.tsv" --out "$TMPD/hits_nr.tsv.gz" --verificar "$VER"
echo " rvdb (sin taxid: solo cabecera)"
{ echo "$COLS_RVDB" | tr ' ' '\t'; cat "$R/diamond_TODOS_rvdb.tsv"; } | gzip -c > "$TMPD/hits_rvdb.tsv.gz"

cp "$VER" "$TMPD/veredicto_contigs.tsv"
cp "$ANN" "$TMPD/annot_viraldb.tsv"
cp "$NRD/MANIFEST.txt" "$TMPD/MANIFEST_nr_jr46.txt"
[[ -s "$NRD/resumen_veredicto.txt" ]] && cp "$NRD/resumen_veredicto.txt" "$TMPD/"
[[ -s "$VAD/MANIFEST.txt" ]] && cp "$VAD/MANIFEST.txt" "$TMPD/MANIFEST_viraldb.txt"
[[ -d "$NRD/cruce_votus" ]] && { mkdir -p "$TMPD/cruce_votus"; find "$NRD/cruce_votus" -maxdepth 1 -type f \
  \( -name '*.tsv' -o -name '*.txt' \) -size -20M -exec cp {} "$TMPD/cruce_votus/" \; ; }

# comprobación final: lo congelado tiene lo mismo que el origen
for nom in viraldb rvdb nr; do
  a=$(zcat "$TMPD/hits_$nom.tsv.gz" | tail -n +2 | wc -l)
  b=$(awk -F'\t' -v n="$nom" '$1 == n { print $2 }' "$TMPC/$nom.stat")
  [[ "$a" -eq "$b" ]] || die "hits_$nom: $a líneas congeladas y $b en el origen"
done

M="$TMPD/MANIFEST.tsv"
{
  echo "# congelado: $(date '+%F %T')  host: $(hostname)  lote: $RESULT_FROM  taxdump: $TAXDUMP"
  echo "# scripts: $(git -C "$SCRIPTS_DIR" log -1 --format='%h %s' 2>/dev/null || echo '?')"
  echo "# contigs: $nv (unión, >= 1000 pb, 22 muestras)"
  printf 'archivo\torigen\tlineas\tmd5\tmodificado_origen\n'
} > "$M"
fila() { printf '%s\t%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$(md5sum "$TMPD/$1" | cut -d' ' -f1)" "$(date -r "$2" '+%F %T')" >> "$M"; }
fila hits_viraldb.tsv.gz "$R/diamond_TODOS_vcustom.tsv" "$(zcat "$TMPD/hits_viraldb.tsv.gz" | wc -l)"
fila hits_rvdb.tsv.gz    "$R/diamond_TODOS_rvdb.tsv"    "$(zcat "$TMPD/hits_rvdb.tsv.gz" | wc -l)"
fila hits_nr.tsv.gz      "$NRD/diamond.tsv"             "$(zcat "$TMPD/hits_nr.tsv.gz" | wc -l)"
fila veredicto_contigs.tsv "$VER" "$(wc -l < "$TMPD/veredicto_contigs.tsv")"
fila annot_viraldb.tsv   "$ANN" "$(wc -l < "$TMPD/annot_viraldb.tsv")"

rm -rf "$DEST"; mv "$TMPD" "$DEST"
echo
echo "congelado en $DEST:"; ls -lh "$DEST" | tail -n +2 | awk '{printf "  %6s  %s\n", $5, $9}'
echo
echo "siguiente:"
echo "  git add $(realpath --relative-to="$RAIZ" "$DEST") && git commit -m 'datos_diamond_$MODO congelados' && git push"
