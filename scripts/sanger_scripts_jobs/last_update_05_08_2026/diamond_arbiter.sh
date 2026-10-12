#!/bin/bash
# ---------------------------------------------------------------------------
# diamond_arbiter.sh — DIAMOND blastx de un FASTA contra una base .dmnd CON
#   taxonomía, arbitraje viral/celular (arbiter_nr.py), lista viral y
#   comparación con el arbitraje anterior (nr feb-2024, job 960825).
#
#   Mismos parámetros que el job 960825: --sensitive -e 1e-5 -k 25
#   -> la única variable que cambia es la base.
#
# Uso:
#   ./diamond_arbiter.sh <query.fasta> <db.dmnd> <etiqueta> [annot_por_votu.tsv|-] [outdir]
#     etiqueta : nombre corto de la base (p.ej. nr_jr46); va al nombre de la salida
#     annot    : tabla de annotate_votus.sh (base viral propia). '-' = sin arbitraje
#     outdir   : def. $RESULTS_DIR/votu_annot/<query>_<etiqueta>
#
# Variables opcionales:
#   TAXDUMP  (def. $REFS_DIR/viral_custom/20260925/taxonomy — fresco, con merged.dmp)
#   PREV     (def. $RESULTS_DIR/votu_annot/votus_nr/veredicto_por_votu.tsv)
#   SENS     (def. --sensitive)   DMND_EXTRA (p.ej. "-b 4 -c 1"; no cambia resultados)
#
# Idempotente: si diamond.tsv ya existe, no repite la búsqueda; rehace tablas.
# Envío: no es array (no usa LSB_JOBINDEX), -J sin corchetes está bien.
# ---------------------------------------------------------------------------
set -uo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

Q="${1:?uso: diamond_arbiter.sh <query.fasta> <db.dmnd> <etiqueta> [annot|-] [outdir]}"
DB="${2:?falta la base .dmnd}"
TAG="${3:?falta la etiqueta (p.ej. nr_jr46)}"
ANNOT="${4:--}"
[[ -s "$Q"  ]] || { echo "ERROR: no existe o vacío: $Q"; exit 1; }
[[ -r "$DB" ]] || { echo "ERROR: no se puede leer la base: $DB"; exit 1; }
Q="$(readlink -f "$Q")"; DB="$(readlink -f "$DB")"
qname="$(basename "${Q%.*}")"
OUT="$(readlink -m "${5:-$RESULTS_DIR/votu_annot/${qname}_${TAG}}")"
TAXDUMP="${TAXDUMP:-$REFS_DIR/viral_custom/20260925/taxonomy}"
PREV="${PREV:-$RESULTS_DIR/votu_annot/votus_nr/veredicto_por_votu.tsv}"
SENS="${SENS:---sensitive}"
DMND_EXTRA="${DMND_EXTRA:-}"
THREADS="${LSB_DJOB_NUMPROC:-8}"
ENV_VIRALDB="${ENV_VIRALDB:-viraldb}"
ARBITER="$SCRIPT_DIR/arbiter_nr.py"
TMPD="${TMPDIR:-/tmp}/dmndarb_$$"; mkdir -p "$TMPD" "$OUT"
trap 'rm -rf "$TMPD"' EXIT

# columnas EXACTAS que espera arbiter_nr.py (lista NR del script) — no tocar el orden
FMT=(qseqid sseqid pident length qlen evalue bitscore staxids sscinames qcovhsp stitle)

if [[ "$ANNOT" != "-" ]]; then
  [[ -s "$ANNOT"   ]] || { echo "ERROR: no existe la tabla annot: $ANNOT"; exit 1; }
  [[ -s "$ARBITER" ]] || { echo "ERROR: falta $ARBITER (junto a config.sh, con viral_db_tools.py)"; exit 1; }
  [[ -s "$TAXDUMP/nodes.dmp" ]] || { echo "ERROR: taxdump incompleto: $TAXDUMP"; exit 1; }
fi

echo "query   : $Q ($(grep -c '^>' "$Q") secuencias)"
echo "base    : $DB"
echo "salida  : $OUT"
echo "params  : blastx $SENS -e 1e-5 -k 25 $DMND_EXTRA | threads $THREADS | host $(hostname)"

# --- 1. la base: versión y metadatos al MANIFEST (no es nuestra: puede cambiar) --
activate_env "$ENV_DIAMOND" >/dev/null
{
  echo "fecha: $(date -Iseconds)"
  echo "query: $Q"
  echo "base: $DB"
  echo "base_ls: $(ls -l --time-style=long-iso "$DB")"
  echo "diamond: $(diamond --version 2>&1 | head -1)"
  echo "params: blastx $SENS -e 1e-5 -k 25 $DMND_EXTRA -f 6 ${FMT[*]}"
  echo "taxdump: $TAXDUMP"
  echo "annot: $ANNOT"
  echo "--- dbinfo"
} > "$OUT/MANIFEST.txt"
if ! diamond dbinfo -d "$DB" >> "$OUT/MANIFEST.txt" 2>&1; then
  tail -5 "$OUT/MANIFEST.txt"
  echo "ERROR: dbinfo falla. Lo más probable: la base se construyó con un DIAMOND más"
  echo "       nuevo que el de \$ENV_DIAMOND. Pregunta a jr46 la versión y actualiza el env."
  exit 1
fi
sed -n '/--- dbinfo/,$p' "$OUT/MANIFEST.txt"

# --- 2. DIAMOND -------------------------------------------------------------------
if [[ ! -s "$OUT/diamond.tsv" ]]; then
  echo "== DIAMOND blastx"
  # shellcheck disable=SC2086
  diamond blastx -d "$DB" -q "$Q" $SENS -e 1e-5 -k 25 $DMND_EXTRA \
      --threads "$THREADS" --tmpdir "$TMPD" \
      -f 6 "${FMT[@]}" -o "$OUT/diamond.tsv.tmp" 2> "$OUT/diamond.log"
  rc=$?
  if (( rc != 0 )); then
    tail -15 "$OUT/diamond.log"
    if grep -qiE 'taxon|taxonomy' "$OUT/diamond.log"; then
      echo ">> La base NO trae taxonomía (makedb sin --taxonmap/--taxonnodes/--taxonnames)."
      echo ">> Sin taxids no sirve para el arbitraje viral/celular. Pide a jr46 los ficheros"
      echo ">> prot.accession2taxid + taxdump con los que se construyó, o una versión con taxonomía."
    fi
    echo "ERROR: diamond salió con $rc"; exit 1
  fi
  mv "$OUT/diamond.tsv.tmp" "$OUT/diamond.tsv"
  echo "  $(cut -f1 "$OUT/diamond.tsv" | sort -u | wc -l) queries con hit | $(wc -l < "$OUT/diamond.tsv") hits"
else
  echo "== DIAMOND: ya está ($OUT/diamond.tsv) — no se repite"
fi

# chequeo de taxonomía en los hits: si casi todo es 0/vacío, el arbitraje no vale
awk -F'\t' '{n++; if($8==""||$8=="0") z++} END{
  p=(n?100*z/n:0); printf "  hits sin taxid: %d de %d (%.1f%%)\n", z, n, p
  if (p>20) print "  AVISO: muchos hits sin taxid; revisa la taxonomía de la base antes de creer el veredicto"
}' "$OUT/diamond.tsv" | tee -a "$OUT/MANIFEST.txt"

[[ "$ANNOT" == "-" ]] && { echo "Sin tabla annot: termina tras DIAMOND. Done ($((SECONDS/60))m)"; exit 0; }

# --- 3. arbitraje (mismo criterio: viral si bits_viral >= 1.1 x bits_celular) -----
activate_env "$ENV_VIRALDB" >/dev/null
echo "== arbitraje"
python3 "$ARBITER" "$TAXDUMP" "$OUT/diamond.tsv" "$ANNOT" "$OUT/veredicto_por_votu.tsv.tmp" \
  > "$OUT/resumen_veredicto.txt" || { echo "fallo arbiter_nr.py"; exit 1; }
mv "$OUT/veredicto_por_votu.tsv.tmp" "$OUT/veredicto_por_votu.tsv"
cat "$OUT/resumen_veredicto.txt"

# --- 4. lista viral (criterio de RESULTADO_votus_vs_nr.md) ------------------------
# veredicto in {viral, ambiguo}  o  (sin_hit_nr y con hit en la base viral propia)
viral_set() {
  awk -F'\t' 'NR>1 && ($18=="viral" || $18=="ambiguo" ||
              ($18=="sin_hit_nr" && $3!="" && $3!="sin_hit_viral")) {print $1}' "$1" | sort -u
}
viral_set "$OUT/veredicto_por_votu.tsv" > "$OUT/votus_virales.txt"
seqkit grep -f "$OUT/votus_virales.txt" "$Q" > "$OUT/votus_virales.fasta" 2>/dev/null
echo
echo "== lista viral: $(wc -l < "$OUT/votus_virales.txt") secuencias -> $OUT/votus_virales.fasta"

# --- 5. comparación con el arbitraje anterior --------------------------------------
if [[ -s "$PREV" && "$(readlink -f "$PREV")" != "$OUT/veredicto_por_votu.tsv" ]]; then
  echo
  echo "== cambio de veredicto: anterior ($PREV) -> esta base"
  awk -F'\t' 'FNR==1{next} NR==FNR{old[$1]=$18; next}
              {print (($1 in old)?old[$1]:"(no_estaba)") " -> " $18}' \
      "$PREV" "$OUT/veredicto_por_votu.tsv" | sort | uniq -c | sort -k1,1nr \
    | tee "$OUT/cambio_veredicto.txt"

  viral_set "$PREV" > "$TMPD/prev_viral.txt"
  comm -13 "$TMPD/prev_viral.txt" "$OUT/votus_virales.txt" > "$OUT/ganan_viral.txt"
  comm -23 "$TMPD/prev_viral.txt" "$OUT/votus_virales.txt" > "$OUT/pierden_viral.txt"
  echo "  entran en la lista viral: $(wc -l < "$OUT/ganan_viral.txt")  | salen: $(wc -l < "$OUT/pierden_viral.txt")"

  # detalle de los que cambian: votu, familia (base viral), veredicto antes/ahora, mejor hit ahora
  awk -F'\t' -v OFS='\t' '
    FILENAME==ARGV[1] {g[$1]="ENTRA"; next}
    FILENAME==ARGV[2] {g[$1]="SALE";  next}
    FILENAME==ARGV[3] {if(FNR>1) old[$1]=$18; next}
    FNR==1 {print "cambio","votu","familia_viral","antes","ahora","nr_nombre","nr_pident","nr_titulo"; next}
    ($1 in g) {print g[$1],$1,$4,old[$1],$18,$12,$7,substr($14,1,80)}' \
    "$OUT/ganan_viral.txt" "$OUT/pierden_viral.txt" "$PREV" "$OUT/veredicto_por_votu.tsv" \
    > "$OUT/cambios_lista_viral.tsv"
  echo "  detalle: $OUT/cambios_lista_viral.tsv"
else
  echo "(sin arbitraje anterior en $PREV: no hay comparación)"
fi

echo
echo "  $OUT/veredicto_por_votu.tsv  <- una fila por secuencia"
echo "  $OUT/diamond.tsv             <- todos los hits"
echo "  $OUT/MANIFEST.txt            <- base, versión, dbinfo, parámetros"
echo "Done ($((SECONDS/60))m $((SECONDS%60))s)"
