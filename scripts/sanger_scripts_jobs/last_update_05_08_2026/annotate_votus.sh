#!/bin/bash
# ---------------------------------------------------------------------------
# annotate_votus.sh — anota un FASTA (vOTUs o contigs) contra la base viral propia
#
#   blastn  (megablast + dc-megablast para los que no dan hit)  -> blast/viral_nt
#   DIAMOND blastx --very-sensitive -e 1e-3 (PLAN_diamond.md §6) -> viral_prot.dmnd
#   tabla: una fila por secuencia, mejor hit de cada uno + familia/orden/realm
#
# Uso:
#   ./annotate_votus.sh <query.fasta> [outdir]
#   outdir por defecto: $RESULTS_DIR/votu_annot/<nombre_fasta>_viraldb<version>
#
# Envío (idempotente: si relanzas, salta lo que ya está):
#   bsub -J vannot -o "$LOGS_DIR/vannot.%J.log" -e "$LOGS_DIR/vannot.%J.err" \
#        -q normal -n 8 -M 16000 -R "select[mem>16000] rusage[mem=16000] span[hosts=1]" \
#        "$PWD/annotate_votus.sh <ruta>/vOTUs.fasta"
#
# Categorías de la tabla (orientativas, NO son demarcación de especie ICTV):
#   conocido_nt         blastn >= 95% id y >= 80% de la query cubierta
#   pariente_nt         hit blastn por debajo de eso
#   divergente_solo_aa  sin hit nt, pero sí DIAMOND: candidato a virus nuevo
#   sin_hit_viral       nada contra la base viral (¿huésped?, ¿oscuro?) -> nr
# ---------------------------------------------------------------------------
set -uo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
TOOLS="$SCRIPT_DIR/viral_db_tools.py"

Q="${1:?uso: annotate_votus.sh <query.fasta> [outdir]}"
[[ -s "$Q" ]] || { echo "ERROR: no existe o vacío: $Q"; exit 1; }
Q="$(readlink -f "$Q")"
VDB="${VDB:-$REFS_DIR/viral_custom/20260925}"
VER="$(basename "$VDB")"
qname="$(basename "${Q%.*}")"
OUT="${2:-$RESULTS_DIR/votu_annot/${qname}_viraldb${VER}}"
THREADS="${LSB_DJOB_NUMPROC:-4}"
ENV_VIRALDB="${ENV_VIRALDB:-viraldb}"
TMPD="${TMPDIR:-/tmp}/vannot_$$"; mkdir -p "$TMPD" "$OUT"
trap 'rm -rf "$TMPD"' EXIT

[[ -s "$VDB/diamond/viral_prot.dmnd" && -s "$VDB/taxonomy/nodes.dmp" ]] || { echo "ERROR: base incompleta en $VDB"; exit 1; }
echo "query : $Q ($(grep -c '^>' "$Q") secuencias)"
echo "base  : $VDB"
echo "salida: $OUT"
{ echo "fecha: $(date -Iseconds)"; echo "query: $Q"; echo "base: $VDB"; } > "$OUT/MANIFEST.txt"

activate_env "$ENV_VIRALDB" >/dev/null
BLASTDB="$VDB/blast" blastdbcmd -db viral_nt -info >/dev/null 2>&1 || { echo "ERROR: no se lee $VDB/blast/viral_nt"; exit 1; }
seqkit fx2tab -n -i -l "$Q" > "$OUT/lengths.tsv"

# --- blastn -----------------------------------------------------------------
export BLASTDB="$VDB/blast"
BNFMT="6 qseqid sacc pident length qlen slen evalue bitscore staxids sscinames qcovs"
if [[ ! -s "$OUT/blastn.tsv" ]]; then
  echo "== blastn megablast"
  blastn -task megablast -db viral_nt -query "$Q" -evalue 1e-10 -max_target_seqs 50 \
         -num_threads "$THREADS" -outfmt "$BNFMT" > "$TMPD/mb.tsv" || { echo "fallo megablast"; exit 1; }
  cut -f1 "$TMPD/mb.tsv" | sort -u > "$TMPD/con_hit.txt"
  seqkit grep -v -f "$TMPD/con_hit.txt" "$Q" > "$TMPD/sin_hit.fa" 2>/dev/null
  n=$(grep -c '^>' "$TMPD/sin_hit.fa" || true)
  echo "== blastn dc-megablast sobre $n sin hit en megablast"
  : > "$TMPD/dc.tsv"
  if (( n > 0 )); then
    blastn -task dc-megablast -db viral_nt -query "$TMPD/sin_hit.fa" -evalue 1e-5 -max_target_seqs 50 \
           -num_threads "$THREADS" -outfmt "$BNFMT" > "$TMPD/dc.tsv" || { echo "fallo dc-megablast"; exit 1; }
  fi
  cat "$TMPD/mb.tsv" "$TMPD/dc.tsv" > "$OUT/blastn.tsv.tmp" && mv "$OUT/blastn.tsv.tmp" "$OUT/blastn.tsv"
  echo "  megablast: $(cut -f1 "$TMPD/mb.tsv" | sort -u | wc -l) con hit | dc-megablast: $(cut -f1 "$TMPD/dc.tsv" | sort -u | wc -l) más"
else
  echo "== blastn: ya está ($OUT/blastn.tsv)"
fi

# --- DIAMOND ----------------------------------------------------------------
if [[ ! -s "$OUT/diamond.tsv" ]]; then
  echo "== DIAMOND blastx --very-sensitive"
  activate_env "$ENV_DIAMOND" >/dev/null
  diamond blastx -d "$VDB/diamond/viral_prot.dmnd" -q "$Q" --very-sensitive -e 1e-3 -k 25 \
      --threads "$THREADS" --tmpdir "$TMPD" --quiet \
      -f 6 qseqid sseqid pident length qlen slen evalue bitscore staxids sscinames qcovhsp stitle \
      -o "$OUT/diamond.tsv.tmp" || { echo "fallo diamond"; exit 1; }
  mv "$OUT/diamond.tsv.tmp" "$OUT/diamond.tsv"
  activate_env "$ENV_VIRALDB" >/dev/null
else
  echo "== DIAMOND: ya está ($OUT/diamond.tsv)"
fi

# --- tabla ------------------------------------------------------------------
echo "== tabla por secuencia"
python3 "$TOOLS" annot --taxdump "$VDB/taxonomy" --lengths "$OUT/lengths.tsv" \
    --blastn "$OUT/blastn.tsv" --diamond "$OUT/diamond.tsv" --out "$OUT/annot_por_votu.tsv" \
  | tee "$OUT/resumen.txt"
echo
echo "  $OUT/annot_por_votu.tsv   <- una fila por vOTU"
echo "  $OUT/blastn.tsv, diamond.tsv <- todos los hits (para ver más allá del mejor)"
echo "Done ($((SECONDS/60))m $((SECONDS%60))s)"
