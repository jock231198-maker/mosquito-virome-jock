#!/bin/bash
# ---------------------------------------------------------------------------
# diamond_contigs.sh — los contigs de la unión (22 muestras MERIDA) contra una
#   base DIAMOND con taxonomía, con el MISMO flujo que los vOTUs:
#     annotate_votus.sh (base viral propia)  +  diamond_arbiter.sh (arbitraje)
#
#   Se corre desde el nodo de login; él mismo envía los trabajos:
#     ./diamond_contigs.sh <db.dmnd> [etiqueta=nr_jr46] [n_trozos=8]
#
#   Lo que envía:
#     1. dmnd_<tag>_c[1-N]%4  DIAMOND por trozos (cola long). Partir no cambia
#                             resultados: -k y el e-value son por query y por base.
#     2. cannot_MERIDA        annotate_votus.sh sobre el FASTA completo (normal)
#     3. arb_<tag>            espera a 1 y 2: une los trozos, arbitraje, tabla
#                             por muestra  (llama a este script en modo 'merge')
#
#   Entrada: $UNION_DIR/<muestra>/contigs.fasta (def. $SCRATCH/union), muestras de
#   $SCRATCH/asm_samples.txt. Cabecera nueva: <muestra>__<ensamblador>__<contig>
#   (sin '|': BLAST+ intenta interpretar los '|' como seq-id del NCBI).
#
#   Comprobaciones antes de enviar (si falla alguna, no envía nada):
#     - nº de contigs == ESPERADO (def. 6004, la unión reconstruida)
#     - ningún contig de spades_meta ni de trinity_rescate (excluidos a propósito)
#     - IDs únicos
#
#   Salida: $RESULTS_DIR/contig_annot/
#     contigs_union_MERIDA.fasta
#     contigs_union_MERIDA_viraldb20260925/   <- annotate_votus.sh
#     contigs_union_MERIDA_<tag>/             <- diamond.tsv, veredicto, por_muestra.tsv
# ---------------------------------------------------------------------------
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

MODE=submit
[[ "${1:-}" == "merge" ]] && { MODE=merge; shift; }
DB="${1:?uso: diamond_contigs.sh [merge] <db.dmnd> [etiqueta] [n_trozos]}"
TAG="${2:-nr_jr46}"
NCH="${3:-8}"
ESPERADO="${ESPERADO:-6004}"
UNION_DIR="${UNION_DIR:-$SCRATCH/union}"
VDB="${VDB:-$REFS_DIR/viral_custom/20260925}"
ENV_VIRALDB="${ENV_VIRALDB:-viraldb}"

BASE="$RESULTS_DIR/contig_annot"
FA="$BASE/contigs_union_MERIDA.fasta"
ANN_OUT="$BASE/contigs_union_MERIDA_viraldb$(basename "$VDB")"
OUT="$BASE/contigs_union_MERIDA_${TAG}"
CH="$OUT/chunks"
die() { echo "ERROR: $*" >&2; exit 1; }

# ===========================================================================
if [[ "$MODE" == "merge" ]]; then
  n=$(ls "$CH"/part_*.fasta 2>/dev/null | wc -l)
  m=$(ls "$CH"/out_*/diamond.tsv 2>/dev/null | wc -l)
  (( n > 0 && n == m )) || die "trozos: $n fasta, $m diamond.tsv terminados"
  [[ -s "$ANN_OUT/annot_por_votu.tsv" ]] || die "falta $ANN_OUT/annot_por_votu.tsv"
  if [[ ! -s "$OUT/diamond.tsv" ]]; then
    cat "$CH"/out_*/diamond.tsv > "$OUT/diamond.tsv.tmp" && mv "$OUT/diamond.tsv.tmp" "$OUT/diamond.tsv"
    echo "unidos $m trozos -> $OUT/diamond.tsv ($(wc -l < "$OUT/diamond.tsv") hits)"
  fi
  # PREV=- : el arbitraje anterior es de vOTUs, los IDs no son comparables
  PREV=- "$SCRIPT_DIR/diamond_arbiter.sh" "$FA" "$DB" "$TAG" "$ANN_OUT/annot_por_votu.tsv" "$OUT" \
    || die "diamond_arbiter.sh falló"

  # tabla por muestra (la muestra es lo que va antes del primer '__')
  {
    printf 'muestra\tcontigs\tviral\tambiguo\tsin_hit_nr\tcelular\n'
    awk -F'\t' -v OFS='\t' 'NR>1 { split($1, a, "__"); s=a[1]; t[s]++; v=$18
        if (v ~ /^celular_/) c[s]++; else k[s, v]++ }
      END { for (s in t) print s, t[s], k[s,"viral"]+0, k[s,"ambiguo"]+0, k[s,"sin_hit_nr"]+0, c[s]+0 }' \
      "$OUT/veredicto_por_votu.tsv" | sort -V
  } > "$OUT/por_muestra.tsv"
  echo; echo "== por muestra"; { column -t "$OUT/por_muestra.tsv" 2>/dev/null || cat "$OUT/por_muestra.tsv"; }
  echo "Done merge"
  exit 0
fi

# ===========================================================================
# submit
for s in diamond_arbiter.sh annotate_votus.sh arbiter_nr.py viral_db_tools.py; do
  [[ -s "$SCRIPT_DIR/$s" ]] || die "falta $SCRIPT_DIR/$s"
done
[[ -r "$DB" ]] || die "no se puede leer la base: $DB"
DB="$(readlink -f "$DB")"
[[ -s "$SCRATCH/asm_samples.txt" ]] || die "falta $SCRATCH/asm_samples.txt"
[[ -s "$VDB/diamond/viral_prot.dmnd" ]] || die "base viral incompleta: $VDB"
activate_env "$ENV_VIRALDB" >/dev/null
mkdir -p "$BASE" "$CH" "$LOGS_DIR"

# --- 1. FASTA único con la muestra en la cabecera ---------------------------
if [[ ! -s "$FA" ]]; then
  : > "$FA.tmp"
  while read -r s; do
    [[ -n "$s" ]] || continue
    f="$UNION_DIR/$s/contigs.fasta"
    [[ -s "$f" ]] || die "falta $f"
    sed "s/^>/>${s}__/" "$f" >> "$FA.tmp"
  done < "$SCRATCH/asm_samples.txt"
  mv "$FA.tmp" "$FA"
fi
n=$(grep -c '^>' "$FA")
echo "FASTA : $FA  ($n contigs, $(grep '^>' "$FA" | cut -d_ -f1-2 | sort -u | wc -l) muestras)"
echo "origen:"; grep '^>' "$FA" | awk -F'__' '{print "  " $2}' | sort | uniq -c

if grep -qE '^>[^ ]*__(meta|trinity-rescate)__' "$FA"; then
  rm -f "$FA"; die "hay contigs de spades_meta o trinity_rescate: $UNION_DIR no es la unión reconstruida (usa UNION_DIR=...)"
fi
if [[ "$n" -ne "$ESPERADO" ]]; then
  rm -f "$FA"; die "$n contigs y se esperaban $ESPERADO. Revisa UNION_DIR; si el número es correcto: ESPERADO=$n $0 $*"
fi
d=$(grep '^>' "$FA" | cut -d' ' -f1 | sort | uniq -d | head -3)
[[ -z "$d" ]] || { rm -f "$FA"; die "IDs repetidos: $d"; }

# --- 2. trozos --------------------------------------------------------------
if ! ls "$CH"/part_*.fasta >/dev/null 2>&1; then
  seqkit split2 -p "$NCH" -O "$CH/tmp" "$FA" 2>/dev/null
  i=0; for f in $(ls "$CH"/tmp/*.fasta | sort); do i=$((i+1)); mv "$f" "$CH/part_$i.fasta"; done
  rmdir "$CH/tmp"
fi
NCH=$(ls "$CH"/part_*.fasta | wc -l)
echo "trozos: $NCH  ($(for f in "$CH"/part_*.fasta; do grep -c '^>' "$f"; done | paste -sd' '))"

# --- 3. envíos ----------------------------------------------------------------
jid() { sed -n 's/^Job <\([0-9]*\)>.*/\1/p'; }

J1=$(bsub -J "dmnd_${TAG}_c[1-$NCH]%4" \
     -o "$LOGS_DIR/dmnd_${TAG}_c.%J.%I.log" -e "$LOGS_DIR/dmnd_${TAG}_c.%J.%I.err" \
     -q long -n 16 -M 24000 -R "select[mem>24000] rusage[mem=24000] span[hosts=1]" \
     "$SCRIPT_DIR/diamond_arbiter.sh $CH/part_\$LSB_JOBINDEX.fasta $DB $TAG - $CH/out_\$LSB_JOBINDEX" | jid)
[[ -n "$J1" ]] || die "no se envió el array de DIAMOND"

J2=$(bsub -J "cannot_MERIDA" \
     -o "$LOGS_DIR/cannot.%J.log" -e "$LOGS_DIR/cannot.%J.err" \
     -q normal -n 8 -M 16000 -R "select[mem>16000] rusage[mem=16000] span[hosts=1]" \
     "$SCRIPT_DIR/annotate_votus.sh $FA $ANN_OUT" | jid)
[[ -n "$J2" ]] || die "no se envió annotate_votus.sh"

J3=$(bsub -J "arb_${TAG}" -w "done($J1) && done($J2)" \
     -o "$LOGS_DIR/arb_${TAG}.%J.log" -e "$LOGS_DIR/arb_${TAG}.%J.err" \
     -q normal -n 2 -M 8000 -R "select[mem>8000] rusage[mem=8000] span[hosts=1]" \
     "$SCRIPT_DIR/diamond_contigs.sh merge $DB $TAG" | jid)
[[ -n "$J3" ]] || die "no se envió el merge"

echo
echo "enviados:  DIAMOND $J1 [1-$NCH]  |  base viral $J2  |  merge $J3 (espera a los dos)"
echo "seguir:    bjobs -w | grep -E '$J1|$J2|$J3'"
echo "si un trozo falla: bkill $J3; relanza solo ese índice y vuelve a enviar el merge:"
echo "  bsub ... \"$SCRIPT_DIR/diamond_contigs.sh merge $DB $TAG\""
