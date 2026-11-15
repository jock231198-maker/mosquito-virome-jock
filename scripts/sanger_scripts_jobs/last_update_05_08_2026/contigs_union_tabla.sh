#!/bin/bash
# ---------------------------------------------------------------------------
# contigs_union_tabla.sh — genera contigs_union.tsv y contigs_union.tsv.cov
#
# POR QUE EXISTE
#   La tabla de la rama de 1 kb se construyo a mano: ningun script del repo la
#   genera, solo la consumen (tabla_final.py, diamond_contigs.sh,
#   virsorter2_union.sh, isv_sanity.sh). Este script documenta ese paso y lo
#   hace reproducible para cualquier rama.
#
# QUE HACE
#   Lee el indice .fai del FASTA de la union, cuyos nombres tienen la forma
#       M11_S3__rnaviral__NODE_1_length_12887_cov_36.291323
#   y escribe:
#
#   contigs_union.tsv      muestra  ensamblador  contig  longitud  cov
#   contigs_union.tsv.cov  ... + cov_fix
#
#   Reglas, deducidas de la rama de 1 kb:
#     muestra      = lo anterior al primer '__'
#     ensamblador  = primer campo de lo que sigue
#     contig       = ensamblador + '__' + nombre     (se conserva el prefijo)
#     longitud     = columna 2 del .fai
#     cov          = de NODE_x_length_y_cov_Z cuando el nombre lo trae
#                    (spades: rnaviral, rna, metaviral), con %.3f;
#                    -1.000 en los demas
#     cov_fix      = spades  -> copia de cov
#                    megahit -> campo multi= de la cabecera del fasta original
#                    trinity -> ausente (la fila se queda en 5 campos)
#
# USO
#   ./contigs_union_tabla.sh check                 # valida contra la rama 1 kb
#   QC_DIR=$RESULTS_DIR/qc_control_min500 \
#   FAI=$RESULTS_DIR/contig_annot_min500/contigs_union_MERIDA.fasta.fai \
#   ./contigs_union_tabla.sh run
#
# VALIDACION
#   'check' no necesita el .fai de 1 kb: reconstruye los nombres a partir del
#   propio contigs_union.tsv.cov que ya existe (muestra + '__' + contig), les
#   aplica el mismo parser y compara el resultado con el original, fila a fila.
#   Si sale IDENTICO sobre los 6,004 contigs, las reglas son correctas.
# ---------------------------------------------------------------------------
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
[[ -f "$HERE/config.sh" ]] && source "$HERE/config.sh"

: "${RESULTS_DIR:?define RESULTS_DIR (o sourcea config.sh)}"
: "${SCRATCH:?define SCRATCH (o sourcea config.sh)}"

QC_DIR="${QC_DIR:-$RESULTS_DIR/qc_control_min500}"
FAI="${FAI:-$RESULTS_DIR/contig_annot_min500/contigs_union_MERIDA.fasta.fai}"
MEGAHIT_DIR="${MEGAHIT_DIR:-$SCRATCH/megahit}"

TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT

# --- tabla multi= de megahit, por muestra y contig --------------------------
# El fasta se llama <muestra>_final.fasta; se admiten otros nombres y se elige
# el primero cuyas cabeceras lleven multi=.
build_lookup() {
  local out="$1"; : > "$out"
  local n=0 sinfa=0
  for d in "$MEGAHIT_DIR"/*/; do
    [[ -d "$d" ]] || continue
    local m fa=""
    m="$(basename "$d")"
    for cand in "${MEGAHIT_FA:-}" "${m}_final.fasta" "contigs.fasta" "final.contigs.fa"; do
      [[ -n "$cand" && -s "$d/$cand" ]] || continue
      if head -400 "$d/$cand" | grep -q '^>.*multi='; then fa="$d/$cand"; break; fi
    done
    if [[ -z "$fa" ]]; then
      sinfa=$((sinfa+1))
      echo "  aviso: en $d no hay fasta con cabeceras multi=" >&2
      continue
    fi
    grep '^>' "$fa" | awk -v s="$m" 'BEGIN{OFS="\t"}
      { id = substr($1, 2); multi = ""
        for (i = 2; i <= NF; i++) if ($i ~ /^multi=/) multi = substr($i, 7)
        if (multi != "") print s, "megahit__" id, multi }' >> "$out"
    n=$((n+1))
  done
  echo "  muestras de megahit leidas: $n   sin fasta valido: $sinfa" >&2
  echo "  contigs con multi=:         $(wc -l < "$out")" >&2
  [[ $n -gt 0 ]] || { echo "  ERROR: ninguna muestra aporto multi=" >&2; return 1; }
}

# --- parser: de (nombre_completo, longitud) a las 6 columnas ----------------
# entrada: TSV de 2 campos  nombre_completo <tab> longitud
parse_rows() {
  local lookup="$1"
  awk -F'\t' -v OFS='\t' -v LK="$lookup" '
    BEGIN {
      while ((getline line < LK) > 0) { split(line, a, "\t"); multi[a[1] SUBSEP a[2]] = a[3] }
      close(LK)
      print "muestra", "ensamblador", "contig", "longitud", "cov", "cov_fix"
    }
    {
      full = $1; len = $2
      i = index(full, "__");            if (i == 0) { bad++; next }
      muestra = substr(full, 1, i - 1)
      resto   = substr(full, i + 2)     # ensamblador__nombre
      j = index(resto, "__");           if (j == 0) { bad++; next }
      asm  = substr(resto, 1, j - 1)
      name = substr(resto, j + 2)

      # cobertura desde el nombre, estilo spades: NODE_x_length_y_cov_Z
      split(name, p, "_")
      cov = (p[3] == "length" && p[5] == "cov") ? sprintf("%.3f", p[6]) : "-1.000"

      if (asm == "rnaviral" || asm == "rna" || asm == "metaviral") {
        print muestra, asm, resto, len, cov, cov
      } else if (asm == "megahit") {
        k = muestra SUBSEP resto
        if (k in multi) print muestra, asm, resto, len, cov, multi[k]
        else { print muestra, asm, resto, len, cov, ""; sinmulti++ }
      } else {
        print muestra, asm, resto, len, cov      # trinity y otros: 5 campos
      }
    }
    END {
      if (bad)      printf("  aviso: %d nombres sin el patron muestra__ensamblador__contig\n", bad) > "/dev/stderr"
      if (sinmulti) printf("  aviso: %d contigs de megahit sin multi=\n", sinmulti) > "/dev/stderr"
    }'
}

reparto() {
  awk -F'\t' 'NR>1{
    e = (NF < 6 || $6 == "") ? "vacio" : (($5 == $6) ? "igual" : "reparado")
    c[$2"\t"e]++
  } END { for (k in c) printf("  %7d  %s\n", c[k], k) }' "$1" | sort -k2
}

case "${1:-}" in
  check)
    echo "== check: reconstruir la tabla de 1 kb y compararla con la existente =="
    REF="$RESULTS_DIR/qc_control/contigs_union.tsv.cov"
    [[ -f "$REF" ]] || { echo "falta $REF"; exit 1; }
    build_lookup "$TMP/multi.tsv"
    # nombre completo sintetico = muestra + '__' + contig ; longitud = col4
    awk -F'\t' 'NR>1 && NF>=4 {print $1 "__" $3 "\t" $4}' "$REF" > "$TMP/fai_sint.tsv"
    echo "  filas reconstruidas: $(wc -l < "$TMP/fai_sint.tsv")"
    parse_rows "$TMP/multi.tsv" < "$TMP/fai_sint.tsv" > "$TMP/repro.tsv"
    echo
    if diff -q <(sed 's/[[:space:]]*$//' "$REF") <(sed 's/[[:space:]]*$//' "$TMP/repro.tsv") >/dev/null; then
      echo "IDENTICO sobre $(( $(wc -l < "$REF") - 1 )) contigs — las reglas son correctas."
      echo "Puedes correr 'run' sobre la rama _min500."
      echo; echo "reparto en la rama de 1 kb:"; reparto "$REF"
    else
      echo "DIFIERE — no corras 'run'. Primeras diferencias:"
      diff <(sed 's/[[:space:]]*$//' "$REF") <(sed 's/[[:space:]]*$//' "$TMP/repro.tsv") | head -24
      echo; echo "lineas: original=$(wc -l < "$REF")  reproducido=$(wc -l < "$TMP/repro.tsv")"
      exit 1
    fi
    ;;

  run)
    echo "== run: generar la tabla en $QC_DIR =="
    [[ -f "$FAI" ]] || { echo "falta el indice $FAI  (samtools faidx sobre el fasta de la union)"; exit 1; }
    mkdir -p "$QC_DIR"
    for f in contigs_union.tsv contigs_union.tsv.cov; do
      [[ -f "$QC_DIR/$f" ]] && { echo "  respaldo: $f -> $f.bak"; cp "$QC_DIR/$f" "$QC_DIR/$f.bak"; }
    done
    build_lookup "$TMP/multi.tsv"
    cut -f1,2 "$FAI" | parse_rows "$TMP/multi.tsv" > "$QC_DIR/contigs_union.tsv.cov"
    # la version de 5 columnas es la misma sin cov_fix
    cut -f1-5 "$QC_DIR/contigs_union.tsv.cov" > "$QC_DIR/contigs_union.tsv"
    echo
    echo "escritos:"
    echo "  $QC_DIR/contigs_union.tsv       $(wc -l < "$QC_DIR/contigs_union.tsv") lineas"
    echo "  $QC_DIR/contigs_union.tsv.cov   $(wc -l < "$QC_DIR/contigs_union.tsv.cov") lineas"
    echo
    echo "reparto de cov_fix por ensamblador:"
    reparto "$QC_DIR/contigs_union.tsv.cov"
    ;;

  *)
    echo "uso: $0 {check|run}"
    echo "  check  reconstruye la tabla de 1 kb desde la existente y compara"
    echo "  run    genera contigs_union.tsv y .cov en QC_DIR a partir de FAI"
    echo
    echo "variables: QC_DIR  FAI  MEGAHIT_DIR  MEGAHIT_FA"
    exit 1
    ;;
esac
