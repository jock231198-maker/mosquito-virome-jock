#!/bin/bash
# ---------------------------------------------------------------------------
# cov_fix.sh — reconstruye la columna cov_fix de contigs_union.tsv
#
# Que hace
#   Toma  $QC_DIR/contigs_union.tsv  (5 columnas: muestra ensamblador contig
#   longitud cov)  y escribe  $QC_DIR/contigs_union.tsv.cov  anadiendo una
#   sexta columna, cov_fix, segun la regla deducida de la rama de 1 kb:
#
#     spades (rnaviral|rna|metaviral) -> cov_fix = cov   (el nombre del contig
#                                                         ya trae cov_, el
#                                                         parser funciona)
#     megahit                         -> cov_fix = campo multi= de la cabecera
#                                                  de final.contigs.fa
#     trinity                         -> cov_fix vacio   (Trinity no escribe
#                                                         cobertura; tabla_final
#                                                         lo maneja con try)
#
# Uso
#   QC_DIR=$RESULTS_DIR/qc_control_min500 ./cov_fix.sh check
#   QC_DIR=$RESULTS_DIR/qc_control_min500 ./cov_fix.sh run
#
#   Por defecto QC_DIR = $RESULTS_DIR/qc_control_min500
#
# Verificacion
#   'check' reproduce la columna sobre la rama de 1 kb y la compara con el
#   contigs_union.tsv.cov que ya existe. Si no cuadra fila a fila, no corras
#   'run': la regla no es la que creemos.
# ---------------------------------------------------------------------------
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
[[ -f "$HERE/config.sh" ]] && source "$HERE/config.sh"

: "${RESULTS_DIR:?define RESULTS_DIR (o sourcea config.sh)}"
: "${SCRATCH:?define SCRATCH (o sourcea config.sh)}"

QC_DIR="${QC_DIR:-$RESULTS_DIR/qc_control_min500}"
MEGAHIT_DIR="${MEGAHIT_DIR:-$SCRATCH/megahit}"
MEGAHIT_FA="${MEGAHIT_FA:-final.contigs.fa}"

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

# --- 1. tabla de multi= por muestra y contig de megahit ---------------------
build_lookup() {
  local out="$1"
  : > "$out"
  local n=0
  for d in "$MEGAHIT_DIR"/*/; do
    [[ -d "$d" ]] || continue
    local m fa
    m="$(basename "$d")"
    fa="$d/$MEGAHIT_FA"
    [[ -f "$fa" ]] || { echo "  aviso: falta $fa" >&2; continue; }
    # cabecera megahit:  >k119_218 flag=1 multi=27.0000 len=1484
    grep '^>' "$fa" | awk -v s="$m" 'BEGIN{OFS="\t"}
      {
        id = substr($1, 2)
        multi = ""
        for (i = 2; i <= NF; i++) if ($i ~ /^multi=/) { multi = substr($i, 7) }
        if (multi != "") print s, "megahit__" id, multi
      }' >> "$out"
    n=$((n+1))
  done
  echo "  muestras de megahit leidas: $n" >&2
  echo "  contigs con multi=:         $(wc -l < "$out")" >&2
}

# --- 2. anadir cov_fix ------------------------------------------------------
add_cov_fix() {
  local union="$1" lookup="$2" out="$3"
  awk -F'\t' -v OFS='\t' -v LK="$lookup" '
    BEGIN {
      while ((getline line < LK) > 0) {
        split(line, a, "\t")
        multi[a[1] SUBSEP a[2]] = a[3]
      }
      close(LK)
    }
    NR == 1 { print $0, "cov_fix"; next }
    {
      asm = $2
      if (asm == "rnaviral" || asm == "rna" || asm == "metaviral") {
        print $0, $5
      } else if (asm == "megahit") {
        k = $1 SUBSEP $3
        if (k in multi) { print $0, multi[k] }
        else { print $0, ""; sinmulti++ }
      } else {
        # trinity y cualquier otro: columna vacia -> queda en 5 campos
        print $0
      }
    }
    END { if (sinmulti) printf("  aviso: %d contigs de megahit sin multi=\n", sinmulti) > "/dev/stderr" }
  ' "$union" > "$out"
}

case "${1:-}" in
  check)
    echo "== check: reproducir la columna sobre la rama de 1 kb =="
    REF_QC="$RESULTS_DIR/qc_control"
    [[ -f "$REF_QC/contigs_union.tsv" ]]     || { echo "falta $REF_QC/contigs_union.tsv"; exit 1; }
    [[ -f "$REF_QC/contigs_union.tsv.cov" ]] || { echo "falta $REF_QC/contigs_union.tsv.cov"; exit 1; }
    build_lookup "$TMP/multi.tsv"
    add_cov_fix "$REF_QC/contigs_union.tsv" "$TMP/multi.tsv" "$TMP/reproducido.tsv"
    echo
    if diff -q <(sed 's/[[:space:]]*$//' "$REF_QC/contigs_union.tsv.cov") \
               <(sed 's/[[:space:]]*$//' "$TMP/reproducido.tsv") >/dev/null; then
      echo "IDENTICO — la regla es correcta, puedes correr 'run'"
    else
      echo "DIFIERE — NO corras 'run'. Primeras diferencias:"
      diff <(sed 's/[[:space:]]*$//' "$REF_QC/contigs_union.tsv.cov") \
           <(sed 's/[[:space:]]*$//' "$TMP/reproducido.tsv") | head -20
      echo
      echo "conteo de lineas: existente=$(wc -l < "$REF_QC/contigs_union.tsv.cov")  reproducido=$(wc -l < "$TMP/reproducido.tsv")"
      exit 1
    fi
    ;;

  run)
    echo "== run: generar $QC_DIR/contigs_union.tsv.cov =="
    UNION="$QC_DIR/contigs_union.tsv"
    [[ -f "$UNION" ]] || { echo "falta $UNION"; exit 1; }
    if [[ -f "$QC_DIR/contigs_union.tsv.cov" ]]; then
      echo "ya existe $QC_DIR/contigs_union.tsv.cov — se respalda como .bak"
      cp "$QC_DIR/contigs_union.tsv.cov" "$QC_DIR/contigs_union.tsv.cov.bak"
    fi
    build_lookup "$TMP/multi.tsv"
    add_cov_fix "$UNION" "$TMP/multi.tsv" "$QC_DIR/contigs_union.tsv.cov"
    echo
    echo "escrito: $QC_DIR/contigs_union.tsv.cov"
    echo "lineas:  $(wc -l < "$QC_DIR/contigs_union.tsv.cov")"
    echo
    echo "reparto de cov_fix por ensamblador:"
    awk -F'\t' 'NR>1{
      estado = (NF < 6 || $6 == "") ? "vacio" : (($5 == $6) ? "igual" : "reparado")
      c[$2"\t"estado]++
    } END { for (k in c) printf("  %7d  %s\n", c[k], k) }' \
      "$QC_DIR/contigs_union.tsv.cov" | sort -k2
    ;;

  *)
    echo "uso: QC_DIR=... $0 {check|run}"
    echo "  check  reproduce la columna sobre la rama de 1 kb y la compara"
    echo "  run    genera contigs_union.tsv.cov en QC_DIR"
    exit 1
    ;;
esac
