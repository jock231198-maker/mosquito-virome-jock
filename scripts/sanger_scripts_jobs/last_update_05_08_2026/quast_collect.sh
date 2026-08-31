#!/bin/bash
# ---------------------------------------------------------------------------
# quast_collect.sh — junta los 22 report.tsv de quast_compare en una sola tabla.
#
# quast_compare_worker.sh deja un informe por muestra. Esto los apila para poder
# responder "que ensamblador da mejor N50 en el conjunto del lote", que es la
# pregunta de la seccion de metodos.
#
# SALIDA:
#   $RESULTS_DIR/qc_control/quast_todas.tsv      una fila por muestra+ensamblador
#   $RESULTS_DIR/qc_control/quast_resumen.tsv    mediana por ensamblador
#
# Uso:  ./quast_collect.sh [directorio_quast_compare]
# ---------------------------------------------------------------------------
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

BASE="${1:-$RESULTS_DIR/quast_compare}"
OUT="$RESULTS_DIR/qc_control"
mkdir -p "$OUT"
TODAS="$OUT/quast_todas.tsv"
RESUMEN="$OUT/quast_resumen.tsv"

[[ -d "$BASE" ]] || { echo "No existe $BASE"; exit 1; }

# report.tsv de QUAST viene TRASPUESTO: metricas en filas, ensamblajes en
# columnas. transposed_report.tsv es el mismo dato con ensamblajes en filas, que
# es lo que hace falta para apilar. Se prefiere ese y report.tsv es el respaldo.
printf 'sample\tensamblador\tcontigs\tbases\tmas_largo\tN50\tL50\tGC_pct\n' > "$TODAS"

n=0
for d in "$BASE"/*/; do
  s=$(basename "$d")
  f="$d/transposed_report.tsv"
  [[ -s "$f" ]] || continue
  n=$((n+1))
  awk -F'\t' -v s="$s" '
    NR==1 {
      for (i = 1; i <= NF; i++) {
        h = $i
        if (h == "Assembly")                 c_asm = i
        else if (h == "# contigs")           c_ctg = i
        else if (h == "Total length")        c_len = i
        else if (h == "Largest contig")      c_max = i
        else if (h == "N50")                 c_n50 = i
        else if (h == "L50")                 c_l50 = i
        else if (h == "GC (%)")              c_gc  = i
      }
      next
    }
    {
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", s, $c_asm,
        (c_ctg ? $c_ctg : "-"), (c_len ? $c_len : "-"), (c_max ? $c_max : "-"),
        (c_n50 ? $c_n50 : "-"), (c_l50 ? $c_l50 : "-"), (c_gc ? $c_gc : "-")
    }' "$f" >> "$TODAS"
done

(( n > 0 )) || { echo "No se encontro ningun transposed_report.tsv en $BASE"; exit 1; }
echo "Muestras leidas: $n   ->  $TODAS"

# --- Mediana por ensamblador ------------------------------------------------
# Mediana y no media: con N50 unas pocas muestras muy buenas distorsionan la
# media y darian una lectura optimista de un ensamblador.
awk -F'\t' '
  NR>1 { n50[$2] = n50[$2] " " $6; ctg[$2] = ctg[$2] " " $3; cnt[$2]++ }
  function mediana(str,   a, k, i) {
    k = split(str, a, " ")
    # split deja a[1] vacio por el espacio inicial; se compacta
    for (i = 1; i <= k; i++) if (a[i] == "") { for (j = i; j < k; j++) a[j] = a[j+1]; k--; i-- }
    for (i = 1; i < k; i++) for (j = i+1; j <= k; j++)
      if (a[i]+0 > a[j]+0) { t = a[i]; a[i] = a[j]; a[j] = t }
    return (k % 2 ? a[int(k/2)+1]+0 : (a[k/2]+a[k/2+1])/2)
  }
  END {
    printf "ensamblador\tmuestras\tN50_mediana\tcontigs_mediana\n"
    for (k in cnt) printf "%s\t%d\t%.0f\t%.0f\n", k, cnt[k], mediana(n50[k]), mediana(ctg[k])
  }' "$TODAS" > "$RESUMEN"

echo
echo "-----------------------------------------------------------"
echo "RESUMEN POR ENSAMBLADOR (mediana sobre las muestras)"
echo "-----------------------------------------------------------"
{ head -1 "$RESUMEN"; tail -n +2 "$RESUMEN" | sort -t$'\t' -k3 -rn; } \
  | column -t -s$'\t' | sed 's/^/  /'

echo
echo "  N50: la mitad del ensamblaje esta en contigs de este tamano o mayores."
echo "  Mas alto = mas contiguo. Pero OJO: un ensamblador que descarta material"
echo "  dudoso sube su N50 quedandose con menos, asi que hay que leerlo JUNTO"
echo "  al numero de contigs y a las bases totales, nunca solo."
echo
echo "Tabla completa: $TODAS"