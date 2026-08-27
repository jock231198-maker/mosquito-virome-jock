#!/bin/bash
# ---------------------------------------------------------------------------
# contig_tables.sh — extrae longitud y cobertura de cada contig, una tabla por
#                    modo de ensamblaje.
#
# Las cabeceras de SPAdes ya traen los dos datos:
#     >NODE_1_length_12898_cov_47.312
# asi que no hace falta ninguna herramienta externa.
#
# La COBERTURA es el dato que separa una molecula abundante -un virus real, o
# rRNA- de un fragmento espurio. Ordenar por cobertura suele ser mas informativo
# que ordenar por longitud.
#
# Uso:
#   ./contig_tables.sh                 # contigs >=1000 bp (por defecto)
#   ./contig_tables.sh contigs.fasta   # todos los contigs
#
# SALIDA: $RESULTS_DIR/qc_control/contigs_<modo>.tsv   (con encabezados)
#         $RESULTS_DIR/qc_control/resumen_modos.tsv
# ---------------------------------------------------------------------------
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

FASTA="${1:-contigs_1000bp.fasta}"
OUT="$RESULTS_DIR/qc_control"
mkdir -p "$OUT"

# Nombre de modo a partir del directorio: spades -> rnaviral, spades_rna -> rna
modo_de() {
  local d; d=$(basename "$1")
  [[ "$d" == "spades" ]] && { echo rnaviral; return; }
  echo "${d#spades_}"
}

echo "==========================================================="
echo "contig_tables.sh   fichero: $FASTA"
echo "==========================================================="

RESUMEN="$OUT/resumen_modos.tsv"
printf 'modo\tmuestras\tcontigs\tlong_min\tlong_mediana\tlong_max\tcov_mediana\tcov_max\n' > "$RESUMEN"

shopt -s nullglob
for dir in "$SCRATCH"/spades "$SCRATCH"/spades_*; do
  [[ -d "$dir" ]] || continue
  modo=$(modo_de "$dir")
  tsv="$OUT/contigs_${modo}.tsv"

  printf 'sample\tnode\tlength_bp\tcoverage\n' > "$tsv"

  n_muestras=0
  for d in "$dir"/*/; do
    [[ -s "$d/$FASTA" ]] || continue
    n_muestras=$((n_muestras+1))
    s=$(basename "$d")
    # NODE_1_length_12898_cov_47.312 -> campos 2, 4 y 6 separando por '_'
    awk -v s="$s" '/^>/ {
        split(substr($0,2), a, "_")
        if (a[3]=="length" && a[5]=="cov") print s "\t" a[2] "\t" a[4] "\t" a[6]
      }' "$d/$FASTA"
  done >> "$tsv"

  n=$(( $(wc -l < "$tsv") - 1 ))
  if (( n > 0 )); then
    read -r lmin lmed lmax cmed cmax < <(
      awk -F'\t' 'NR>1 {L[++i]=$3+0; C[i]=$4+0}
        END{ n=i
             asort(L); asort(C)
             printf "%d %d %d %.1f %.1f\n", L[1], L[int((n+1)/2)], L[n], C[int((n+1)/2)], C[n] }' "$tsv" 2>/dev/null \
      || awk -F'\t' 'NR>1{print $3"\t"$4}' "$tsv" | sort -k1,1n | awk -F'\t' '
          {L[NR]=$1; C[NR]=$2} END{
            asort_l=0
            printf "%d %d %d %s %s\n", L[1], L[int((NR+1)/2)], L[NR], "-", "-"}'
    )
    printf '%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\n' \
      "$modo" "$n_muestras" "$n" "$lmin" "$lmed" "$lmax" "$cmed" "$cmax" >> "$RESUMEN"
  else
    printf '%s\t%d\t0\t-\t-\t-\t-\t-\n' "$modo" "$n_muestras" >> "$RESUMEN"
  fi

  echo
  echo "-----------------------------------------------------------"
  echo "MODO: $modo    ->  $tsv"
  echo "  muestras con salida: $n_muestras     contigs: $n"
  echo "-----------------------------------------------------------"

  (( n > 0 )) || { echo "  (sin contigs)"; continue; }

  echo
  echo "  10 mas LARGOS:"
  { head -1 "$tsv"; tail -n +2 "$tsv" | sort -t$'\t' -k3 -rn | head -10; } \
    | column -t -s$'\t' | sed 's/^/    /'

  echo
  echo "  10 de mayor COBERTURA (los mas abundantes):"
  { head -1 "$tsv"; tail -n +2 "$tsv" | sort -t$'\t' -k4 -rn | head -10; } \
    | column -t -s$'\t' | sed 's/^/    /'

  echo
  echo "  contigs por muestra:"
  tail -n +2 "$tsv" | cut -f1 | sort | uniq -c | sort -rn \
    | awk '{printf "    %-10s %6d\n", $2, $1}' | head -25
done
shopt -u nullglob

echo
echo "==========================================================="
echo "RESUMEN POR MODO   ->  $RESUMEN"
echo "==========================================================="
column -t -s$'\t' "$RESUMEN" | sed 's/^/  /'

echo
echo "Consultas utiles sobre cualquiera de las tablas:"
echo "  # los de ~12.9 kb, el contig compartido por once muestras"
echo "  awk -F'\\t' 'NR==1 || (\$3>=12500 && \$3<=13200)' $OUT/contigs_rnaviral.tsv | column -t"
echo
echo "  # los de cobertura por encima de 100"
echo "  awk -F'\\t' 'NR==1 || \$4>100' $OUT/contigs_rnaviral.tsv | sort -t\$'\\t' -k4 -rn | column -t"
