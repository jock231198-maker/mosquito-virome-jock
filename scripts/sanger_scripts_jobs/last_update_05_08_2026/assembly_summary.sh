#!/bin/bash
# ---------------------------------------------------------------------------
# assembly_summary.sh — una sola tabla con TODOS los ensambladores y modos.
#
# Recorre los directorios de ensamblaje que encuentre y saca, por muestra:
#   contigs totales, contigs >=1000 bp, bases, contig mas largo
# y -- si existe el conteo de pares -- CONTIGS POR MILLON DE PARES.
#
# POR QUE LA NORMALIZACION IMPORTA
#   Comparar "numero de contigs" entre muestras sin dividir por profundidad es
#   enganoso: una muestra pequena da pocos contigs por ser pequena, no por
#   ensamblar mal. M79_S28 tiene la mejor diversidad de lecturas de las 22 y solo
#   27 contigs; con 681k pares eso es normal, no una excepcion.
#   Este era un fallo conocido del analisis anterior. Aqui queda corregido.
#
# SALIDA:
#   $RESULTS_DIR/qc_control/ensambladores_largo.tsv   una fila por muestra+modo
#   $RESULTS_DIR/qc_control/ensambladores_ancho.tsv   muestras x modos (>=1kb)
#
# Uso:
#   ./assembly_summary.sh                  # contigs_1000bp.fasta
#   ./assembly_summary.sh contigs.fasta    # todos los contigs
# ---------------------------------------------------------------------------
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

FASTA="${1:-contigs_1000bp.fasta}"
OUT="$RESULTS_DIR/qc_control"
PAIRS="$OUT/unmapped_pairs.tsv"      # sample<TAB>pares. Opcional.
mkdir -p "$OUT"

LARGO="$OUT/ensambladores_largo.tsv"
ANCHO="$OUT/ensambladores_ancho.tsv"

# Nombre de directorio -> etiqueta legible.
# El sufijo _norrna marca los ensamblajes hechos tras la depleción de rRNA, para
# que el antes y el despues convivan en la misma tabla sin confundirse.
etiqueta() {
  local d; d=$(basename "$1")
  case "$d" in
    spades)            echo "spades:rnaviral" ;;
    spades_*)          echo "spades:${d#spades_}" ;;
    megahit)           echo "megahit:default" ;;
    megahit_*)         echo "megahit:${d#megahit_}" ;;
    trinity)           echo "trinity:default" ;;
    trinity_*)         echo "trinity:${d#trinity_}" ;;
    *)                 echo "$d" ;;
  esac
}

echo "==========================================================="
echo "assembly_summary.sh   fichero: $FASTA"
echo "==========================================================="

printf 'sample\tensamblador\tcontigs_total\tcontigs_filtrados\tbases\tmas_largo\tpares\tcontigs_por_Mpar\n' > "$LARGO"

shopt -s nullglob
dirs=()
for d in "$SCRATCH"/spades "$SCRATCH"/spades_* "$SCRATCH"/megahit "$SCRATCH"/megahit_* \
         "$SCRATCH"/trinity "$SCRATCH"/trinity_*; do
  [[ -d "$d" ]] || continue
  # trinity_rescate y similares no son directorios de ensamblaje por muestra
  compgen -G "$d/*/$FASTA" >/dev/null || continue
  dirs+=( "$d" )
done

(( ${#dirs[@]} )) || { echo "No hay ningun directorio de ensamblaje con $FASTA"; exit 1; }

for d in "${dirs[@]}"; do
  lab=$(etiqueta "$d")
  for s in "$d"/*/; do
    [[ -s "$s/$FASTA" ]] || continue
    sample=$(basename "$s")

    n_filt=$(grep -c '^>' "$s/$FASTA" 2>/dev/null || echo 0)
    n_all=$(grep -c '^>' "$s/contigs.fasta" 2>/dev/null || echo "$n_filt")

    read -r bases largo < <(
      awk '/^>/{if(l>m)m=l; l=0; next}{l+=length($0); t+=length($0)}
           END{if(l>m)m=l; print t+0, m+0}' "$s/$FASTA")

    pares=""; porM=""
    if [[ -s "$PAIRS" ]]; then
      pares=$(awk -F'\t' -v s="$sample" '$1==s{print $2; exit}' "$PAIRS")
      [[ -n "$pares" ]] && (( pares > 0 )) && \
        porM=$(awk -v c="$n_filt" -v p="$pares" 'BEGIN{printf "%.1f", c/(p/1000000)}')
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$sample" "$lab" "$n_all" "$n_filt" "$bases" "$largo" "${pares:--}" "${porM:--}" \
      >> "$LARGO"
  done
done
shopt -u nullglob

# --- Tabla ancha: muestras x modos ------------------------------------------
awk -F'\t' '
  NR==1 { next }
  { v[$1"\t"$2] = $4; s[$1]; m[$2] }
  END {
    n = 0; for (k in m) mods[++n] = k
    # orden alfabetico de modos, estable entre ejecuciones
    for (i = 1; i < n; i++) for (j = i+1; j <= n; j++)
      if (mods[i] > mods[j]) { t = mods[i]; mods[i] = mods[j]; mods[j] = t }
    printf "sample"
    for (i = 1; i <= n; i++) printf "\t%s", mods[i]
    printf "\n"
    for (k in s) {
      printf "%s", k
      for (i = 1; i <= n; i++) {
        key = k "\t" mods[i]
        printf "\t%s", (key in v ? v[key] : "-")
      }
      printf "\n"
    }
  }' "$LARGO" | { read -r h; echo "$h"; sort; } > "$ANCHO"

echo
echo "-----------------------------------------------------------"
echo "CONTIGS >=1000 bp POR MUESTRA Y ENSAMBLADOR   -> $ANCHO"
echo "-----------------------------------------------------------"
column -t -s$'\t' "$ANCHO" | sed 's/^/  /'

echo
echo "-----------------------------------------------------------"
echo "RESUMEN POR ENSAMBLADOR"
echo "-----------------------------------------------------------"
awk -F'\t' 'NR>1 {
    n[$2]++; tot[$2]+=$4
    if ($6+0 > max[$2]) max[$2] = $6+0
  }
  END {
    printf "  %-22s %8s %10s %12s\n", "ensamblador", "muestras", "media>=1kb", "mas_largo"
    for (k in n) printf "  %-22s %8d %10.1f %12d\n", k, n[k], tot[k]/n[k], max[k]
  }' "$LARGO" | { read -r h; echo "$h"; sort -k2; }

if [[ ! -s "$PAIRS" ]]; then
  echo
  echo "  NOTA: falta $PAIRS, asi que no hay normalizacion por profundidad."
  echo "  Sin ella, comparar numero de contigs entre MUESTRAS es enganoso."
  echo "  Generalo con:"
  echo "    seqkit stats -T -j 8 \$SCRATCH/unmapped_fastq/*_unmapped_R1.fastq.gz \\"
  echo "      | awk -F'\\t' 'NR>1{n=split(\$1,p,\"/\"); s=p[n];"
  echo "          sub(/_unmapped_R1.fastq.gz/,\"\",s); print s\"\\t\"\$4}' \\"
  echo "      | sort -k1,1 > $PAIRS"
fi

echo
echo "Tabla larga (una fila por muestra+modo): $LARGO"