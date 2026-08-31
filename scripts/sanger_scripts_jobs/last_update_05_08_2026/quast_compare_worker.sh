#!/bin/bash
# ---------------------------------------------------------------------------
# quast_compare_worker.sh — QUAST comparando TODOS los ensambladores de una
#                           muestra en una sola corrida (un elemento por muestra).
#
# POR QUE ASI Y NO UNO POR ENSAMBLADOR
#   QUAST admite varios ensamblajes en la misma llamada con -l y produce un
#   informe lado a lado: N50, L50, contigs, longitud total y mas largo, todo en
#   la misma tabla y con las mismas reglas. Correrlo cinco veces por separado da
#   cinco informes que luego hay que cuadrar a mano, con el riesgo de comparar
#   numeros calculados con umbrales distintos.
#
# SIN GENOMA DE REFERENCIA
#   No se le pasa -r a proposito: no hay una referencia contra la que evaluar un
#   viroma. QUAST entra en modo sin referencia y da metricas de contiguidad
#   (N50, L50, longitudes). NO da correccion ni errores de ensamblaje, porque sin
#   referencia no se pueden calcular. Eso no es una limitacion del script.
#
# --min-contig 1000 para que las metricas correspondan al mismo conjunto que
# alimenta el pipeline (contigs_1000bp.fasta). Con el defecto de QUAST (500) los
# numeros no serian comparables con el resto de tablas del proyecto.
#
# SALIDA:
#   $RESULTS_DIR/quast_compare/<sample>/report.tsv     <- la tabla comparativa
#   $RESULTS_DIR/quast_compare/<sample>/report.html    <- para mirar
#   $RESULTS_DIR/quast_compare/<sample>/transposed_report.tsv
#
# Variables de entorno:
#   ASM_DIRS     directorios a comparar, separados por espacios
#                (por defecto: todos los que encuentre bajo $SCRATCH)
#   QUAST_MINLEN --min-contig                                  (1000)
#   QUAST_FASTA  que fichero coger de cada ensamblaje          (contigs.fasta)
#
# Prep + envio:
#   source .../config.sh
#   N=$(wc -l < "$SCRATCH/asm_samples.txt")
#   bsub -J "quastcmp[1-$N]" \
#        -o "$LOGS_DIR/quastcmp.%J.%I.log" -e "$LOGS_DIR/quastcmp.%J.%I.err" \
#        -q normal -n 2 -M 8000 \
#        -R "select[mem>8000] rusage[mem=8000] span[hosts=1]" \
#        "$SCRIPTS_DIR/quast_compare_worker.sh $SCRATCH/asm_samples.txt"
# ---------------------------------------------------------------------------
set -euo pipefail
SECONDS=0

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

samples="${1:?Uso: quast_compare_worker.sh <samples.txt> [outbase]}"
outbase="${2:-$RESULTS_DIR/quast_compare}"

THREADS="${LSB_DJOB_NUMPROC:-2}"
QUAST_MINLEN="${QUAST_MINLEN:-1000}"
QUAST_FASTA="${QUAST_FASTA:-contigs.fasta}"

idx="${LSB_JOBINDEX:?Este script debe enviarse como job array}"
sample=$(sed -n "${idx}p" "$samples")
[[ -n "$sample" ]] || { echo "ERROR: linea $idx vacia en $samples"; exit 1; }

final_dir="$outbase/$sample"
mkdir -p "$outbase"

if [[ -s "$final_dir/report.tsv" ]]; then
  echo "Ya existe $final_dir/report.tsv. No se rehace."
  echo "Para rehacerlo:  rm -rf $final_dir"
  exit 0
fi

etiqueta() {
  local d; d=$(basename "$1")
  case "$d" in
    spades)     echo "rnaviral" ;;
    spades_*)   echo "${d#spades_}" ;;
    megahit)    echo "megahit" ;;
    megahit_*)  echo "megahit-${d#megahit_}" ;;
    trinity)    echo "trinity" ;;
    trinity_*)  echo "trinity-${d#trinity_}" ;;
    union)      echo "union" ;;
    *)          echo "$d" ;;
  esac
}

# --- Reunir los ensamblajes de esta muestra ---------------------------------
if [[ -n "${ASM_DIRS:-}" ]]; then
  read -r -a dirs <<< "$ASM_DIRS"
else
  dirs=()
  shopt -s nullglob
  for d in "$SCRATCH"/spades "$SCRATCH"/spades_* "$SCRATCH"/megahit "$SCRATCH"/megahit_* \
           "$SCRATCH"/trinity "$SCRATCH"/trinity_* "$SCRATCH"/union; do
    [[ -s "$d/$sample/$QUAST_FASTA" ]] && dirs+=( "$d" )
  done
  shopt -u nullglob
fi

files=(); labels=()
for d in "${dirs[@]}"; do
  f="$d/$sample/$QUAST_FASTA"
  [[ -s "$f" ]] || continue
  files+=( "$f" )
  labels+=( "$(etiqueta "$d")" )
done
(( ${#files[@]} )) || { echo "ERROR: ningun $QUAST_FASTA para $sample"; exit 1; }

# QUAST separa las etiquetas por comas, asi que no pueden llevar comas dentro.
lab=$(IFS=,; echo "${labels[*]}")

work="${TMPDIR:-/tmp}/quastcmp_${sample}_${LSB_JOBID:-$$}"
rm -rf "$work"; mkdir -p "$work"
trap 'rm -rf "$work"' EXIT

activate_env "$ENV_QUAST"
command -v quast.py >/dev/null || command -v quast >/dev/null \
  || { echo "ERROR: quast no esta en el PATH"; exit 1; }
QUAST=$(command -v quast.py || command -v quast)

echo "==========================================================="
echo "QUAST comparativo  $sample   (host: $(hostname))"
echo "  ensamblajes : ${#files[@]}"
echo "  etiquetas   : $lab"
echo "  min-contig  : $QUAST_MINLEN"
"$QUAST" --version 2>&1 | head -1 || true
echo "==========================================================="
for i in "${!files[@]}"; do printf '  %-14s %s\n' "${labels[$i]}" "${files[$i]}"; done

"$QUAST" \
    -o "$work/out" \
    -l "$lab" \
    --min-contig "$QUAST_MINLEN" \
    --threads "$THREADS" \
    "${files[@]}"

[[ -s "$work/out/report.tsv" ]] || { echo "ERROR: QUAST no dejo report.tsv"; ls -la "$work/out" || true; exit 1; }

rm -rf "$final_dir"; mkdir -p "$(dirname "$final_dir")"
mv "$work/out" "$final_dir"

echo
echo "-----------------------------------------------------------"
column -t -s$'\t' "$final_dir/report.tsv" | sed 's/^/  /'
echo "-----------------------------------------------------------"
echo "  salida: $final_dir"
echo "Done: $sample  ($((SECONDS/60))m $((SECONDS%60))s)"