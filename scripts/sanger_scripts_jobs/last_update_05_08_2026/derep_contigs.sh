#!/bin/bash
# ---------------------------------------------------------------------------
# derep_contigs.sh — union de todos los ensambladores por muestra, sin redundancia.
#
# Junta los contigs >=1000 bp de cada ensamblaje disponible, les pone el origen
# en la cabecera, y agrupa al 95% con CD-HIT-EST para quitar lo que es la misma
# molecula ensamblada por varias vias.
#
# POR QUE
#   M71_S19 tiene el mismo contig de 15,437 bp en metaviral, megahit y rnaviral.
#   Sin dereplicar, geNomad y CheckV lo clasificarian TRES veces y la tabla final
#   contaria tres virus donde hay uno. Con dereplicacion queda un representante,
#   y la cabecera dice de que ensamblador salio.
#
#   Ademas es lo que permite aprovechar los cinco ensamblajes sin multiplicar por
#   cinco el coste de geNomad (32 GB/elemento) y CheckV (16 GB/elemento).
#
# ESTO NO SUSTITUYE A CORRER CADA ENSAMBLADOR POR SEPARADO
#   La union va a $SCRATCH/union. Los directorios originales siguen intactos, y
#   make_fasta_lists.sh genera una lista por conjunto, asi que puedes lanzar el
#   downstream sobre la union Y sobre cada ensamblador para comparar.
#
# SALIDA (mismo contrato que los workers de ensamblaje):
#   $SCRATCH/union/<sample>/contigs.fasta          <- quast
#   $SCRATCH/union/<sample>/contigs_1000bp.fasta   <- diamond
#   $SCRATCH/union/<sample>/<sample>_final.fasta   <- genomad, checkv
#   $SCRATCH/union/<sample>/<sample>.clstr         <- los grupos de CD-HIT
#   $RESULTS_DIR/derep_stats/<sample>_derep.tsv    <- cuanto aporto cada uno
#
# Variables de entorno:
#   ASM_DIRS    directorios a unir, separados por espacios
#               (por defecto: todos los que encuentre bajo $SCRATCH)
#   DEREP_ID    identidad minima para agrupar          (0.95)
#   DEREP_COV   cobertura minima del contig corto      (0.90)
#   CDHIT_MEM   memoria de CD-HIT en MB, 0 = sin tope  (4000)
#
# Prep + envio:
#   source .../config.sh
#   ls "$SCRATCH/unmapped_fastq"/*_unmapped_R1.fastq.gz | xargs -n1 basename \
#     | sed 's/_unmapped_R1.fastq.gz//' | sort > "$SCRATCH/asm_samples.txt"
#   N=$(wc -l < "$SCRATCH/asm_samples.txt")
#
#   bsub -J "derep[1-$N]" \
#        -o "$LOGS_DIR/derep.%J.%I.log" -e "$LOGS_DIR/derep.%J.%I.err" \
#        -q normal -n 2 -M 6000 \
#        -R "select[mem>6000] rusage[mem=6000] span[hosts=1]" \
#        "$SCRIPTS_DIR/derep_contigs.sh $SCRATCH/asm_samples.txt"
# ---------------------------------------------------------------------------
set -euo pipefail
SECONDS=0

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

samples="${1:?Uso: derep_contigs.sh <samples.txt> [outbase]}"
outbase="${2:-$SCRATCH/union}"

THREADS="${LSB_DJOB_NUMPROC:-2}"
DEREP_ID="${DEREP_ID:-0.95}"
DEREP_COV="${DEREP_COV:-0.90}"
CDHIT_MEM="${CDHIT_MEM:-4000}"

idx="${LSB_JOBINDEX:?Este script debe enviarse como job array}"
sample=$(sed -n "${idx}p" "$samples")
[[ -n "$sample" ]] || { echo "ERROR: linea $idx vacia en $samples"; exit 1; }

statsdir="${DEREP_STATS_DIR:-$RESULTS_DIR/derep_stats}"
final_dir="$outbase/$sample"
mkdir -p "$outbase" "$statsdir"

if [[ -s "$final_dir/contigs.fasta" ]]; then
  echo "Ya existe $final_dir/contigs.fasta. No se rehace."
  echo "Para rehacerlo:  rm -rf $final_dir"
  exit 0
fi

# --- Que ensamblajes hay ----------------------------------------------------
etiqueta() {
  local d; d=$(basename "$1")
  case "$d" in
    spades)     echo "rnaviral" ;;
    spades_*)   echo "${d#spades_}" ;;
    megahit)    echo "megahit" ;;
    megahit_*)  echo "megahit-${d#megahit_}" ;;
    trinity)    echo "trinity" ;;
    trinity_*)  echo "trinity-${d#trinity_}" ;;
    *)          echo "$d" ;;
  esac
}

if [[ -n "${ASM_DIRS:-}" ]]; then
  read -r -a dirs <<< "$ASM_DIRS"
else
  dirs=()
  shopt -s nullglob
  for d in "$SCRATCH"/spades "$SCRATCH"/spades_* "$SCRATCH"/megahit "$SCRATCH"/megahit_* \
           "$SCRATCH"/trinity "$SCRATCH"/trinity_*; do
    [[ -s "$d/$sample/contigs_1000bp.fasta" ]] && dirs+=( "$d" )
  done
  shopt -u nullglob
fi
(( ${#dirs[@]} )) || { echo "ERROR: no hay ningun ensamblaje con contigs_1000bp.fasta para $sample"; exit 1; }

work="${TMPDIR:-/tmp}/derep_${sample}_${LSB_JOBID:-$$}"
rm -rf "$work"; mkdir -p "$work"
trap 'rm -rf "$work"' EXIT

echo "==========================================================="
echo "derep $sample   (host: $(hostname))"
echo "  identidad : $DEREP_ID     cobertura: $DEREP_COV"
echo "  fuentes   : ${#dirs[@]}"
echo "==========================================================="

# --- Concatenar con el origen en la cabecera --------------------------------
# MEGAHIT mete espacios en la cabecera (>k141_1 flag=1 multi=...). Se corta en el
# primer espacio: muchas herramientas de aguas abajo lo hacen igual, y asi el
# nombre que ve geNomad es el mismo que el que aparece en el .clstr.
: > "$work/todos.fasta"
declare -A aporte
for d in "${dirs[@]}"; do
  f="$d/$sample/contigs_1000bp.fasta"
  [[ -s "$f" ]] || continue
  tag=$(etiqueta "$d")
  n=$(grep -c '^>' "$f" || true)
  aporte[$tag]=$n
  printf '  %-14s %6d contigs\n' "$tag" "$n"
  awk -v t="$tag" '/^>/ { sub(/^>/,""); split($0, a, /[ \t]/); print ">" t "__" a[1]; next } { print }' \
      "$f" >> "$work/todos.fasta"
done

n_antes=$(grep -c '^>' "$work/todos.fasta" || true)
(( n_antes > 0 )) || { echo "ERROR: nada que dereplicar"; exit 1; }

# --- Dereplicar -------------------------------------------------------------
# -c identidad, -aS cobertura del MAS CORTO (asi un fragmento contenido en un
# contig largo se absorbe), -n 10 es la longitud de palabra que toca para -c>=0.9,
# -d 0 conserva la cabecera entera en el .clstr, que es donde esta la trazabilidad.
activate_env "$ENV_CDHIT"
command -v cd-hit-est >/dev/null || { echo "ERROR: cd-hit-est no esta en el PATH"; exit 1; }
cd-hit-est --version 2>&1 | head -1 || true

cd-hit-est \
    -i "$work/todos.fasta" \
    -o "$work/derep.fasta" \
    -c "$DEREP_ID" \
    -aS "$DEREP_COV" \
    -n 10 -d 0 \
    -T "$THREADS" \
    -M "$CDHIT_MEM"

n_despues=$(grep -c '^>' "$work/derep.fasta" || true)
(( n_despues > 0 )) || { echo "ERROR: CD-HIT no dejo nada"; exit 1; }

# --- Publicar ---------------------------------------------------------------
stage="$work/publish"; mkdir -p "$stage"
cp "$work/derep.fasta"       "$stage/contigs.fasta"
cp "$work/derep.fasta"       "$stage/contigs_1000bp.fasta"   # ya venian filtrados
cp "$work/derep.fasta"       "$stage/${sample}_final.fasta"
cp "$work/derep.fasta.clstr" "$stage/${sample}.clstr"

rm -rf "$final_dir"; mkdir -p "$(dirname "$final_dir")"
mv "$stage" "$final_dir"

# --- Estadisticas -----------------------------------------------------------
# De que ensamblador salio cada representante. Si un ensamblador aporta muchos
# representantes unicos, esta encontrando cosas que los demas no.
{
  printf 'sample\torigen\trepresentantes\n'
  awk '/^>/ { sub(/^>/,""); split($0, a, "__"); print a[1] }' "$final_dir/contigs.fasta" \
    | sort | uniq -c | awk -v s="$sample" '{print s "\t" $2 "\t" $1}'
} > "$statsdir/${sample}_origen.tsv"

{
  printf 'sample\tcontigs_antes\tcontigs_despues\tredundancia_pct\tfuentes\n'
  awk -v s="$sample" -v a="$n_antes" -v b="$n_despues" -v f="${#dirs[@]}" \
    'BEGIN{printf "%s\t%d\t%d\t%.1f\t%d\n", s, a, b, 100*(a-b)/a, f}'
} > "$statsdir/${sample}_derep.tsv"

echo
echo "  contigs antes   : $n_antes"
echo "  representantes  : $n_despues"
awk -F'\t' 'NR>1{printf "    %-14s %6d representantes\n", $2, $3}' "$statsdir/${sample}_origen.tsv"
echo "  salida          : $final_dir"
echo "Done: $sample  ($((SECONDS/60))m $((SECONDS%60))s)"