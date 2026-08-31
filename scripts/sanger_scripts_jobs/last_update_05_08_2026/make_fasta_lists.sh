#!/bin/bash
# ---------------------------------------------------------------------------
# make_fasta_lists.sh — genera una lista de FASTA por CONJUNTO de ensamblaje.
#
# genomad_worker.sh, checkv_worker.sh, quast_worker.sh y el de diamond toman una
# LISTA como argumento 1 y un OUTBASE como argumento 2. Este script fabrica una
# lista por conjunto y te imprime el bsub correspondiente, para que puedas correr
# el downstream sobre la union Y sobre cada ensamblador por separado sin tocar
# ningun worker.
#
# SALIDA:
#   $SCRATCH/lists/fasta_<conjunto>.txt       *_final.fasta   -> genomad, checkv
#   $SCRATCH/lists/asm_<conjunto>.txt         contigs.fasta   -> quast
#   $SCRATCH/lists/filt_<conjunto>.txt        contigs_1000bp  -> diamond
#
# Uso:
#   ./make_fasta_lists.sh              # todos los conjuntos que encuentre
#   ./make_fasta_lists.sh union        # solo uno
# ---------------------------------------------------------------------------
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

solo="${1:-}"
LISTS="$SCRATCH/lists"
mkdir -p "$LISTS"

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

echo "==========================================================="
echo "make_fasta_lists.sh   ->  $LISTS"
echo "==========================================================="
printf '%-14s %8s %8s %8s   %s\n' conjunto final contigs filt directorio

shopt -s nullglob
hay=0
for d in "$SCRATCH"/union "$SCRATCH"/spades "$SCRATCH"/spades_* \
         "$SCRATCH"/megahit "$SCRATCH"/megahit_* "$SCRATCH"/trinity "$SCRATCH"/trinity_*; do
  [[ -d "$d" ]] || continue
  tag=$(etiqueta "$d")
  [[ -n "$solo" && "$tag" != "$solo" ]] && continue

  find "$d" -mindepth 2 -maxdepth 2 -name "*_final.fasta"    -size +0 | sort > "$LISTS/fasta_$tag.txt"
  find "$d" -mindepth 2 -maxdepth 2 -name "contigs.fasta"    -size +0 | sort > "$LISTS/asm_$tag.txt"
  find "$d" -mindepth 2 -maxdepth 2 -name "contigs_1000bp.fasta" -size +0 | sort > "$LISTS/filt_$tag.txt"

  n1=$(wc -l < "$LISTS/fasta_$tag.txt")
  n2=$(wc -l < "$LISTS/asm_$tag.txt")
  n3=$(wc -l < "$LISTS/filt_$tag.txt")
  (( n1 + n2 + n3 == 0 )) && { rm -f "$LISTS"/{fasta,asm,filt}_"$tag".txt; continue; }
  hay=1
  printf '%-14s %8d %8d %8d   %s\n' "$tag" "$n1" "$n2" "$n3" "$d"
done
shopt -u nullglob

(( hay )) || { echo "No se encontro ningun conjunto${solo:+ llamado '$solo'}"; exit 1; }

cat <<EOF

-----------------------------------------------------------
COMO USARLAS
-----------------------------------------------------------
Cada worker de aguas abajo toma <lista> y <outbase>. Dando un outbase distinto
por conjunto, los resultados no se pisan y se pueden comparar. Ejemplo con la
union y con rnaviral:

  C=union     # o rnaviral, metaviral, megahit, rna, meta, trinity
  N=\$(wc -l < "$LISTS/fasta_\$C.txt")

  bsub -J "genomad\$C[1-\$N]%4" \\
       -o "$LOGS_DIR/genomad\$C.%J.%I.log" -e "$LOGS_DIR/genomad\$C.%J.%I.err" \\
       -q normal -n 8 -M 32000 \\
       -R "select[mem>32000] rusage[mem=32000] span[hosts=1]" \\
       "$SCRIPTS_DIR/genomad_worker.sh $LISTS/fasta_\$C.txt \\
        $RESULTS_DIR/genomad_\$C \$GENOMAD_DB"

  bsub -J "checkv\$C[1-\$N]%6" \\
       -o "$LOGS_DIR/checkv\$C.%J.%I.log" -e "$LOGS_DIR/checkv\$C.%J.%I.err" \\
       -q normal -n 8 -M 16000 \\
       -R "select[mem>16000] rusage[mem=16000] span[hosts=1]" \\
       "$SCRIPTS_DIR/checkv_worker.sh $LISTS/fasta_\$C.txt \\
        $RESULTS_DIR/checkv_\$C \$CHECKV_DB"

  bsub -J "quast\$C[1-\$N]%10" \\
       -o "$LOGS_DIR/quast\$C.%J.%I.log" -e "$LOGS_DIR/quast\$C.%J.%I.err" \\
       -q normal -n 8 -M 8000 \\
       -R "select[mem>8000] rusage[mem=8000] span[hosts=1]" \\
       "$SCRIPTS_DIR/quast_worker.sh $LISTS/asm_\$C.txt $RESULTS_DIR/quast_\$C"

AVISO: antes de genomad/checkv, arregla el desajuste de nombres de check_step.sh
(deriva 'M1_S1' pero los workers escriben 'M1_S1_final') o contara las 22 mal.
EOF