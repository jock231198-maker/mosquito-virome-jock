#!/bin/bash
# ---------------------------------------------------------------------------
# assembly_spades.sh — ensamblaje con SPAdes (una muestra por elemento)
#
# Adaptacion a la farm de spades_assembly.sh (version Mac). Cambios:
#   - un elemento de array por muestra, en vez de un bucle secuencial
#   - rutas de config.sh
#   - staging + idempotencia: un job matado no deja un ensamblaje a medias
#   - tmp en disco LOCAL del nodo: SPAdes hace mucha E/S pequena y Lustre sufre
#   - normaliza el nombre del ensamblaje a contigs.fasta (ver abajo)
#   - integra el filtro >=1000 bp que en la version Mac iba suelto al final
#   - borra los intermedios de SPAdes (K*/ corrected/ tmp/), que son el 80%
#
# POR QUE SE NORMALIZA EL NOMBRE
#   No todos los modos de SPAdes escriben contigs.fasta: --rna produce
#   transcripts.fasta. Pero check_step.sh, quast_worker.sh y read_accounting.sh
#   buscan "contigs.fasta" literalmente. El worker detecta lo que SPAdes haya
#   producido y lo deja tambien como contigs.fasta, asi cualquiera de los ocho
#   modos encaja con el resto del pipeline.
#
# SALIDA (la que espera todo lo demas):
#   $SCRATCH/spades/<sample>/contigs.fasta          <- quast, read_accounting
#   $SCRATCH/spades/<sample>/contigs_1000bp.fasta   <- diamond
#   $SCRATCH/spades/<sample>/<sample>_final.fasta   <- genomad, checkv
#   $SCRATCH/spades/<sample>/spades.log, params.txt
#   $RESULTS_DIR/assembly_stats/<sample>_assembly.tsv
#
# Variables de entorno:
#   SPADES_MODE    rnavirus|rna|meta|metaviral|plasmid|bio|corona|isolate  (rnavirus)
#   SPADES_MEM     tope de memoria en GB que se le pasa a SPAdes           (60)
#   MIN_LEN        longitud minima del filtro                              (1000)
#   USE_UNPAIRED   1 para meter tambien los *_unpaired del trimming        (0)
#
# OJO: el -M del bsub debe ser MAYOR que SPADES_MEM. SPAdes trata --memory como
# un tope propio; si LSF le da menos, lo mata antes de que SPAdes reaccione.
# Con SPADES_MEM=60 pide -M 68000.
#
# Prep + envio:
#   source .../config.sh
#   ls "$SCRATCH/unmapped_fastq"/*_unmapped_R1.fastq.gz | xargs -n1 basename \
#     | sed 's/_unmapped_R1.fastq.gz//' | sort > "$SCRATCH/asm_samples.txt"
#   N=$(wc -l < "$SCRATCH/asm_samples.txt")
#
#   bsub -J "spades[1-$N]%6" \
#        -o "$LOGS_DIR/spades.%J.%I.log" -e "$LOGS_DIR/spades.%J.%I.err" \
#        -q normal -n 16 -M 68000 \
#        -R "select[mem>68000] rusage[mem=68000] span[hosts=1]" \
#        "$PWD/assembly_spades.sh $SCRATCH/asm_samples.txt $SCRATCH/unmapped_fastq"
# ---------------------------------------------------------------------------
set -euo pipefail
SECONDS=0

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

samples="${1:?Uso: assembly_spades.sh <samples.txt> <input_dir> [modo] [outbase]}"
indir="${2:?Falta input_dir}"
mode="${3:-${SPADES_MODE:-rnavirus}}"
outbase="${4:-$SCRATCH/spades}"

THREADS="${LSB_DJOB_NUMPROC:-8}"
SPADES_MEM="${SPADES_MEM:-60}"      # GB
MIN_LEN="${MIN_LEN:-1000}"
USE_UNPAIRED="${USE_UNPAIRED:-0}"

# Mapa de modos (de la version Mac, tal cual)
declare -A SPADES_FLAGS=(
  [rnavirus]="--rnaviral"
  [rna]="--rna"
  [meta]="--meta"
  [metaviral]="--metaviral"
  [plasmid]="--plasmid"
  [bio]="--bio"
  [corona]="--corona"
  [isolate]="--isolate"
)
[[ -n "${SPADES_FLAGS[$mode]:-}" ]] || {
  echo "ERROR: modo '$mode' no reconocido. Disponibles: ${!SPADES_FLAGS[*]}"; exit 1; }
MODE_FLAG="${SPADES_FLAGS[$mode]}"

idx="${LSB_JOBINDEX:?Este script debe enviarse como job array}"
sample=$(sed -n "${idx}p" "$samples")
[[ -n "$sample" ]] || { echo "ERROR: linea $idx vacia en $samples"; exit 1; }

R1="$indir/${sample}_unmapped_R1.fastq.gz"
R2="$indir/${sample}_unmapped_R2.fastq.gz"
[[ -f "$R1" && -f "$R2" ]] || { echo "ERROR: faltan R1/R2 para $sample en $indir"; exit 1; }

final_dir="$outbase/$sample"
statsdir="$RESULTS_DIR/assembly_stats"
mkdir -p "$outbase" "$statsdir"

# --- Idempotencia -----------------------------------------------------------
if [[ -s "$final_dir/contigs.fasta" ]]; then
  echo "Ya existe $final_dir/contigs.fasta ($(grep -c '^>' "$final_dir/contigs.fasta") contigs). No se rehace."
  echo "Para rehacerlo:  rm -rf $final_dir"
  exit 0
fi

# --- Staging ----------------------------------------------------------------
# Fuera de $outbase: check_step.sh busca contigs.fasta con find RECURSIVO, y un
# ensamblaje a medias dentro contaminaria el recuento.
stage="${SPADES_STAGE:-$SCRATCH/.spades_stage}/$sample"
rm -rf "$stage"; mkdir -p "$stage"

# tmp en disco local del nodo. SPAdes hace mucha E/S de ficheros pequenos y en
# Lustre eso es lento y ademas molesta al resto del cluster.
tmpdir="${TMPDIR:-/tmp}/spades_${sample}_${LSB_JOBID:-$$}"
mkdir -p "$tmpdir"
cleanup() { rm -rf "$tmpdir"; }
trap cleanup EXIT

activate_env "$ENV_SPADES"
command -v spades.py >/dev/null || { echo "ERROR: spades.py no esta en el PATH"; exit 1; }

echo "==========================================================="
echo "SPAdes $sample   (host: $(hostname))"
echo "  modo    : $mode ($MODE_FLAG)"
echo "  hilos   : $THREADS      memoria: ${SPADES_MEM} GB"
echo "  entrada : $R1"
echo "  tmp     : $tmpdir"
spades.py --version
echo "==========================================================="

# --- Reads sueltas del trimming (opcional) ----------------------------------
extra=()
if [[ "$USE_UNPAIRED" == "1" ]]; then
  u1="$SCRATCH/trimmed/$RESULT_FROM/${sample}_R1_001_unpaired.fastq.gz"
  u2="$SCRATCH/trimmed/$RESULT_FROM/${sample}_R2_001_unpaired.fastq.gz"
  # OJO: estas reads NO pasaron por la depleción de huesped. Solo tienen sentido
  # si antes se mapean y filtran igual que las pareadas.
  for u in "$u1" "$u2"; do
    [[ -f "$u" ]] && extra+=( -s "$u" )
  done
  (( ${#extra[@]} )) && echo "  + ${#extra[@]} ficheros unpaired"
fi

# --- Ensamblaje -------------------------------------------------------------
spades.py \
    $MODE_FLAG \
    -1 "$R1" -2 "$R2" \
    ${extra[@]+"${extra[@]}"} \
    -o "$stage" \
    --threads "$THREADS" \
    --memory "$SPADES_MEM" \
    --tmp-dir "$tmpdir"

# --- Normalizar el nombre del ensamblaje ------------------------------------
asm=""
for cand in contigs.fasta transcripts.fasta scaffolds.fasta; do
  if [[ -s "$stage/$cand" ]]; then asm="$stage/$cand"; break; fi
done
[[ -n "$asm" ]] || { echo "ERROR: SPAdes no dejo ningun fasta de ensamblaje en $stage"; ls -la "$stage"; exit 1; }
echo "  ensamblaje encontrado: $(basename "$asm")"
[[ "$(basename "$asm")" == "contigs.fasta" ]] || cp "$asm" "$stage/contigs.fasta"

n_all=$(grep -c '^>' "$stage/contigs.fasta" || true)
(( n_all > 0 )) || { echo "ERROR: contigs.fasta vacio"; exit 1; }

# --- Filtro por longitud ----------------------------------------------------
# seqkit si esta (conserva el ajuste de linea); si no, awk, que produce una
# linea por secuencia. Las dos salidas son FASTA valido.
filt="$stage/contigs_${MIN_LEN}bp.fasta"
if command -v seqkit >/dev/null 2>&1; then
  seqkit seq -m "$MIN_LEN" "$stage/contigs.fasta" > "$filt"
else
  awk -v L="$MIN_LEN" '
    /^>/ { if (h != "" && length(s) >= L) print h ORS s; h=$0; s=""; next }
         { s = s $0 }
    END  { if (h != "" && length(s) >= L) print h ORS s }' \
    "$stage/contigs.fasta" > "$filt"
fi
n_filt=$(grep -c '^>' "$filt" || true)

# genomad_worker.sh y checkv_worker.sh buscan *_final.fasta y sacan el nombre de
# muestra del basename, asi que el fichero debe llevar el prefijo de la muestra.
cp "$filt" "$stage/${sample}_final.fasta"

# --- Limpiar intermedios ----------------------------------------------------
# K21/ K33/ K55/ corrected/ tmp/ misc/ son el grueso del directorio y solo
# sirven para --continue. Se conservan log y params para trazabilidad.
rm -rf "$stage"/{K??,K???,corrected,tmp,misc,split_input,pipeline_state} 2>/dev/null || true

# --- Publicar de forma atomica ----------------------------------------------
rm -rf "$final_dir"
mkdir -p "$(dirname "$final_dir")"
mv "$stage" "$final_dir"

# --- Estadisticas -----------------------------------------------------------
longest=$(awk '/^>/{if(l)print l; l=0; next}{l+=length($0)}END{if(l)print l}' \
          "$final_dir/contigs.fasta" | sort -rn | head -1)
total_bp=$(awk '!/^>/{t+=length($0)}END{print t+0}' "$final_dir/contigs.fasta")

printf 'sample\tmode\tcontigs_total\tcontigs_ge%s\ttotal_bp\tlongest_bp\n' "$MIN_LEN" \
  > "$statsdir/${sample}_assembly.tsv"
printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$sample" "$mode" "$n_all" "$n_filt" "$total_bp" "${longest:-0}" \
  >> "$statsdir/${sample}_assembly.tsv"

echo
echo "  contigs totales   : $n_all"
echo "  contigs >=${MIN_LEN} bp : $n_filt"
echo "  bases totales     : $total_bp"
echo "  contig mas largo  : ${longest:-0} bp"
echo "  salida            : $final_dir"
echo "Done: $sample  ($((SECONDS/3600))h $((SECONDS%3600/60))m $((SECONDS%60))s)"
