#!/bin/bash
# ---------------------------------------------------------------------------
# rrna_worker.sh — depleción de rRNA con SortMeRNA (una muestra por elemento)
#
# DONDE ENCAJA EN EL PIPELINE
#   trimmed -> bowtie2 (huesped) -> unmapped_fastq -> [ESTE PASO] -> ensamblaje
#
#   Va DESPUES de quitar el huesped, no antes: hay menos lecturas que procesar
#   (SortMeRNA es lento) y es exactamente el material que entra al ensamblador.
#
# EL DETALLE QUE LO HACE INDOLORO
#   La salida conserva el MISMO nombre de fichero que la entrada:
#       $SCRATCH/unmapped_fastq/<sample>_unmapped_R1.fastq.gz   (entrada)
#       $SCRATCH/norrna/<sample>_unmapped_R1.fastq.gz           (salida)
#   Asi assembly_spades.sh, assembly_megahit.sh y assembly_trinity.sh funcionan
#   sin tocar una linea: basta cambiarles el directorio de entrada.
#   Y los resultados actuales siguen intactos para comparar.
#
# SALIDA:
#   $SCRATCH/norrna/<sample>_unmapped_R{1,2}.fastq.gz    lecturas SIN rRNA
#   $SCRATCH/rrna/<sample>_rrna_R{1,2}.fastq.gz          el rRNA (opcional)
#   $RESULTS_DIR/rrna_stats/<sample>_rrna.tsv            cuanto se quito
#
# TRES TRAMPAS DE SORTMERNA
#
# 1. kvdb. El almacen de resultados TIENE que estar vacio antes de cada
#    ejecucion, o aborta. En la 7.x un kvdb no vacio activa la reanudacion y, si
#    las opciones no coinciden, aborta con un mensaje sobre el kvdb. Aqui cada
#    elemento usa un kvdb propio en disco local, recien creado. No hay colision.
#
# 2. El indice. Se construye en la primera ejecucion. Si lanzas 22 elementos sin
#    indice previo, los 22 lo construyen a la vez sobre el mismo directorio.
#    Por eso existe build_sortmerna_db.sh: aqui el indice se usa SOLO LECTURA.
#
# 3. paired_in vs paired_out. Con lecturas pareadas y solo una de las dos
#    alineando al rRNA:
#       --paired_in   las DOS van a 'aligned'   -> 'other' queda 100% libre de rRNA
#       --paired_out  las DOS van a 'other'     -> 'aligned' es 100% rRNA
#    Para depleción antes de ensamblar interesa la primera: preferimos perder
#    alguna lectura buena a dejar rRNA dentro. Es el defecto de este worker.
#
# Variables de entorno:
#   SMR_DB        fasta de referencia   ($REFS_DIR/sortmerna/smr_v4.3_default_db.fasta)
#   SMR_REFS      lista de fastas separada por espacios, admite varias  ($SMR_DB)
#                 ej: SMR_REFS="$REFS_DIR/sortmerna/mosquito_rrna.fasta"
#                 ej: SMR_REFS="$SMR_DB $REFS_DIR/sortmerna/mosquito_rrna.fasta"
#                 CADA conjunto de refs necesita su PROPIO SMR_IDX.
#   SMR_IDX       indice compartido     ($REFS_DIR/sortmerna/idx_default)
#   SMR_PAIRED    paired_in | paired_out                    (paired_in)
#   KEEP_RRNA     1 conserva las lecturas de rRNA           (1)
#   SMR_MEM       memoria en MB para el indexado            (4096)
#
# Prep + envio:
#   source .../config.sh
#   ls "$SCRATCH/unmapped_fastq"/*_unmapped_R1.fastq.gz | xargs -n1 basename \
#     | sed 's/_unmapped_R1.fastq.gz//' | sort > "$SCRATCH/asm_samples.txt"
#   N=$(wc -l < "$SCRATCH/asm_samples.txt")
#
#   bsub -J "rrna[1-$N]%8" \
#        -o "$LOGS_DIR/rrna.%J.%I.log" -e "$LOGS_DIR/rrna.%J.%I.err" \
#        -q long -n 8 -M 16000 \
#        -R "select[mem>16000] rusage[mem=16000] span[hosts=1]" \
#        "$SCRIPTS_DIR/rrna_worker.sh $SCRATCH/asm_samples.txt $SCRATCH/unmapped_fastq"
# ---------------------------------------------------------------------------
set -euo pipefail
SECONDS=0

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

samples="${1:?Uso: rrna_worker.sh <samples.txt> <input_dir> [outbase]}"
indir="${2:?Falta input_dir}"
outbase="${3:-$SCRATCH/norrna}"

THREADS="${LSB_DJOB_NUMPROC:-8}"
SMR_DB="${SMR_DB:-$REFS_DIR/sortmerna/smr_v4.3_default_db.fasta}"
SMR_IDX="${SMR_IDX:-$REFS_DIR/sortmerna/idx_default}"
# Varias referencias, separadas por espacios. SortMeRNA admite --ref repetido.
# OJO: cada CONJUNTO de referencias necesita su PROPIO SMR_IDX. Reutilizar el
# indice de otro conjunto da resultados silenciosamente equivocados.
SMR_REFS="${SMR_REFS:-$SMR_DB}"
SMR_PAIRED="${SMR_PAIRED:-paired_in}"
KEEP_RRNA="${KEEP_RRNA:-1}"
SMR_MEM="${SMR_MEM:-4096}"

[[ "$SMR_PAIRED" == "paired_in" || "$SMR_PAIRED" == "paired_out" ]] \
  || { echo "ERROR: SMR_PAIRED debe ser paired_in o paired_out"; exit 1; }

idx="${LSB_JOBINDEX:?Este script debe enviarse como job array}"
sample=$(sed -n "${idx}p" "$samples")
[[ -n "$sample" ]] || { echo "ERROR: linea $idx vacia en $samples"; exit 1; }

R1="$indir/${sample}_unmapped_R1.fastq.gz"
R2="$indir/${sample}_unmapped_R2.fastq.gz"
[[ -f "$R1" && -f "$R2" ]] || { echo "ERROR: faltan R1/R2 para $sample en $indir"; exit 1; }

ref_args=()
for r in $SMR_REFS; do
  [[ -s "$r" ]] || { echo "ERROR: no existe la referencia $r"; exit 1; }
  ref_args+=( --ref "$r" )
done
(( ${#ref_args[@]} )) || { echo "ERROR: SMR_REFS vacio"; exit 1; }
[[ -d "$SMR_IDX" ]] || { echo "ERROR: no existe el indice $SMR_IDX. Corre build_sortmerna_db.sh"; exit 1; }

rrnadir="${RRNA_DIR:-$SCRATCH/rrna}"
statsdir="${RRNA_STATS_DIR:-$RESULTS_DIR/rrna_stats}"
mkdir -p "$outbase" "$statsdir"
[[ "$KEEP_RRNA" == "1" ]] && mkdir -p "$rrnadir"

# --- Idempotencia -----------------------------------------------------------
O1="$outbase/${sample}_unmapped_R1.fastq.gz"
O2="$outbase/${sample}_unmapped_R2.fastq.gz"
if [[ -s "$O1" && -s "$O2" ]]; then
  echo "Ya existen las salidas de $sample. No se rehace."
  echo "Para rehacerlo:  rm -f $O1 $O2"
  exit 0
fi

# --- Trabajo en disco local -------------------------------------------------
work="${TMPDIR:-/tmp}/smr_${sample}_${LSB_JOBID:-$$}"
rm -rf "$work"; mkdir -p "$work/out"
cleanup() { rm -rf "$work"; }
trap cleanup EXIT

# kvdb propio y recien creado. Ver la trampa #1 de la cabecera.
kvdb="$work/kvdb"
rm -rf "$kvdb"

activate_env "$ENV_SORTMERNA"
command -v sortmerna >/dev/null || { echo "ERROR: sortmerna no esta en el PATH"; exit 1; }

echo "==========================================================="
echo "SortMeRNA $sample   (host: $(hostname))"
echo "  refs    : $SMR_REFS"
echo "  indice  : $SMR_IDX   (solo lectura)"
echo "  modo    : --$SMR_PAIRED"
echo "  hilos   : $THREADS"
echo "  trabajo : $work   (disco local)"
sortmerna --version 2>&1 | head -1
echo "==========================================================="

# --num_alignments 1 -> se para en el primer alineamiento que pasa el umbral.
# Es la opcion rapida y la correcta cuando solo se quiere FILTRAR, no clasificar.
sortmerna \
    "${ref_args[@]}" \
    --reads "$R1" --reads "$R2" \
    --workdir "$work" \
    --idx-dir "$SMR_IDX" \
    --kvdb "$kvdb" \
    --aligned "$work/out/rrna" \
    --other   "$work/out/norrna" \
    --fastx \
    --out2 \
    --"$SMR_PAIRED" \
    --num_alignments 1 \
    --threads "$THREADS" \
    -m "$SMR_MEM" \
    -v

# --- Localizar las salidas --------------------------------------------------
# Con --out2 escribe <base>_fwd.<ext> y <base>_rev.<ext>. La extension varia
# entre versiones (.fq / .fastq), asi que se busca en vez de darla por hecha.
find_out() {  # find_out <base> <fwd|rev>
  local b="$1" s="$2" f
  for f in "${b}_${s}".fq "${b}_${s}".fastq "${b}_${s}".fq.gz "${b}_${s}".fastq.gz; do
    [[ -s "$f" ]] && { echo "$f"; return 0; }
  done
  return 1
}

N1=$(find_out "$work/out/norrna" fwd) || {
  echo "ERROR: no aparece la salida sin rRNA. Contenido de $work/out:"; ls -la "$work/out"; exit 1; }
N2=$(find_out "$work/out/norrna" rev) || {
  echo "ERROR: falta el _rev sin rRNA."; ls -la "$work/out"; exit 1; }

# --- Publicar ---------------------------------------------------------------
gz_to() {  # gz_to <fichero> <destino.gz>
  local src="$1" dst="$2"
  if [[ "$src" == *.gz ]]; then cp "$src" "$dst.part"
  elif command -v pigz >/dev/null 2>&1; then pigz -p "$THREADS" -c "$src" > "$dst.part"
  else gzip -c "$src" > "$dst.part"; fi
  mv "$dst.part" "$dst"      # atomico: nadie ve un .gz a medias
}

gz_to "$N1" "$O1"
gz_to "$N2" "$O2"

if [[ "$KEEP_RRNA" == "1" ]]; then
  A1=$(find_out "$work/out/rrna" fwd) && gz_to "$A1" "$rrnadir/${sample}_rrna_R1.fastq.gz" || true
  A2=$(find_out "$work/out/rrna" rev) && gz_to "$A2" "$rrnadir/${sample}_rrna_R2.fastq.gz" || true
fi

# --- Estadisticas -----------------------------------------------------------
# El .log de SortMeRNA da el recuento verdadero de rRNA, que NO es el numero de
# lecturas en aligned.fasta cuando se usa paired_in. Se copia para trazabilidad.
smrlog=$(find "$work/out" -name "rrna*.log" | head -1 || true)
[[ -n "$smrlog" ]] && cp "$smrlog" "$statsdir/${sample}_sortmerna.log"

count() { echo $(( $(zcat "$1" | wc -l) / 4 )); }
in_pairs=$(count "$R1")
out_pairs=$(count "$O1")
rm_pairs=$(( in_pairs - out_pairs ))
pct=$(awk -v a="$rm_pairs" -v b="$in_pairs" 'BEGIN{printf "%.2f", (b>0? 100*a/b : 0)}')

printf 'sample\tpairs_in\tpairs_out\tpairs_rrna\tpct_rrna\tmodo\n' > "$statsdir/${sample}_rrna.tsv"
printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
  "$sample" "$in_pairs" "$out_pairs" "$rm_pairs" "$pct" "$SMR_PAIRED" \
  >> "$statsdir/${sample}_rrna.tsv"

echo
echo "  pares de entrada  : $in_pairs"
echo "  pares sin rRNA    : $out_pairs"
echo "  pares retirados   : $rm_pairs  (${pct}%)"
echo "  salida            : $O1"
echo "Done: $sample  ($((SECONDS/3600))h $((SECONDS%3600/60))m $((SECONDS%60))s)"
