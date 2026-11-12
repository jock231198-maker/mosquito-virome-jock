#!/bin/bash
# ---------------------------------------------------------------------------
# assembly_megahit.sh — ensamblaje con MEGAHIT (una muestra por elemento)
#
# Mismo contrato de salida que assembly_spades.sh, para que geNomad, CheckV,
# DIAMOND, QUAST y check_step.sh funcionen sin tocar nada:
#
#   $SCRATCH/megahit/<sample>/contigs.fasta          <- quast, read_accounting
#   $SCRATCH/megahit/<sample>/contigs_1000bp.fasta   <- diamond
#   $SCRATCH/megahit/<sample>/<sample>_final.fasta   <- genomad, checkv
#   $SCRATCH/megahit/<sample>/megahit.log, options.json
#   $RESULTS_DIR/assembly_stats_megahit/<sample>_assembly.tsv
#
# TRES MANIAS DE MEGAHIT QUE ESTE SCRIPT RESUELVE
#
# 1. Se NIEGA a arrancar si el directorio de salida ya existe. No hay bandera
#    equivalente a un "sobrescribe y calla" fiable entre versiones, asi que se
#    le da un subdirectorio que garantizamos inexistente.
#
# 2. La bandera -m, si le pasas un numero entre 0 y 1, lo interpreta como
#    FRACCION DE LA RAM TOTAL DE LA MAQUINA. En un nodo compartido de LSF eso es
#    una trampa: MEGAHIT ve los 500 GB del nodo, pide el 90%, y LSF lo mata por
#    pasarse del -M. Aqui SIEMPRE se le pasan bytes absolutos.
#
# 3. Las cabeceras NO son las de SPAdes:
#       SPAdes   >NODE_1_length_12898_cov_47.312
#       MEGAHIT  >k141_0 flag=1 multi=2.0000 len=335
#    contig_tables.sh parsea el formato de SPAdes y con esto devolveria una
#    tabla vacia. El worker escribe ademas contigs_table.tsv ya parseado.
#    'multi' es la multiplicidad del k-mer, el analogo de la cobertura.
#
# Todo el ensamblaje ocurre en el disco LOCAL del nodo y solo se publica el
# resultado a Lustre. MEGAHIT genera muchos ficheros intermedios y Lustre es
# lento para eso, ademas de que gastan inodes de la cuota.
#
# Variables de entorno:
#   MEGAHIT_PRESET   ""|meta-sensitive|meta-large          ("")
#   MEGAHIT_MEM      tope de memoria en GB (bytes absolutos) (16)
#   MEGAHIT_MINLEN   --min-contig-len                        (200)
#   MIN_LEN          longitud del filtro posterior           (1000)
#   KEEP_INTERMEDIATE 1 para conservar intermediate_contigs/ (0)
#
# OJO: MEGAHIT_MEM por debajo del -M del bsub, igual que con SPAdes.
#
# Prep + envio:
#   source .../config.sh
#   ls "$SCRATCH/unmapped_fastq"/*_unmapped_R1.fastq.gz | xargs -n1 basename \
#     | sed 's/_unmapped_R1.fastq.gz//' | sort > "$SCRATCH/asm_samples.txt"
#   N=$(wc -l < "$SCRATCH/asm_samples.txt")
#
#   bsub -J "megahit[1-$N]%6" \
#        -o "$LOGS_DIR/megahit.%J.%I.log" -e "$LOGS_DIR/megahit.%J.%I.err" \
#        -q normal -n 4 -M 20000 \
#        -R "select[mem>20000] rusage[mem=20000] span[hosts=1]" \
#        "$SCRIPTS_DIR/assembly_megahit.sh $SCRATCH/asm_samples.txt $SCRATCH/unmapped_fastq"
# ---------------------------------------------------------------------------
set -euo pipefail
SECONDS=0

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

samples="${1:?Uso: assembly_megahit.sh <samples.txt> <input_dir> [preset] [outbase]}"
indir="${2:?Falta input_dir}"
preset="${3:-${MEGAHIT_PRESET:-}}"
outbase="${4:-$SCRATCH/megahit}"

THREADS="${LSB_DJOB_NUMPROC:-4}"
MEGAHIT_MEM="${MEGAHIT_MEM:-16}"          # GB
MEGAHIT_MINLEN="${MEGAHIT_MINLEN:-200}"
MIN_LEN="${MIN_LEN:-1000}"
KEEP_INTERMEDIATE="${KEEP_INTERMEDIATE:-0}"

idx="${LSB_JOBINDEX:?Este script debe enviarse como job array}"
sample=$(sed -n "${idx}p" "$samples")
[[ -n "$sample" ]] || { echo "ERROR: linea $idx vacia en $samples"; exit 1; }

R1="$indir/${sample}_unmapped_R1.fastq.gz"
R2="$indir/${sample}_unmapped_R2.fastq.gz"
[[ -f "$R1" && -f "$R2" ]] || { echo "ERROR: faltan R1/R2 para $sample en $indir"; exit 1; }

final_dir="$outbase/$sample"
statsdir="${ASM_STATS_DIR:-$RESULTS_DIR/assembly_stats_megahit}"
mkdir -p "$outbase" "$statsdir"

# --- Idempotencia -----------------------------------------------------------
if [[ -s "$final_dir/contigs.fasta" ]]; then
  echo "Ya existe $final_dir/contigs.fasta ($(grep -c '^>' "$final_dir/contigs.fasta" || true) contigs). No se rehace."
  echo "Para rehacerlo:  rm -rf $final_dir"
  exit 0
fi

# --- Trabajo en disco local del nodo ----------------------------------------
work="${TMPDIR:-/tmp}/megahit_${sample}_${LSB_JOBID:-$$}"
rm -rf "$work"; mkdir -p "$work"
cleanup() { rm -rf "$work"; }
trap cleanup EXIT

# MEGAHIT exige que este directorio NO exista. No lo creamos.
mh_out="$work/out"

activate_env "$ENV_MEGAHIT"
command -v megahit >/dev/null || { echo "ERROR: megahit no esta en el PATH"; exit 1; }

# --- Preset -----------------------------------------------------------------
# Sin preset: k de 21 a 141 en pasos de 12. Es el defecto y va bien.
# meta-sensitive: paso de 10, mas k-mers, mas lento pero mejor con cobertura baja
#                 -- que es justo el caso de S17-S21.
# meta-large:     pensado para metagenomas de suelo, complejos. Aqui sobra.
pre=()
[[ -n "$preset" ]] && pre=( --presets "$preset" )

echo "==========================================================="
echo "MEGAHIT $sample   (host: $(hostname))"
echo "  preset  : ${preset:-<defecto: k 21-141 paso 12>}"
echo "  hilos   : $THREADS      memoria: ${MEGAHIT_MEM} GB"
echo "  entrada : $R1"
echo "  trabajo : $work   (disco local)"
megahit --version
echo "==========================================================="

# -m en BYTES absolutos. Ver la mania #2 de la cabecera.
mem_bytes=$(( MEGAHIT_MEM * 1000 * 1000 * 1000 ))

megahit \
    -1 "$R1" -2 "$R2" \
    -o "$mh_out" \
    ${pre[@]+"${pre[@]}"} \
    -t "$THREADS" \
    -m "$mem_bytes" \
    --min-contig-len "$MEGAHIT_MINLEN"

# --- Normalizar el nombre ---------------------------------------------------
[[ -s "$mh_out/final.contigs.fa" ]] || {
  echo "ERROR: MEGAHIT no dejo final.contigs.fa en $mh_out"; ls -la "$mh_out" || true; exit 1; }

stage="$work/publish"
mkdir -p "$stage"
cp "$mh_out/final.contigs.fa" "$stage/contigs.fasta"
[[ -f "$mh_out/log" ]]           && cp "$mh_out/log" "$stage/megahit.log"
[[ -f "$mh_out/options.json" ]]  && cp "$mh_out/options.json" "$stage/options.json"
if [[ "$KEEP_INTERMEDIATE" == "1" && -d "$mh_out/intermediate_contigs" ]]; then
  cp -r "$mh_out/intermediate_contigs" "$stage/"
fi

n_all=$(grep -c '^>' "$stage/contigs.fasta" || true)
(( n_all > 0 )) || { echo "ERROR: contigs.fasta vacio"; exit 1; }

# --- Filtro por longitud ----------------------------------------------------
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

cp "$filt" "$stage/${sample}_final.fasta"

# --- Tabla de contigs, formato MEGAHIT --------------------------------------
# >k141_1234 flag=1 multi=47.3120 len=12898
# Se emite con las MISMAS columnas que contig_tables.sh produce para SPAdes,
# para poder comparar ensambladores sin reescribir nada aguas abajo.
{
  printf 'sample\tnode\tlength_bp\tcoverage\n'
  awk -v s="$sample" '/^>/ {
      node=""; len=""; mul=""
      n = split(substr($0,2), f, /[ \t]+/)
      node = f[1]
      for (i = 2; i <= n; i++) {
        if (f[i] ~ /^len=/)   { len = substr(f[i], 5) }
        if (f[i] ~ /^multi=/) { mul = substr(f[i], 7) }
      }
      if (len != "") print s "\t" node "\t" len "\t" mul
    }' "$stage/contigs.fasta"
} > "$stage/contigs_table.tsv"

# --- Publicar de forma atomica ----------------------------------------------
rm -rf "$final_dir"
mkdir -p "$(dirname "$final_dir")"
mv "$stage" "$final_dir"

# --- Estadisticas -----------------------------------------------------------
# El maximo se calcula DENTRO de awk. Con 'sort -rn | head -1' y set -o pipefail,
# head cierra la tuberia, sort recibe SIGPIPE y el script muere justo despues de
# haber publicado -- que es exactamente como se perdio un TSV con SPAdes.
longest=$(awk '/^>/{if(l>m)m=l; l=0; next}{l+=length($0)}END{if(l>m)m=l; print m+0}' \
          "$final_dir/contigs.fasta")
total_bp=$(awk '!/^>/{t+=length($0)}END{print t+0}' "$final_dir/contigs.fasta")

printf 'sample\tassembler\tmode\tcontigs_total\tcontigs_ge%s\ttotal_bp\tlongest_bp\n' "$MIN_LEN" \
  > "$statsdir/${sample}_assembly.tsv"
printf '%s\tmegahit\t%s\t%s\t%s\t%s\t%s\n' \
  "$sample" "${preset:-default}" "$n_all" "$n_filt" "$total_bp" "$longest" \
  >> "$statsdir/${sample}_assembly.tsv"

echo
echo "  contigs totales   : $n_all"
echo "  contigs >=${MIN_LEN} bp : $n_filt"
echo "  bases totales     : $total_bp"
echo "  contig mas largo  : $longest bp"
echo "  salida            : $final_dir"
echo "Done: $sample  ($((SECONDS/3600))h $((SECONDS%3600/60))m $((SECONDS%60))s)"
